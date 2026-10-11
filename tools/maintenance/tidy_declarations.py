#!/usr/bin/env python3
"""Remove the local variables a routine no longer uses and order the declarations that remain.

    python3 tools/maintenance/tidy_declarations.py FILE ROUTINE            # report only
    python3 tools/maintenance/tidy_declarations.py FILE ROUTINE --apply

1. A local variable of ROUTINE that is referenced nowhere (not in the routine and not in an internal
   procedure that uses it from the routine) is removed: from its declaration, from the PRIVATE / FIRSTPRIVATE /
   REDUCTION clauses of the OpenMP directives of the routine, and from the lines that mark variables as not
   set. With `implicit none` a variable that is removed by mistake is a compilation error.

2. The declarations that remain are written in groups, each under a heading:
     - arguments, in the order of the argument list, one per line with its intent;
     - named constants (kept as they were, in their order: other declarations may depend on them);
     - switches and options read from the input   (logical and character variables read by MIO_InputParameter);
     - model parameters read from the input        (numbers read by MIO_InputParameter);
     - counters summed over the parallel loop      (REDUCTION clauses);
     - variables with one copy per thread          (PRIVATE and FIRSTPRIVATE clauses);
     - other local variables.
   Within a group: logical, integer, real, complex, character; alphabetical; one type per statement, lines
   broken below 110 columns. Declarations with attributes other than a dimension (pointer, allocatable,
   target, save), with an initial value, with a trailing comment, or inside preprocessor conditionals are
   kept as they are, after the groups.

Only the specification part of ROUTINE is rewritten; no executable statement is touched except the removal
described in 1.
"""
import argparse
import re
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import outline_block as ob          # noqa: E402
from outline_branch import DECL, split_top, decl_match    # noqa: E402

TYPE_ORDER = ["logical", "integer", "real", "complex", "character", "double"]


def word(name):
    return re.compile(r"(?<![\w%])" + re.escape(name) + r"(?!\w)", re.I)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("file")
    p.add_argument("routine")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--no-reorder", action="store_true", help="only remove the unused variables")
    a = p.parse_args()

    lines = open(a.file).read().split("\n")
    r0 = next(i for i, l in enumerate(lines) if re.match(r"^\s*(?:subroutine|function)\s+" + a.routine + r"\b", l, re.I))
    r1 = next(i for i, l in enumerate(lines) if i > r0 and re.match(r"^\s*end\s+(?:subroutine|function)\s+" + a.routine + r"\b", l, re.I))
    contains = next((i for i in range(r0, r1) if re.match(r"^\s*contains\s*$", ob.code(lines[i]), re.I)), r1)

    # ---- specification part: from the first declaration to the last one before the first executable statement
    h = r0
    while ob.code(lines[h]).rstrip().endswith("&"):       # the routine statement may be continued
        h += 1
    st = ob.statements(lines, h + 1, contains)
    spec_first = spec_last = None
    cpp = 0
    entries = []                    # (first line, last line, kind, payload)
    for i, j, text, is_omp in st:
        if is_omp:
            break
        if text.startswith("#"):
            cpp += 1 if re.match(r"#\s*if", text) else (-1 if re.match(r"#\s*endif", text) else 0)
            entries.append((i, j, "cpp", None))
            continue
        if not text:
            entries.append((i, j, "blank", None))
            continue
        if re.match(r"^(use|implicit)\b", text, re.I):
            entries.append((i, j, "head", None))
            continue
        m = decl_match(text)
        if not m and not re.match(r"^(type\s*\(|class\s*\(|external\b|intrinsic\b|parameter\b|dimension\b|save\b|data\b|interface\b|end\s+interface|procedure\b)", text, re.I):
            break
        raw = "\n".join(lines[i:j + 1])
        has_comment = any(ob.code(l).rstrip() != l.rstrip() for l in lines[i:j + 1])
        entries.append((i, j, "decl" if m else "other", (m, raw, has_comment, cpp > 0)))
        if spec_first is None:
            spec_first = i
        spec_last = j
    if spec_first is None:
        print(f"{a.routine}: no declarations")
        return
    # drop trailing non-declaration entries
    while entries and entries[-1][0] > spec_last:
        entries.pop()
    entries = [e for e in entries if e[0] >= spec_first]

    # ---- variables
    var = {}                        # lower name -> dict
    order = []
    for i, j, kind, payload in entries:
        if kind != "decl":
            continue
        m, raw, has_comment, in_cpp = payload
        typ = re.sub(r"\s+", "", m.group(1) + (m.group(2) or ""))
        attrs = split_top(m.group(3).lstrip(","))
        for item in split_top(m.group(4)):
            mm = re.match(r"^(\w+)\s*(\(.*?\))?\s*(=.*)?$", item)
            if not mm:
                sys.exit(f"line {i + 1}: declaration not understood: {item}")
            name = mm.group(1)
            is_arg = any(x.lower().startswith("intent") for x in attrs)
            var[name.lower()] = dict(name=name, typ=typ, attrs=attrs, dims=mm.group(2) or "", init=mm.group(3) or "",
                                     verbatim=(bool(attrs) and not is_arg) or bool(mm.group(3)) or has_comment or in_cpp,
                                     entry=(i, j), arg=is_arg and not in_cpp and not has_comment)
            order.append(name.lower())

    # ---- where is each variable used?
    body_lines = []
    marker = False
    for n in range(spec_last + 1, contains):
        l = lines[n]
        if "! >>> unset markers" in l:
            marker = True
        if marker or l.lstrip().lower().startswith("!$omp"):
            if "! <<< unset markers" in l:
                marker = False
            continue
        body_lines.append(ob.code(l))
    host_text = "\n".join(body_lines)
    # dimensions and initial values of other declarations count as uses
    spec_text = "\n".join(v["dims"] + " " + v["init"] + " " + " ".join(v["attrs"]) for v in var.values())
    spec_text += "\n" + "\n".join(ob.code(l) for i, j, kind, payload in entries if kind == "other" for l in lines[i:j + 1])
    # internal procedures: text of those that do not declare the name themselves
    procs = []
    k = contains + 1
    while k < r1:
        m = re.match(r"^\s*(?:subroutine|function)\s+(\w+)", lines[k], re.I)
        if m:
            e = next(i for i in range(k + 1, r1 + 1) if re.match(r"^\s*end\s+(subroutine|function)\b", lines[i], re.I))
            text = "\n".join(ob.code(l) for l in lines[k:e + 1])
            own = set()
            for i, j, t, is_omp in ob.statements(lines, k, e):
                md = decl_match(t) if not is_omp and t and not t.startswith("#") else None
                if md:
                    for item in split_top(md.group(4)):
                        mm = re.match(r"^(\w+)", item)
                        if mm:
                            own.add(mm.group(1).lower())
            # the dummy arguments are declared too, so they are in `own`
            procs.append((own, text))
            k = e + 1
        else:
            k += 1

    unused = []
    for n in order:
        if any(x.lower().startswith("intent") for x in var[n]["attrs"]):
            continue
        w = word(n)
        if w.search(host_text) or w.search(spec_text.replace(var[n]["dims"], "", 1) if False else ""):
            continue
        # used in the dimension / initial value of ANOTHER declaration?
        other_spec = "\n".join(v["dims"] + " " + v["init"] for q, v in var.items() if q != n)
        if w.search(other_spec) or w.search("\n".join(ob.code(l) for i, j, kind, payload in entries if kind == "other" for l in lines[i:j + 1])):
            continue
        if any(n not in own and w.search(text) for own, text in procs):
            continue
        unused.append(n)
    print(f"{a.routine}: {len(order)} local variables, {len(unused)} unused")

    # ---- OpenMP clauses of the routine
    reduction, private = set(), set()
    for i, j, text, is_omp in st:
        if is_omp:
            t = text.lower()
            for clause in re.findall(r"reduction\s*\(\s*[^:]*:([^)]*)\)", t):
                reduction |= {x.strip() for x in clause.split(",") if x.strip()}
            for clause in re.findall(r"(?:first)?private\s*\(([^)]*)\)", t):
                private |= {x.strip() for x in clause.split(",") if x.strip()}

    # ---- variables read from the input
    read = set()
    for own_text in [host_text] + [t for o, t in procs]:
        for m in re.finditer(r"mio_inputparameter\s*\(\s*[^,]*,\s*(\w+)", own_text, re.I):
            read.add(m.group(1).lower())

    out = list(lines)
    removed = set(unused)

    # 1a. executable part: marker lines and OpenMP clauses
    def clean_directive(first, last):
        text = " ".join(l.strip() for l in out[first:last + 1])
        joined = re.sub(r"&\s*!\$omp&?", " ", text, flags=re.I)
        joined = re.sub(r"\s*&\s*$", "", joined)
        m = re.match(r"^(!\$omp\s+[a-z ]+?)\s*((?:(?:first)?private|reduction|shared|default|schedule|collapse|num_threads|if)\s*\(.*)$",
                     joined, re.I)
        if not m:
            return None
        head, rest = m.group(1).strip(), m.group(2)
        clauses, depth, cur = [], 0, ""
        for ch in rest:
            if ch == "(":
                depth += 1
            if ch == ")":
                depth -= 1
                cur += ch
                if depth == 0:
                    clauses.append(cur.strip(" ,"))
                    cur = ""
                continue
            cur += ch
        new = []
        changed = False
        for c in clauses:
            mc = re.match(r"^(\w+)\s*\((.*)\)$", c, re.S)
            if not mc:
                return None
            key, inside = mc.group(1), mc.group(2)
            prefix = ""
            if key.lower() == "reduction":
                prefix, inside = inside.split(":", 1)
                prefix = prefix.strip() + ":"
            if key.lower() in ("private", "firstprivate", "reduction"):
                names = [x.strip() for x in inside.split(",") if x.strip()]
                kept = [x for x in names if x.lower() not in removed]
                changed |= len(kept) != len(names)
                if kept:
                    new.append((key.upper(), prefix, kept))
            else:
                new.append((key.upper(), "", [inside.strip()]))
        if not changed:
            return None
        ind = re.match(r"\s*", out[first]).group(0)
        res = []
        pieces = []
        for key, prefix, names in new:
            cur, chunk = [], f"{key}({prefix}"
            for nme in names:
                if len(chunk) + len(nme) > 95:
                    pieces.append(chunk.rstrip(", ") + ")")
                    chunk = f"{key}({prefix}"
                chunk += nme + ", "
            pieces.append(chunk.rstrip(", ") + ")")
        for q, piece in enumerate(pieces):
            start = f"{ind}{head.upper().replace('!$OMP', '!$OMP')} " if q == 0 else f"{ind}!$OMP& "
            res.append(start + piece + (", &" if q < len(pieces) - 1 else ""))
        if not pieces:
            res = [f"{ind}{head}"]
        return res

    new_body = []
    n = spec_last + 1
    marker = False
    while n < contains:
        l = out[n]
        if "! >>> unset markers" in l:
            marker = True
        if marker:
            m = re.match(r"^\s*(\w+)\s*=", l)
            if m and m.group(1).lower() in removed:
                n += 1
                continue
            if "! <<< unset markers" in l:
                marker = False
            new_body.append(l)
            n += 1
            continue
        if re.match(r"^\s*!\$omp\s+parallel", l, re.I):
            e = n
            while out[e].rstrip().endswith("&") and re.match(r"^\s*!\$omp&", out[e + 1], re.I):
                e += 1
            repl = clean_directive(n, e)
            new_body += repl if repl is not None else out[n:e + 1]
            n = e + 1
            continue
        new_body.append(l)
        n += 1

    # 1b / 2. specification part
    keep = [x for x in order if x not in removed]
    spec = []
    verbatim_entries, seen_entry = [], set()
    const_entries = []
    for i, j, kind, payload in entries:
        if kind in ("head",):
            spec += out[i:j + 1]
    groups = {"switch": [], "param": [], "count": [], "thread": [], "other": []}
    arguments = []
    for nme in keep:
        v = var[nme]
        if v.get("arg") and not a.no_reorder:
            arguments.append(nme)
            continue
        is_const = any(x.lower() == "parameter" for x in v["attrs"])
        if v["verbatim"] or a.no_reorder:
            if v["entry"] not in seen_entry:
                seen_entry.add(v["entry"])
                (const_entries if is_const else verbatim_entries).append(v["entry"])
            continue
        base = re.match(r"[a-z]+", v["typ"].lower()).group(0)
        if nme in reduction:
            groups["count"].append(nme)
        elif nme in private:
            groups["thread"].append(nme)
        elif nme in read:
            groups["switch" if base in ("logical", "character") else "param"].append(nme)
        else:
            groups["other"].append(nme)

    def rewrite_entry(entry):
        """A declaration kept as it is, without the variables that were removed."""
        i, j = entry
        names_here = [q for q in order if var[q]["entry"] == entry]
        if not any(q in removed for q in names_here):
            return out[i:j + 1]
        kept_here = [q for q in names_here if q not in removed]
        if not kept_here:
            return []
        v0 = var[kept_here[0]]
        m = decl_match(ob.statements(lines, i, j + 1)[0][2])
        left = m.group(1) + (m.group(2) or "") + m.group(3)
        comment = ""
        mc = re.search(r"(!.*)$", lines[j]) if ob.code(lines[j]).rstrip() != lines[j].rstrip() else None
        if mc:
            comment = "   " + lines[j][len(ob.code(lines[j]).rstrip()):].strip()
        ind = re.match(r"\s*", lines[i]).group(0)
        return emit(ind + left.strip() + " :: ", [var[q]["name"] + var[q]["dims"] + (" " + var[q]["init"] if var[q]["init"] else "")
                                                   for q in kept_here], comment)

    def emit(prefix, items, comment=""):
        res, cur = [], prefix
        for q, it in enumerate(items):
            piece = it + (", " if q < len(items) - 1 else "")
            if len(cur) + len(piece) > 110 and cur.strip() != prefix.strip():
                res.append(cur.rstrip() + " &")
                cur = " " * (len(prefix) - len(prefix.lstrip()) + 6)
            cur += piece
        res.append(cur.rstrip() + comment)
        return res

    ind = "   "
    titles = [("switch", "Switches and options read from the input"),
              ("param", "Model parameters read from the input"),
              ("count", "Counters summed over the parallel loop"),
              ("thread", "Variables with one copy per thread in the parallel loops"),
              ("other", "Other local variables")]
    if arguments:
        # in the order of the argument list of the routine
        headtext = " ".join(ob.code(l) for l in lines[r0:h + 1])
        pos = {q: (re.search(r"(?<![\w%])" + re.escape(q) + r"(?!\w)", headtext.lower()) or re.search("$", headtext)).start()
               for q in arguments}
        spec += ["", ind + "! Arguments"]
        for q in sorted(arguments, key=lambda x: pos[x]):
            v = var[q]
            m = decl_match(ob.statements(lines, v["entry"][0], v["entry"][1] + 1)[0][2])
            written = (m.group(1) + (m.group(2) or "")).strip()
            spec.append(f"{ind}{written}, {', '.join(v['attrs'])} :: {v['name']}{v['dims']}")
    if const_entries:
        spec += ["", ind + "! Named constants"]
        for e in const_entries:
            spec += rewrite_entry(e)
    for key, title in titles:
        if not groups[key]:
            continue
        spec += ["", ind + "! " + title]
        bytype = {}
        for nme in groups[key]:
            bytype.setdefault(var[nme]["typ"], []).append(nme)
        for typ in sorted(bytype, key=lambda t: (next((q for q, x in enumerate(TYPE_ORDER) if t.lower().startswith(x)), 9), t.lower())):
            names = sorted(bytype[typ], key=str.lower)
            orig = next(var[q] for q in names)
            # the type as it was written in the source
            m = decl_match(ob.statements(lines, orig["entry"][0], orig["entry"][1] + 1)[0][2])
            written = (m.group(1) + (m.group(2) or "")).strip()
            spec += emit(f"{ind}{written} :: ", [var[q]["name"] + var[q]["dims"] for q in names])
    if verbatim_entries:
        if not a.no_reorder:
            spec += ["", ind + "! Declarations with attributes, initial values or their own comments"]
        for e in verbatim_entries:
            spec += rewrite_entry(e)
    others = [out[i:j + 1] for i, j, kind, payload in entries if kind == "other"]
    for o in others:
        spec += o
    # preprocessor lines of the specification part are kept only around verbatim declarations: refuse otherwise
    if any(kind == "cpp" for i, j, kind, payload in entries) and not a.no_reorder:
        sys.exit("the specification part contains preprocessor conditionals: use --no-reorder for this routine")

    if a.no_reorder:
        # only delete: rebuild the specification part entry by entry, in place
        spec = []
        for i, j, kind, payload in entries:
            if kind == "decl":
                first_name = next(q for q in order if var[q]["entry"] == (i, j))
                if var[first_name]["entry"] == (i, j):
                    spec += rewrite_entry((i, j))
            else:
                spec += out[i:j + 1]

    result = out[:spec_first] + spec + [""] + new_body + out[contains:] if not a.no_reorder else \
        out[:spec_first] + spec + new_body + out[contains:]
    # avoid a growing run of blank lines after the declarations
    print("   unused: " + ", ".join(var[q]["name"] for q in unused[:40]) + (" ..." if len(unused) > 40 else ""))
    if not a.no_reorder:
        print("   groups: " + ", ".join(f"{t.split()[0].lower()} {len(groups[k])}" for k, t in titles)
              + f", kept as declared {sum(1 for q in keep if var[q]['verbatim'])}")
    if a.apply:
        open(a.file, "w").write("\n".join(result))
    else:
        print("   (report only; --apply edits)")


if __name__ == "__main__":
    main()
