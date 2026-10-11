#!/usr/bin/env python3
"""Move one branch of a parallel loop into an internal procedure of the routine.

    python3 tools/maintenance/outline_branch.py FILE ROUTINE LINE NAME "what the branch does" [--apply]

LINE is the header of a branch (`if (...) then`, `else if (...) then` or `else`) inside an OpenMP parallel
loop of ROUTINE. The statements of that branch, up to the next `else` / `end if` of the same block, are
moved unchanged into an internal procedure NAME, and the branch keeps one line that calls it.

Inside a parallel region an internal procedure sees the variables of its routine, but for a variable that
is PRIVATE to the loop it would see the shared original and not the copy of the thread. Therefore every
variable of the branch that has a copy per thread (the PRIVATE, FIRSTPRIVATE and REDUCTION clauses, and the
loop indices) becomes an argument of the procedure, declared as in the routine with intent(inout) (the loop indices with intent(in)); the
procedure then works on the copies of the calling thread. Everything else (the model parameters, the
arrays of the Hamiltonian) is used from the routine as before.

A `cycle` of the surrounding loop cannot stay in a procedure. If the branch contains one, NAME becomes a
logical function that returns .true. where the branch said `cycle`, and the call reads
`if (NAME(...)) cycle`.

The tool refuses a branch that contains `return`, `exit` of an outer loop, `goto`, statement labels or
unbalanced preprocessor conditionals, and a branch that assigns a local variable of the routine that is
neither private nor an argument (such an assignment would be a data race).
"""
import argparse
import re
import sys

import outline_block as ob

DECL = re.compile(r"^\s*(integer|real|complex|logical|character|double\s+precision)\s*(\([^)]*\)|\*\s*\d+)?\s*((?:,[^:]*)?)::\s*(.*)$", re.I)


OLD_DECL = re.compile(r"^\s*(integer|real|complex|logical|character|double\s+precision)\s*(\([^)]*\)|\*\s*\d+)?\s+([A-Za-z_]\w*.*)$", re.I)


def decl_match(text):
    """Match object of a type declaration; the old form without `::` (integer u) is accepted too."""
    m = DECL.match(text)
    if m:
        return m
    m = OLD_DECL.match(text)
    if m and "::" not in text and not re.match(r"^\s*\w+(\s*\([^)]*\))?\s+function\b", text, re.I):
        return DECL.match(f"{m.group(1)}{m.group(2) or ''} :: {m.group(3)}")
    return None


def split_top(text):
    depth, item, items = 0, "", []
    for ch in text:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            items.append(item)
            item = ""
        else:
            item += ch
    items.append(item)
    return [x.strip() for x in items if x.strip()]


def declarations(lines, a, b):
    """{name: (type, attributes, dimensions)} of the variables declared in lines a..b-1."""
    out = {}
    for i, j, text, is_omp in ob.statements(lines, a, b):
        if is_omp or not text or text.startswith("#"):
            continue
        m = decl_match(text)
        if not m:
            continue
        typ = m.group(1) + (m.group(2) or "")
        attrs = [x for x in split_top(m.group(3).lstrip(",")) if x.lower() not in ("save",)]
        for item in split_top(m.group(4)):
            mm = re.match(r"^(\w+)\s*(\(.*\))?\s*(=.*)?$", item)
            if mm:
                out[mm.group(1).lower()] = (typ, attrs, mm.group(2) or "", mm.group(1))
    return out


def clause_names(directive):
    """Variables that have a copy per thread: PRIVATE, FIRSTPRIVATE and REDUCTION clauses of a directive."""
    names = set()
    for clause in re.findall(r"(?:first)?private\s*\(([^)]*)\)", directive):
        names |= {x.strip() for x in clause.split(",") if x.strip()}
    for clause in re.findall(r"reduction\s*\(\s*[^:]*:([^)]*)\)", directive):
        names |= {x.strip() for x in clause.split(",") if x.strip()}
    return names


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("file")
    p.add_argument("routine")
    p.add_argument("line", type=int)
    p.add_argument("name")
    p.add_argument("doc")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--locals", action="store_true",
                   help="real and complex variables that the routine uses in this branch only become local "
                        "variables of the procedure, marked as not set on entry")
    p.add_argument("--unset", default="hamUnset", help="variable of the routine that holds the 'not set' value")
    p.add_argument("--shared-ok", nargs="*", default=[],
                   help="shared variables the branch may assign (written by one thread only, e.g. a print-once flag)")
    p.add_argument("--host", default="",
                   help="ROUTINE is itself an internal procedure of HOST: the declarations are those of HOST and the "
                        "new procedure becomes another internal procedure of HOST, placed after ROUTINE")
    a = p.parse_args()

    lines = open(a.file).read().split("\n")
    r0 = next(i for i, l in enumerate(lines) if re.match(r"^\s*subroutine\s+" + a.routine + r"\b", l, re.I))
    r1 = next(i for i, l in enumerate(lines) if i > r0 and re.match(r"^\s*end\s+subroutine\s+" + a.routine + r"\b", l, re.I))
    if a.host:
        x0, x1 = r0, r1                  # the internal procedure that holds the loop
        r0 = next(i for i, l in enumerate(lines) if re.match(r"^\s*subroutine\s+" + a.host + r"\b", l, re.I))
        r1 = next(i for i, l in enumerate(lines) if i > r0 and re.match(r"^\s*end\s+subroutine\s+" + a.host + r"\b", l, re.I))
        contains = next((i for i in range(r0, r1) if re.match(r"^\s*contains\s*$", ob.code(lines[i]), re.I)), None)
        if contains is None or not contains < x0 < x1 < r1:
            sys.exit(f"{a.routine} is not an internal procedure of {a.host}")
        if any(decl_match(ob.code(l)) for l in lines[x0 + 1:x1]):
            sys.exit(f"{a.routine} has declarations of its own: not supported with --host")
        body_lo, body_hi, insert_at, else_hi = x0 + 1, x1, x1 + 1, r1
    else:
        contains = next((i for i in range(r0, r1) if re.match(r"^\s*contains\s*$", ob.code(lines[i]), re.I)), None)
        if contains is None:
            sys.exit("the routine has no internal procedures yet (use outline_block.py first)")
        body_lo, body_hi, insert_at, else_hi = r0 + 1, contains, r1, contains
    if any(re.match(r"^\s*(subroutine|function)\s+" + a.name + r"\b", l, re.I) for l in lines):
        sys.exit(f"a procedure named {a.name} exists already")
    decl = declarations(lines, r0 + 1, contains)

    st = ob.statements(lines, body_lo, body_hi)
    stack = []                       # open blocks: [kind, index of the header statement]
    private, pending_private = set(), None
    target = None
    for k, (i, j, text, is_omp) in enumerate(st):
        if is_omp:
            t = text.lower()
            if re.match(r"!\$omp\s+parallel\s+do\b", t):
                pending_private = clause_names(t)
            elif re.match(r"!\$omp&", t) and pending_private is not None:
                pending_private |= clause_names(t)
            continue
        if not text or text.startswith("#"):
            continue
        text = re.sub(r"^\d+\s+", "", text)
        if i == a.line - 1:
            target = k
            break
        if ob.OPEN.match(text):
            kind = "do" if re.match(r"^(\w+\s*:\s*)?do\b", text, re.I) else "if"
            entry = [kind, k, None, None]
            if kind == "do":
                m = re.match(r"^(?:\w+\s*:\s*)?do\s+(\w+)\s*=", text, re.I)
                entry[2] = m.group(1).lower() if m else None
                entry[3] = pending_private
                pending_private = None
            stack.append(entry)
        elif ob.CLOSE.match(text):
            stack.pop()
    if target is None:
        sys.exit(f"line {a.line} is not the start of a statement of {a.routine}")
    head = re.sub(r"^\d+\s+", "", st[target][2])
    is_if = bool(re.match(r"^if\s*\(.*\)\s*then$", head, re.I))
    is_else = bool(re.match(r"^else\b", head, re.I))
    if not (is_if or is_else):
        sys.exit(f"line {a.line} is not the header of a branch: {head[:60]}")
    loops = [e for e in stack if e[0] == "do"]
    par = [e for e in loops if e[3] is not None]
    if not par:
        sys.exit("the branch is not inside an OpenMP parallel loop (use outline_block.py)")
    indices = {e[2] for e in loops if e[2]}
    priv = set(par[-1][3]) | indices

    # extent of the branch body
    depth, end_k = 0, None
    for k in range(target + (1 if True else 0), len(st)):
        i, j, text, is_omp = st[k]
        if is_omp or not text or text.startswith("#"):
            continue
        text = re.sub(r"^\d+\s+", "", text)
        if k == target:
            continue
        if ob.OPEN.match(text):
            depth += 1
        elif ob.CLOSE.match(text):
            if depth == 0:
                end_k = k
                break
            depth -= 1
        elif depth == 0 and re.match(r"^else\b", text, re.I):
            end_k = k
            break
    if end_k is None:
        sys.exit("end of the branch not found")
    first, last = st[target][1] + 1, st[end_k][0] - 1
    while last >= first and not lines[last].strip():
        last -= 1
    if last < first:
        sys.exit("the branch is empty")
    body = lines[first:last + 1]

    # checks, and the places where the branch leaves the surrounding loop
    cpp, loop, cycles = 0, 0, []
    for i, j, text, is_omp in ob.statements(lines, first, last + 1):
        if text.startswith("#"):
            cpp += 1 if re.match(r"#\s*if", text) else (-1 if re.match(r"#\s*endif", text) else 0)
            continue
        if is_omp:
            sys.exit(f"line {i + 1}: OpenMP directive inside the branch")
        if not text:
            continue
        t = text.lower()
        if re.match(r"^(\w+\s*:\s*)?do\b", t):
            loop += 1
        elif re.match(r"^end\s*do\b", t):
            loop -= 1
        if re.search(r"(^|[\s)])return\s*$", t) or re.search(r"\bgo\s*to\b", t) or re.match(r"^\d+\s", t):
            sys.exit(f"line {i + 1}: return, goto or statement label in the branch")
        if re.search(r"(^|[\s)])exit\b", t) and loop <= 0:
            sys.exit(f"line {i + 1}: exit of a loop outside the branch")
        if re.search(r"(^|[\s)])cycle\b", t) and loop <= 0:
            if i != j or not (re.match(r"^cycle$", t) or re.match(r"^if\s*\(.*\)\s*cycle$", t)):
                sys.exit(f"line {i + 1}: cycle in a form the tool does not handle: {text[:60]}")
            cycles.append(i)
    if cpp != 0:
        sys.exit("unbalanced preprocessor conditionals in the branch")

    text_body = "\n".join(ob.code(l) for l in body).lower()
    used = []
    for m in re.finditer(r"(?<![\w%])([a-z_]\w*)", text_body):
        n = m.group(1)
        if n in priv and n not in used:
            used.append(n)
    missing = [n for n in used if n not in decl]
    if missing:
        sys.exit("private variables without a declaration found in the routine: " + ", ".join(missing))
    # assignments to local variables that are not private: a data race, and not something to move
    for i, j, text, is_omp in ob.statements(lines, first, last + 1):
        if is_omp or not text or text.startswith("#"):
            continue
        m = re.match(r"^(?:if\s*\(.*?\)\s+)?(\w+)\s*(\(.*?\))?\s*=(?!=)", text)
        if m and m.group(1).lower() in decl and m.group(1).lower() not in priv and not m.group(2):
            if m.group(1).lower() not in ("if",) and m.group(1).lower() not in {x.lower() for x in a.shared_ok}:
                sys.exit(f"line {i + 1}: the branch assigns '{m.group(1)}', which is not private to the loop")

    # variables that the routine uses nowhere else: candidates for local variables of the procedure
    local = []
    if a.locals:
        marker = False
        elsewhere = []
        for n in range(r0 + 1, else_hi):
            l = lines[n]
            if "! >>> unset markers" in l:
                marker = True
            if first <= n <= last or marker or l.lstrip().lower().startswith("!$omp") or decl_match(ob.code(l)):
                if "! <<< unset markers" in l:
                    marker = False
                continue
            elsewhere.append(ob.code(l).lower())
        elsewhere = "\n".join(elsewhere)
        for n in used:
            typ = decl[n][0].lower()
            if (typ.startswith("real") or typ.startswith("complex")) and not decl[n][1] \
                    and not re.search(r"(?<![\w%])" + re.escape(n) + r"(?!\w)", elsewhere):
                local.append(n)
        used = [n for n in used if n not in local]

    function = bool(cycles)
    new_body = []
    for idx, l in enumerate(body):
        n = first + idx
        if n in cycles:
            ind = re.match(r"\s*", l).group(0)
            c = ob.code(l).strip()
            mm = re.match(r"^if\s*\((.*)\)\s*cycle$", c, re.I)
            if mm:
                new_body += [f"{ind}if ({mm.group(1)}) then", f"{ind}   skipBond = .true.", f"{ind}   return", f"{ind}end if"]
            else:
                new_body += [f"{ind}skipBond = .true.", f"{ind}return"]
        else:
            new_body.append(l)

    args = [decl[n][3] for n in used]
    arglist = ", ".join(args)
    decls = []
    for n in used:
        typ, attrs, dims, name = decl[n]
        extra = ", ".join(attrs)
        # a loop index must not be modified inside its loop: it is passed for reading only
        intent = "intent(in)" if n in indices else "intent(inout)"
        decls.append(f"   {typ}{', ' + extra if extra else ''}, {intent} :: {name}{dims}")

    def wrapped(prefix, text, cont_indent):
        out, cur = [], prefix
        for k, piece in enumerate(text.split(", ")):
            piece = piece + (", " if k < text.count(", ") else "")
            if len(cur) + len(piece) > 118:
                out.append(cur.rstrip() + " &")
                cur = cont_indent
            cur += piece
        return out + [cur]

    ind = re.match(r"\s*", lines[first]).group(0)
    if function:
        call = wrapped(f"{ind}if ({a.name}(", arglist + ")) cycle", ind + "      ")
        header = wrapped(f"function {a.name}(", arglist + ") result(skipBond)", "      ")
    else:
        call = wrapped(f"{ind}call {a.name}(", arglist + ")", ind + "      ")
        header = wrapped(f"subroutine {a.name}(", arglist + ")", "      ")
    kind = "function" if function else "subroutine"
    proc = ["!> " + a.doc,
            "!! Internal procedure of " + (a.host or a.routine) + ", called inside the parallel loop" + (" of " + a.routine if a.host else " for one bond") + ". The variables that are",
            "!! private to the loop are its arguments; the model parameters and the arrays of the Hamiltonian are those",
            "!! of " + (a.host or a.routine) + "." + (" The result is .true. where the bond is to be skipped (cycle)." if function else "")]
    ldecls, linit = [], []
    for n in local:
        typ, attrs, dims, name = decl[n]
        ldecls.append(f"   {typ} :: {name}{dims}")
        linit.append(f"   {name} = {a.unset}")
    if local:
        proc += ["!! The variables declared after the arguments are used by this model only; they are marked as not set on",
                 "!! entry, so that a value left from another bond cannot be used by mistake."]
    proc += header + [""] + decls + (["   logical :: skipBond"] if function else []) + ([""] + ldecls if ldecls else [])
    proc += [""] + (["   skipBond = .false."] if function else []) + linit + [""]
    proc += new_body + ["", f"end {kind} {a.name}", ""]

    new = lines[:first] + call + lines[last + 1:]
    shift = len(call) - (last - first + 1)
    new[insert_at + shift:insert_at + shift] = ([""] if a.host else []) + proc
    print(f"{a.routine}: lines {first + 1}-{last + 1} ({last - first + 1} lines) -> internal {kind} {a.name} with "
          f"{len(args)} arguments" + (f", {len(local)} local variables" if local else "")
          + (f", {len(cycles)} cycle statement(s)" if cycles else "")
          + ("" if a.apply else "   (report only; --apply edits)"))
    print("   arguments: " + arglist)
    if a.apply:
        open(a.file, "w").write("\n".join(new))


if __name__ == "__main__":
    sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
    main()
