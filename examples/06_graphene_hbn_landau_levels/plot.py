#!/usr/bin/env python3
"""Density of states of aligned graphene on hBN against energy and magnetic field.

    python3 plot.py                 figure landau_levels.png (from the run, or from reference/dos_map.dat)
    python3 plot.py --table FILE    only collect the run into one table (what reference/dos_map.dat holds)

Table: first line '# B[T] ...' with the fields; then one row per energy: E [eV], DOS at every field
[states per eV per atom].
"""

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference" / "dos_map.dat"
OUTPUT = HERE / "landau_levels.png"
ENERGY_MIN = -0.3          # Kubo.EnergyMin of Gendata.in [eV]: fixes the energy unit of the DOS files
HBAR_V = 1.5 * (2.46 / math.sqrt(3.0)) * 2.9887      # hbar v_F of the model [eV A]


def collect():
    """(fields, energies, dos[energy][field]) from the directories m<m>_B<B>/ of a finished run."""
    runs = sorted(HERE.glob("m[0-9]*_B*"), key=lambda p: int(p.name[1:].split("_")[0]))
    runs = [r for r in runs if (r / "generate.DOS").exists()]
    if not runs:
        return None
    fields, columns, energies = [], [], None
    for run in runs:
        rows = [line.split() for line in (run / "generate.DOS").read_text().splitlines() if line.strip()]
        # the DOS file is in units of the hopping g0: E/g0 and states per g0
        g0 = ENERGY_MIN / float(rows[0][0])
        energies = [float(r[0]) * g0 for r in rows]
        columns.append([float(r[1]) / g0 for r in rows])
        fields.append(float(run.name.split("_B")[1]))
    return fields, energies, [[c[i] for c in columns] for i in range(len(energies))]


def write_table(path, fields, energies, dos):
    with open(path, "w") as handle:
        handle.write("# B[T] " + " ".join(f"{b:.4f}" for b in fields) + "\n")
        for energy, row in zip(energies, dos):
            handle.write(f"{energy:.6f} " + " ".join(f"{value:.6e}" for value in row) + "\n")


def read_table(path):
    lines = Path(path).read_text().splitlines()
    fields = [float(x) for x in lines[0].split()[2:]]
    rows = [[float(x) for x in line.split()] for line in lines[1:] if line.strip()]
    return fields, [r[0] for r in rows], [r[1:] for r in rows]


def main():
    data = collect()
    if len(sys.argv) > 2 and sys.argv[1] == "--table":
        if data is None:
            raise SystemExit("no finished run in this directory")
        write_table(sys.argv[2], *data)
        print(f"Wrote {sys.argv[2]}")
        return
    source = "this run" if data is not None else str(REFERENCE)
    fields, energies, dos = data if data is not None else read_table(REFERENCE)

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("plot.py requires matplotlib. Install it with: python3 -m pip install matplotlib") from exc

    figure, axis = plt.subplots(figsize=(7.0, 5.0))
    top = sorted(v for row in dos for v in row)[int(0.97 * len(dos) * len(fields))]
    mesh = axis.pcolormesh(fields, energies, dos, shading="nearest", cmap="magma", vmin=0.0, vmax=top)
    figure.colorbar(mesh, ax=axis, label="Density of states (states / eV / atom)")
    # Landau levels of massless Dirac electrons, E_N = sgn(N) v_F sqrt(2 e hbar B |N|)
    for level in range(-4, 5):
        curve = [math.copysign(HBAR_V * math.sqrt(2.0 * abs(level) * b / 6.582119569e4), level) for b in fields]
        axis.plot(fields, curve, color="white", linewidth=0.6, linestyle=":")
    axis.set_xlabel("Magnetic field (T)")
    axis.set_ylabel("Energy (eV)")
    axis.set_ylim(energies[0], energies[-1])
    axis.set_title("Graphene on hBN, both valleys: Landau levels\n(dotted: pristine graphene)")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source}")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
