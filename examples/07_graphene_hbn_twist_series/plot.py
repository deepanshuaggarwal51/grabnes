#!/usr/bin/env python3
"""Density of states of graphene on hBN for a series of twist angles.

    python3 plot.py                 figure twist_series_dos.png (from the runs, or from reference/dos.dat)
    python3 plot.py --table FILE    only collect the runs into one table (what reference/dos.dat holds)

Table: first line '# twist[deg] ...'; then one row per energy: E [eV], DOS at every angle
[states per eV per atom].
"""

import cmath
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference" / "dos.dat"
OUTPUT = HERE / "twist_series_dos.png"
ENERGY_MIN = -0.6          # Kubo.EnergyMin of the inputs [eV]: fixes the energy unit of the DOS files
A_G = 2.46
HBAR_V = 1.5 * (A_G / math.sqrt(3.0)) * 2.9887       # hbar v_F of the model [eV A]
MISMATCH = -1.0 / 55.0


def collect():
    runs = sorted(HERE.glob("angle_*"), key=lambda p: float(p.name.split("_")[1]))
    runs = [r for r in runs if (r / "generate.DOS").exists()]
    if not runs:
        return None
    angles, columns, energies = [], [], None
    for run in runs:
        rows = [line.split() for line in (run / "generate.DOS").read_text().splitlines() if line.strip()]
        g0 = ENERGY_MIN / float(rows[0][0])       # the DOS file is in units of the hopping g0
        energies = [float(r[0]) * g0 for r in rows]
        columns.append([float(r[1]) / g0 for r in rows])
        angles.append(float(run.name.split("_")[1]))
    return angles, energies, [[c[i] for c in columns] for i in range(len(energies))]


def read_table(path):
    lines = Path(path).read_text().splitlines()
    angles = [float(x) for x in lines[0].split()[2:]]
    rows = [[float(x) for x in line.split()] for line in lines[1:] if line.strip()]
    return angles, [r[0] for r in rows], [r[1:] for r in rows]


def moire_length(angle):
    """Moire period [A] for the lattice mismatch of the series and a twist angle in degrees."""
    return A_G / abs((1.0 + MISMATCH) * cmath.exp(1j * math.radians(angle)) - 1.0)


def main():
    data = collect()
    if len(sys.argv) > 2 and sys.argv[1] == "--table":
        if data is None:
            raise SystemExit("no finished run in this directory")
        angles, energies, dos = data
        with open(sys.argv[2], "w") as handle:
            handle.write("# twist[deg] " + " ".join(f"{a:.2f}" for a in angles) + "\n")
            for energy, row in zip(energies, dos):
                handle.write(f"{energy:.6f} " + " ".join(f"{value:.6e}" for value in row) + "\n")
        print(f"Wrote {sys.argv[2]}")
        return
    source = "the runs" if data is not None else str(REFERENCE)
    angles, energies, dos = data if data is not None else read_table(REFERENCE)

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("plot.py requires matplotlib. Install it with: python3 -m pip install matplotlib") from exc

    figure, axis = plt.subplots(figsize=(7.0, 5.0))
    colors = ["#2457a6", "#2a9d8f", "#e9a03b", "#c0392b", "#6c3483", "#555555"]
    step = 0.004
    for k, angle in enumerate(angles):
        color = colors[k % len(colors)]
        curve = [row[k] + k * step for row in dos]
        axis.plot(energies, curve, color=color, linewidth=1.3,
                  label=f"{angle:.2f} deg, moire {moire_length(angle) / 10:.1f} nm")
        # energy at which the Dirac cone reaches the edge of the moire Brillouin zone (its M point),
        # hbar v_F |G| / 2 with |G| = 4 pi / (sqrt(3) L): where the moire potential opens the secondary gaps
        edge = HBAR_V * 2.0 * math.pi / (math.sqrt(3.0) * moire_length(angle))
        for sign in (-1.0, 1.0):
            if edge < energies[-1]:
                axis.plot([sign * edge], [k * step], marker="^", color=color, markersize=6)
    axis.set_xlabel("Energy (eV)")
    axis.set_ylabel("Density of states (states / eV / atom), curves shifted")
    axis.set_xlim(energies[0], energies[-1])
    axis.set_ylim(bottom=-0.0005)
    axis.legend(frameon=False, fontsize=8, loc="upper center")
    axis.set_title("Graphene on hBN: twist-angle series (triangles: +- hbar v_F |G| / 2)")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source}")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
