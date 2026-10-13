#!/usr/bin/env python3
"""Plot the bands of one valley around the flat bands from this example."""

from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "plot.py requires matplotlib. Install it with: python3 -m pip install matplotlib"
    ) from exc


HERE = Path(__file__).resolve().parent
CALCULATED = HERE / "generate.bands"
REFERENCE = HERE / "reference" / "bands.dat"
OUTPUT = HERE / "trilayer_bands.png"

# Rows of the path vertices K, Gamma, M, K': the legs get points in proportion to their length
# (Bands.NumPoints 20 on the first leg; Gamma-M is sqrt(3)/2 and M-K' half of K-Gamma).
VERTICES = [0, 20, 37, 47]
LABELS = [r"$K$", r"$\Gamma$", r"$M$", r"$K'$"]
WINDOW = 0.12          # half-width of the energy window around neutrality [eV]


def read_bands(path: Path) -> tuple[list[float], list[list[float]]]:
    """Read a GRABNES .bands file: three header lines, 'bands spins k-points', then for every
    k-point its position along the path followed by the levels (possibly wrapped over lines)."""
    lines = path.read_text().splitlines()
    number_of_bands, _, number_of_k = (int(value) for value in lines[3].split()[:3])
    values = [float(value) for line in lines[4:] for value in line.split()]
    width = number_of_bands + 1
    if len(values) != number_of_k * width:
        raise ValueError(f"{path} does not hold {number_of_k} k-points of {number_of_bands} levels")
    rows = [values[i * width:(i + 1) * width] for i in range(number_of_k)]
    return [row[0] for row in rows], [[row[1 + b] for row in rows] for b in range(number_of_bands)]


def main() -> None:
    source = CALCULATED if CALCULATED.exists() else REFERENCE
    x_values, energies = read_bands(source)

    # Reference energy: the flat bands at the first k-point (K), where they are the two middle levels
    # of the valley. Along the path the Dirac cone of the trilayer crosses them, so the bands are not
    # coloured by their index.
    half = len(energies) // 2
    flat_energy = 0.5 * (energies[half - 1][0] + energies[half][0])

    figure, axis = plt.subplots(figsize=(6.4, 4.6))
    for band in energies:
        shifted = [1000.0 * (energy - flat_energy) for energy in band]
        if max(shifted) < -1000.0 * WINDOW or min(shifted) > 1000.0 * WINDOW:
            continue
        axis.plot(x_values, shifted, color="#2457a6", linewidth=1.3)

    positions = [x_values[i] for i in VERTICES]
    for position in positions:
        axis.axvline(position, color="0.82", linewidth=0.8, zorder=0)
    axis.set_xticks(positions, LABELS)
    axis.set_xlim(positions[0], positions[-1])
    axis.set_ylim(-1000.0 * WINDOW, 1000.0 * WINDOW)
    axis.set_xlabel("Wave-vector path in the moire Brillouin zone")
    axis.set_ylabel("Energy - flat bands at K (meV)")
    axis.set_title("Alternating-twist trilayer graphene at 1.54 degrees, one valley")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source}")
    print(f"Flat bands at K: {flat_energy:.4f} eV")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
