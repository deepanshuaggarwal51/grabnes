#!/usr/bin/env python3
"""Plot the moire bands of aligned graphene on hBN (one valley) from this example."""

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
OUTPUT = HERE / "graphene_hbn_bands.png"

# Rows of the path vertices Gamma, K, M, Gamma: the legs get points in proportion to their length
# (Bands.NumPoints 40 on the first leg; K-M is half and M-Gamma sqrt(3)/2 of Gamma-K).
VERTICES = [0, 40, 60, 95]
LABELS = [r"$\Gamma$", r"$K$", r"$M$", r"$\Gamma$"]
WINDOW = (-0.35, 0.35)


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

    figure, axis = plt.subplots(figsize=(6.4, 4.6))
    for band in energies:
        if max(band) < WINDOW[0] or min(band) > WINDOW[1]:
            continue
        axis.plot(x_values, band, color="#2457a6", linewidth=1.5)

    positions = [x_values[i] for i in VERTICES]
    for position in positions:
        axis.axvline(position, color="0.82", linewidth=0.8, zorder=0)
    axis.set_xticks(positions, LABELS)
    axis.set_xlim(positions[0], positions[-1])
    axis.set_ylim(*WINDOW)
    axis.set_xlabel("Wave-vector path in the moire Brillouin zone")
    axis.set_ylabel("Energy (eV)")
    axis.set_title("Aligned graphene on hBN, effective moire model, one valley")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source}")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
