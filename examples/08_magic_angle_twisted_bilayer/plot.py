#!/usr/bin/env python3
"""Plot the bands of one valley around charge neutrality from this example."""

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
OUTPUT = HERE / "magic_angle_bilayer_bands.png"

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

    # Charge neutrality of one valley: half of the plane-wave states are filled. The two bands next
    # to it are the flat bands.
    half = len(energies) // 2
    flat = [energies[half - 1], energies[half]]
    neutrality = 0.5 * (max(flat[0]) + min(flat[1]))
    width = 1000.0 * (max(flat[1]) - min(flat[0]))

    figure, axis = plt.subplots(figsize=(6.4, 4.6))
    for index, band in enumerate(energies):
        shifted = [1000.0 * (energy - neutrality) for energy in band]
        if max(shifted) < -1000.0 * WINDOW or min(shifted) > 1000.0 * WINDOW:
            continue
        is_flat = index in (half - 1, half)
        axis.plot(x_values, shifted, color="#c0392b" if is_flat else "#2457a6", linewidth=1.8 if is_flat else 1.2)

    positions = [x_values[i] for i in VERTICES]
    for position in positions:
        axis.axvline(position, color="0.82", linewidth=0.8, zorder=0)
    axis.set_xticks(positions, LABELS)
    axis.set_xlim(positions[0], positions[-1])
    axis.set_ylim(-1000.0 * WINDOW, 1000.0 * WINDOW)
    axis.set_xlabel("Wave-vector path in the moire Brillouin zone")
    axis.set_ylabel("Energy - neutrality (meV)")
    axis.set_title("Twisted bilayer graphene at 1.08 degrees, one valley")
    figure.tight_layout()
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source}")
    print(f"Neutrality at {neutrality:.4f} eV; the two flat bands span {width:.1f} meV along this path")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
