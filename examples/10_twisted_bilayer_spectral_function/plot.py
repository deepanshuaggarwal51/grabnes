#!/usr/bin/env python3
"""Plot the unfolded spectral function of this example, one panel per layer.

    python3 plot.py                 draw spectral_function.png
    python3 plot.py --table FILE    only write the per-layer summary (what reference/summary.dat holds)
"""

from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "reference" / "spectral.dat"
OUTPUT = HERE / "spectral_function.png"

HOPPING = 3.5           # TB.Hopping [eV]: the energy column of generate.spectral* is in this unit
TITLE = "Twisted bilayer graphene at 6.01 degrees"
LAYERS = [(1, 'layer 1 (rotated)'), (2, 'layer 2 (aligned)')]
KPOINTS = {20: r"$K_1$", 40: r"$K_2$"}     # rows of the path vertices K1 and K2
SMALLEST = 1e-2         # weights below this are not kept in reference/spectral.dat


def read_run() -> dict | None:
    """Weight on the (k, energy) grid per layer, A + B sublattices, from the files of a run."""
    if not (HERE / "generate.spectral1A").exists():
        return None
    run = {}
    for layer, _ in LAYERS:
        a_file = np.loadtxt(HERE / f"generate.spectral{layer}A")
        b_file = np.loadtxt(HERE / f"generate.spectral{layer}B")
        number_of_energies = len(np.unique(a_file[:, 1]))
        shape = (len(a_file) // number_of_energies, number_of_energies)
        run["distance"] = a_file[:, 0].reshape(shape)[:, 0]
        run["energy"] = a_file[:, 1].reshape(shape)[0] * HOPPING
        run[layer] = (a_file[:, 2] + b_file[:, 2]).reshape(shape)
    return run


def read_reference() -> dict:
    """The same grid from reference/spectral.dat: its rows hold the non-negligible weights."""
    header = REFERENCE.read_text().splitlines()[:2]
    distance = np.array([float(value) for value in header[0].split()[2:]])
    energy = np.array([float(value) for value in header[1].split()[2:]])
    run = {"distance": distance, "energy": energy}
    for layer, _ in LAYERS:
        run[layer] = np.zeros((len(distance), len(energy)))
    for layer, k_row, energy_row, weight in np.loadtxt(REFERENCE):
        run[int(layer)][int(k_row), int(energy_row)] = weight
    return run


def write_reference(run: dict, path: Path) -> None:
    """Keep the non-negligible weights of a run: the figure is drawn from these."""
    with open(path, "w") as out:
        out.write("# distance " + " ".join(f"{value:.6f}" for value in run["distance"]) + "\n")
        out.write("# energy " + " ".join(f"{value:.6f}" for value in run["energy"]) + "\n")
        for layer, _ in LAYERS:
            for k_row, energy_row in zip(*np.nonzero(run[layer] > SMALLEST)):
                out.write(f"{layer} {k_row} {energy_row} {run[layer][k_row, energy_row]:.4f}\n")


def write_summary(run: dict, path: Path) -> None:
    """Per k-point and layer: the centre [eV] and the total of the weight on the energy points.
    One row per k and layer, so that it can be compared between runs whatever the number of peaks."""
    with open(path, "w") as out:
        for layer, _ in LAYERS:
            for k_row, weights in enumerate(run[layer]):
                total = weights.sum()
                centre = (weights * run["energy"]).sum() / total if total > 0.0 else 0.0
                out.write(f"{layer} {k_row} {centre:.5f} {total:.3f}\n")


def main() -> None:
    run = read_run()
    if len(sys.argv) > 2 and sys.argv[1] == "--table":
        write_summary(run, Path(sys.argv[2]))
        return
    if len(sys.argv) > 2 and sys.argv[1] == "--reference":
        write_reference(run, Path(sys.argv[2]))
        return
    source = "this run" if run is not None else str(REFERENCE)
    if run is None:
        run = read_reference()

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "plot.py requires matplotlib. Install it with: python3 -m pip install matplotlib"
        ) from exc

    figure, axes = plt.subplots(
        1, len(LAYERS), figsize=(3.6 * len(LAYERS) + 1.0, 4.4), sharey=True, constrained_layout=True
    )
    for axis, (layer, label) in zip(axes, LAYERS):
        image = axis.pcolormesh(
            run["distance"], run["energy"], np.log10(1.0 + run[layer].T),
            cmap="Blues", vmin=0.0, vmax=2.3, shading="nearest", rasterized=True,
        )
        positions = [run["distance"][row] for row in KPOINTS]
        for position in positions:
            axis.axvline(position, color="0.45", linewidth=0.7, linestyle=":")
        axis.set_xticks(positions, list(KPOINTS.values()))
        axis.set_xlabel("Cut through the two Dirac points")
        axis.set_title(label, fontsize=11)
    axes[0].set_ylabel("Energy (eV)")
    figure.colorbar(image, ax=list(axes), shrink=0.85, label=r"$\log_{10}(1 + \mathrm{unfolded\ weight})$")
    figure.suptitle(TITLE)
    figure.savefig(OUTPUT, dpi=180)
    print(f"Read {source}")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
