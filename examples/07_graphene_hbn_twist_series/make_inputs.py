#!/usr/bin/env python3
"""Write the inputs of the twist-angle series: Gendata_<angle>.in, one per twist angle.

The effective moire model needs a moire pattern that is periodic in the simulation cell. With graphene
lattice vectors a1, a2 (60 degrees apart) and N x N graphene cells, the pattern of an hBN layer with
lattice vectors b1, b2 is periodic if

    N b1 = p a1 + q a2        (p, q integers).

This fixes the lattice mismatch and the twist angle of the hBN layer:

    1 + eps = |p + q exp(i pi/3)| / N ,      theta = arg(p + q exp(i pi/3)).

The series below keeps eps within 1 % of -1/55 (the aligned case of example 05) and takes the twist
angles this allows in cells of at most 110 x 110 graphene cells. `Structure.SuperCell` repeats the cell
so that every run has about 3.5 million atoms.
"""
import cmath
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
A_G = 2.46
#        N    p   q  SuperCell
CASES = [(55, 54, 0, 24),
         (83, 81, 1, 16),
         (110, 107, 2, 12),
         (27, 26, 1, 48)]


def number(value, digits):
    """The solver reads a negative number as '0-x'."""
    return f"0-{abs(value):.{digits}f}" if value < 0 else f"{value:.{digits}f}"


def main():
    template = (HERE / "Gendata.template").read_text()
    for cells, p, q, supercell in CASES:
        w = p + q * cmath.exp(1j * math.pi / 3)
        eps, theta = abs(w) / cells - 1.0, math.degrees(cmath.phase(w))
        moire = A_G / abs((1.0 + eps) * cmath.exp(1j * math.radians(theta)) - 1.0)
        text = (template.replace("@CELLSIZE@", str(cells)).replace("@SUPERCELL@", str(supercell))
                .replace("@EPS@", number(eps, 15)).replace("@THETA@", number(theta, 12))
                .replace("@TWISTED@", ".true." if q else ".false."))
        name = f"Gendata_{theta:.2f}.in"
        (HERE / name).write_text(text)
        print(f"{name}: {cells} x {cells} graphene cells, N b1 = {p} a1 + {q} a2, mismatch {100 * eps:.3f} %, "
              f"twist {theta:.4f} degrees, moire length {moire / 10:.2f} nm, "
              f"{2 * (cells * supercell) ** 2} atoms")


if __name__ == "__main__":
    main()
