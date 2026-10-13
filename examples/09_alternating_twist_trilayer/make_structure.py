#!/usr/bin/env python3
"""Rigid alternating-twist trilayer graphene in a commensurate cell, for Run.TypeOfSystem ReadXYZ.

    python3 make_structure.py [m n]          (default 22 21: twist 1.538 degrees, 8322 atoms)

The cell is spanned by L1 = m a1 + n a2 of the outer layers and L2 = L1 rotated by 60 degrees. The middle layer
is the other commensurate orientation of the same cell (L1 = n a1 + m a2 of its own lattice), which is the outer
lattice rotated by the twist angle theta, with cos(theta) = (m^2 + 4 m n + n^2) / (2 (m^2 + m n + n^2)).
All three layers have an atom at the origin (AA stacking there), interlayer distance 3.35 A.

Files written: generate.xyz (cell vectors, number of atoms, element x y z), layerIndex.dat (layer of every
atom, 1 = bottom), sublatticesSorted.dat (atom number and sublattice 1 or 2).
"""
import math
import sys

A_G = 2.46          # lattice constant of graphene [A]
DZ = 3.35           # interlayer distance [A]
Z0 = 10.0           # height of the bottom layer in the cell [A]


def layer(mm, nn, L1, L2, rot):
    """Atoms of a honeycomb layer, rotated by rot, inside the cell (L1, L2)."""
    c, s = math.cos(rot), math.sin(rot)
    turn = lambda v: (c * v[0] - s * v[1], s * v[0] + c * v[1])
    a1, a2 = turn((A_G, 0.0)), turn((A_G / 2, A_G * math.sqrt(3) / 2))
    basis = [(0.0, 0.0), turn((0.0, A_G / math.sqrt(3)))]
    det = L1[0] * L2[1] - L1[1] * L2[0]
    out, seen, rng = [], set(), 3 * (mm + nn)
    for i in range(-rng, rng + 1):
        for j in range(-rng, rng + 1):
            for ib, b in enumerate(basis):
                x, y = i * a1[0] + j * a2[0] + b[0], i * a1[1] + j * a2[1] + b[1]
                f1, f2 = (x * L2[1] - y * L2[0]) / det, (-x * L1[1] + y * L1[0]) / det
                key = (round(f1 % 1.0, 6) % 1.0, round(f2 % 1.0, 6) % 1.0)
                if key in seen or not (-1e-7 <= f1 < 1 - 1e-7 and -1e-7 <= f2 < 1 - 1e-7):
                    continue
                seen.add(key)
                out.append((x, y, ib + 1))
    assert len(out) == 2 * (mm * mm + mm * nn + nn * nn), len(out)
    return out


def main():
    m, n = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (22, 21)
    a1, a2 = (A_G, 0.0), (A_G / 2, A_G * math.sqrt(3) / 2)
    L1 = (m * a1[0] + n * a2[0], m * a1[1] + n * a2[1])
    L2 = (L1[0] / 2 - L1[1] * math.sqrt(3) / 2, L1[0] * math.sqrt(3) / 2 + L1[1] / 2)
    theta = math.acos((m * m + 4 * m * n + n * n) / (2.0 * (m * m + m * n + n * n)))
    # orientation of L1 in a lattice whose a1 is along x: atan2 of (p a1 + q a2)
    ang = lambda p, q: math.atan2(q * math.sqrt(3) / 2, p + q / 2)
    rot_mid = ang(m, n) - ang(n, m)          # the middle layer sees the cell as (n, m)
    atoms = []
    for il, (mm, nn, rot) in enumerate([(m, n, 0.0), (n, m, rot_mid), (m, n, 0.0)]):
        atoms += [(x, y, Z0 + il * DZ, il + 1, sub) for x, y, sub in layer(mm, nn, L1, L2, rot)]
    with open("generate.xyz", "w") as f:
        f.write(f"{L1[0]:.10f} {L1[1]:.10f} 0.0\n{L2[0]:.10f} {L2[1]:.10f} 0.0\n0.0 0.0 {2 * Z0 + 2 * DZ:.4f}\n{len(atoms)}\n")
        f.write("".join(f"C {a[0]:.10f} {a[1]:.10f} {a[2]:.10f}\n" for a in atoms))
    with open("layerIndex.dat", "w") as f:
        f.write("".join(f"{a[3]}\n" for a in atoms))
    with open("sublatticesSorted.dat", "w") as f:
        f.write("".join(f"{i + 1} {a[4]}\n" for i, a in enumerate(atoms)))
    print(f"(m, n) = ({m}, {n}): twist {math.degrees(theta):.4f} degrees (rotation of the middle layer "
          f"{math.degrees(rot_mid):.4f}), moire length {math.hypot(*L1):.2f} A, {len(atoms)} atoms")


if __name__ == "__main__":
    main()
