# 07: Graphene on hBN, a series of twist angles

The density of states of graphene on hBN for four twist angles between the
two crystals, 0, 0.61, 0.92 and 1.87 degrees, with the effective moire model
of examples 05 and 06. The moire period shrinks from 13.5 nm to 6.6 nm, and
the dip that the moire potential makes in the density of states on the hole
side moves away from the Dirac point accordingly.

## Expected result

![07: Graphene on hBN, a series of twist angles](reference/twist_series_dos.png)

This is the figure `plot.py` draws from the reference data. After `./run.sh`, `python3 plot.py`
writes `twist_series_dos.png` in this folder from your own run; the two figures should look the same.

## Commensurate cells

The effective model adds a moire potential to a graphene cell, so the moire
pattern must be periodic in that cell. For N x N graphene cells this holds
when the hBN lattice vector satisfies `N b1 = p a1 + q a2` with integers p
and q, which fixes the lattice mismatch (`Moire.LatticePercent`) and the twist
angle (`Moire.TwistAngle`, in degrees, with `Moire.Twisted .true.`):

| Input | N | (p, q) | mismatch | twist | moire period | SuperCell | atoms |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Gendata_0.00.in` | 55 | (54, 0) | -1.818 % | 0 | 13.53 nm | 24 | 3 484 800 |
| `Gendata_0.61.in` | 83 | (81, 1) | -1.802 % | 0.609 deg | 11.79 nm | 16 | 3 527 168 |
| `Gendata_0.92.in` | 110 | (107, 2) | -1.806 % | 0.919 deg | 10.23 nm | 12 | 3 484 800 |
| `Gendata_1.87.in` | 27 | (26, 1) | -1.799 % | 1.872 deg | 6.64 nm | 48 | 3 359 232 |

`make_inputs.py` computes these numbers and writes the four inputs from
`Gendata.template`; edit its table to add angles. The two keys must be given
with all their digits: a rounded value makes the potential discontinuous at
the cell boundary.

## Method

The recursion (Lanczos) method of example 06 at zero field: 3000 steps, a
broadening of 5 meV, one random vector. Without a magnetic field the spectrum
is continuous, and a smooth density of states from one random vector needs a
large system: `Structure.SuperCell` repeats each cell to about 3.5 million
atoms. With a tenth of that, the curve is dominated by the discrete levels of
the supercell.

## What to look for

- A dip on the hole side at -0.17 eV for zero twist, moving to -0.20, -0.23
  and -0.35 eV as the twist grows.
- The triangles in the figure mark the energy at which the Dirac cone of
  pristine graphene reaches the edge of the moire Brillouin zone,
  hbar v_F |G| / 2 with |G| = 4 pi / (sqrt(3) L) for a moire period L. The
  dips follow them.
- On the electron side the feature is much weaker, as in example 05.

## Run

```sh
OMP_NUM_THREADS=4 ./run.sh
python3 plot.py
```

With four threads the four runs take about five minutes and 1.5 GB of memory.
Each run is written to its own directory `angle_<degrees>/`. The plotting
step writes `twist_series_dos.png`; without a run it uses `reference/dos.dat`
(energy in eV, then the density of states at every angle).
`python3 plot.py --table FILE` writes the same table from the runs.
