# 09: Alternating-twist trilayer graphene, structure read from files

Three graphene layers in which the middle one is twisted by 1.54 degrees
against the two outer ones (8322 atoms). The structure is not built by the
solver: it is read from files (`Run.TypeOfSystem ReadXYZ`), which is how a
relaxed structure from a molecular-dynamics code is used. Here a small script
writes a rigid structure, so the example needs no other program.

## Expected result

![09: Alternating-twist trilayer graphene, structure read from files](reference/trilayer_bands.png)

This is the figure `plot.py` draws from the reference data. After `./run.sh`, `python3 plot.py`
writes `trilayer_bands.png` in this folder from your own run; the two figures should look the same.

## The structure files

`make_structure.py` writes, for the commensurate cell `(m, n) = (22, 21)`:

| File | Content |
| --- | --- |
| `generate.xyz` | three lines with the cell vectors [A], the number of atoms, then `element x y z` per atom |
| `layerIndex.dat` | the layer of every atom (1 = bottom), read with `Structure.ReadLayerIndex .true.` |
| `sublatticesSorted.dat` | atom number and sublattice (1 or 2) |

The atoms keep the order of the file. `python3 make_structure.py m n` writes
another commensurate cell; `run.sh` calls the script when the files are
missing.

## The model

The model of example 08 (F2G2 intralayer, two-center interlayer with factor
1.008 on a rigid lattice), with `Stack.ThreeLayers` and `Stack.MiddleTwist`:
single-layer intralayer parameters in every layer. The twist is sqrt(2)
times the angle of example 08, the relation between the flat-band angles of
the trilayer and of the bilayer.

This is the model family of Leconte, Park, An, Samudrala and Jung, 2D Mater.
9, 044002 (2022), run here on a rigid lattice.

## What to look for

- Two flat bands at 0.191 eV, 3 meV wide along the path, as in the bilayer.
- In addition the Dirac cone of a single graphene layer, which the trilayer
  keeps because of its mirror symmetry. For this valley it sits at the corner
  `K'` = (1/3, 2/3) of the moire zone; the path ends there.
- In this rigid model the Dirac point lies 26 meV BELOW the flat bands, at
  0.165 eV, the Dirac energy of the single-layer parameters. The offset is a
  property of this parameter set on a rigid lattice, not a general statement
  about the material.

## The plane-wave reduction

`Diag.N_G 8` gives 330 states (55 moire reciprocal vectors, three layers, two
sublattices). `Diag.MoireAngle 0.0` because the first layer of the file has
its lattice vector along x and the cell is not rotated; a structure whose
cell was rotated needs the angle of that rotation.

Compared with the direct diagonalization of the same Hamiltonian in both
directions (October 2026): at six k-points, 135 of 135 levels within 0.3 eV;
differences up to 0.02 meV within 50 meV and 0.11 meV within 0.3 eV.

## Run

```sh
./run.sh
python3 plot.py
```

The calculation takes about three minutes. It writes `generate.bands` and
`job.out`; the plotting step writes `trilayer_bands.png`. If no calculated
file exists, the plotter uses `reference/bands.dat`.
