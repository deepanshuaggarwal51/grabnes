# 08: Twisted bilayer graphene near the magic angle

Twisted bilayer graphene at 1.08 degrees: the commensurate cell
`(m, n) = (31, 30)` with 11 164 atoms, built by the solver, and the bands of
one valley along `K - Gamma - M - K'` of the moire Brillouin zone. The two
bands at charge neutrality are flat: they span 2.2 meV along the path.

## Expected result

![08: Twisted bilayer graphene near the magic angle](reference/magic_angle_bilayer_bands.png)

This is the figure `plot.py` draws from the reference data. After `./run.sh`, `python3 plot.py`
writes `magic_angle_bilayer_bands.png` in this folder from your own run; the two figures should look the same.

## The model

- Intralayer: the F2G2 model of graphene with five neighbor shells
  (`Intralayer.F2G2Model`, `TB.NeighLevels 5`), with the first-neighbor
  hopping set to 3.5 eV and the site energy to 0.4914 eV (`TB.Hopping`,
  `TB.OnSiteC`). The site energy only sets the zero of energy: neutrality
  lies at 0.178 eV.
- Interlayer: the two-center form with a decay that depends on the distance
  (`Interlayer.TypeOfBL Koshino`, `Interlayer.KoshinoSR`), multiplied by one
  factor (`Interlayer.RenormalizeCoupling`, `Interlayer.CouplingFactor`).
- The structure is rigid: flat layers 3.35 A apart, no lattice relaxation.
  The factor 1.008 puts the flattest bands of the RIGID lattice near 1.08
  degrees. With relaxed structures the same model is used with 0.895.

This is the model family of Leconte, Javvaji, An, Samudrala and Jung, Phys.
Rev. B 106, 115410 (2022), run here on a rigid lattice. Lattice relaxation
changes the result qualitatively in one respect: it opens gaps between the
flat bands and the bands above and below. In this rigid calculation those
bands come within 3 meV of the flat bands at `Gamma`.

## The plane-wave reduction

As in example 05, the Hamiltonian is projected on plane waves around one
graphene valley (`Diag.N_G 8`: 55 moire reciprocal vectors, two layers, two
sublattices, 220 states). `Diag.MoireAngle 0.0` is the orientation of the
first layer in the frame of the cell that the solver builds.

Compared with the direct diagonalization of the same 11 164-atom Hamiltonian
in both directions (October 2026): at six k-points the levels of the two
valleys together are the exact levels within 0.3 eV of neutrality, 171 of
171; differences up to 0.02 meV within 50 meV and 0.18 meV within 0.3 eV.

## Run

```sh
./run.sh
python3 plot.py
```

The calculation takes about two minutes. It writes `generate.bands` and
`job.out`; the plotting step writes `magic_angle_bilayer_bands.png` and prints
the width of the flat bands. If no calculated file exists, the plotter uses
`reference/bands.dat`.
