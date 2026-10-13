# Restructuring of `HamHopping` and `HamOnSite`

State and plan of the work on the branch `refactor/ham-structure`. The aim is
the convention of the original code: routines that fit on a few screens, each
with one task.

## How a change is verified

- The model survey (`tests/regression/model_survey.py`, about 430 cases) stores
  a fingerprint of the complete Hamiltonian of every case; a restructuring step
  must leave every fingerprint unchanged, with the checked gfortran build and
  with the Intel build.
- Coverage of the tests, measured with a `--coverage` build of gfortran on the
  survey and the physics checks (October 2026):

  | File | Executable lines run |
  | --- | --- |
  | `ham.F90` | 82 % (`HamHopping` 87 %, `HamOnSite` 89 %) |
  | `atoms.F90` | 74 % |
  | `neigh.F90` | 21 % (the legacy search routines are never run) |
  | `diag.F90` | 14 % (spectral functions, Chern numbers and `DiagH0TAPW_withBlockH` are not run at all) |

  The parts of `ham.F90` that no test reaches: `diagotBG`, the
  `FrankMagneticField` branch, `manyBubbles` / `bigBubble`, `TBG.DiagPRB`,
  `BNBN.Diag`. They are moved, where they are moved, without any change.
  **`diag.F90` needs tests before it can be restructured.**
- Four production inputs (graphene/hBN at zero field and in a field, a
  four-layer cell, a strained cell) are rerun and compared byte for byte.

## Step 1 (done): self-contained blocks become internal procedures

`tools/maintenance/outline_block.py` moves a block that is not inside a loop
into an internal procedure (`contains`) of the same routine. The statements are
not changed: an internal procedure uses the variables of its routine, so
nothing is passed, renamed or declared again. The tool refuses blocks with
`return`, `goto` or statement labels, blocks inside loops, and blocks inside an
OpenMP region that is opened outside them.

- `HamOnSite`: the body went from 2575 to about 490 lines; 15 procedures
  (`HamOnSiteSpeciesEnergies`, `HamOnSiteMoirePotential`,
  `HamOnSiteLayerPotential`, `HamOnSitePNP`, `HamOnSiteSquareFunction`,
  `HamOnSiteTrigonalCDW`, `HamOnSiteMoireCDW`, ...). What is left reads the
  switches and calls one procedure per term. Not moved: the `Bubbles` block
  (it contains a statement label).
- `HamHopping`: 11 procedures (`HamHopShellsFromRigid`,
  `HamHopF2G2sParameters`, `HamHopInterlayerParameters`, `HamHopSecondMoire`,
  `HamHopRandomStrain`, `HamHopRealStrain`, `HamHopPeriodicStrain`,
  `HamHopRealisticBubbles`, `HamHopMagneticField`, `HamHopHaldane`,
  `HamHopWriteTables`).

`HamOnSiteMoirePotential` and `HamHopRealStrain` were split later (step 3).

## Step 2 (done): the main loop of `HamHopping`

The loop over atoms and neighbors was one OpenMP loop of 4500 lines in which the model was decided again for
every bond. It is now a loop of about 400 lines that dispatches to one procedure per model.

1. **Constants out of the loop.** 192 assignments of literal constants of the Srivani model stood inside the
   parallel loop, to shared variables (every thread wrote the same values for every bond). They are set once
   in `HamHopSrivaniConstants`; `aCC` and `aBN1` likewise. The loop no longer assigns any shared local
   variable except two print-once flags.
2. **A defect found on the way.** Six switches of the Srivani model were read for `TypeOfBL BLSrivani` only
   but tested for `BLKaxiras` too; in a Kaxiras run they were undefined and the result depended on the
   compiler and on the build. They are read for both now, and the 22 Kaxiras survey cases have left the list
   of compiler-dependent cases.
3. **One procedure per model** (`tools/maintenance/outline_branch.py`). A branch of the loop is moved,
   unchanged, into an internal procedure. The variables that have a copy per thread (PRIVATE, FIRSTPRIVATE,
   REDUCTION, loop indices) are its arguments; the model parameters and the arrays of the Hamiltonian are
   used from `HamHopping`. A `cycle` becomes the result of a logical function. Temporaries that only one model
   uses are local to its procedure and start as "not set" for every bond, so a value left from another bond
   cannot be used by mistake.

   | Part of the loop | Procedures |
   | --- | --- |
   | bond terms of the moire models, per atom | `HamHopBondTermsMoire`, `HamHopBondTermsTBG`, `HamHopBondTermsGBN` |
   | intralayer hopping | `HamHopIntralayerGraphene` (TypeOfBL Koshino, Mayou, HTC: F2G2, Koshino or Mayou intralayer models), `HamHopIntralayerGBN`, `HamHopIntralayerBNBN`, `HamHopIntralayerEncapsulated3` ... `7`, `HamHopIntralayerT3GwithBN`, `HamHopIntralayerBNt2GBN`, `HamHopIntralayerT2GBN`, `HamHopIntralayerKaxirasSrivani`, `HamHopIntralayerSrivani`, `HamHopIntralayerJeil`, `HamHopIntralayerDefault` |
   | interlayer hopping | `HamHopInterlayerTwoCenter` (TypeOfBL Koshino, HTC; with `KoshinoSR` and `couplingFactor` this is the model for twisted graphene layers), `HamHopInterlayerGBN`, `HamHopInterlayerBNBN`, `HamHopInterlayerMayou`, `HamHopInterlayerKaxirasSrivani`, `HamHopInterlayerJeil` |

   The body of `HamHopping` went from 7390 to about 1900 lines, most of which is now the reading of the
   switches and parameters.
4. **Declarations** (`tools/maintenance/tidy_declarations.py`): 215 of the 917 local variables of
   `HamHopping` were no longer used and are removed; the others are grouped under headings (switches,
   model parameters, counters, per-thread variables, other).

## Step 3, first part (done)

- **`HamOnSiteMoirePotential`**: from 1310 to about 530 lines. The loop over atoms of the structure-based stacks
  calls one procedure per system (`HamOnSiteAtomHarmonicTBG`, `HamOnSiteAtomHarmonicBNBN`,
  `HamOnSiteAtomHarmonicGBN`); the effective models are `HamOnSiteMoireEffective` (graphene/hBN, one moire),
  `HamOnSiteMoireEffectiveBilayer`, `HamOnSiteSecondMoireTerms` and `HamOnSiteMoireKekule`. Both tools take
  `--host HOST` for a block that lies in an internal procedure: the new procedure becomes another internal
  procedure of the host (an internal procedure cannot contain procedures).
- **`HamHopRealStrain`**: from 316 to 118 lines. The rule that gives the unstrained length of a bond from its
  neighbor shell stood 19 times in the loop (once per stack and kind of layer); the 19 copies were compared
  with one template by a script, found identical, and replaced by calls of one module function,
  `HamStrainReferenceDistance(d, a1, a, a1Ref, aRef)`.
- The remaining long procedure is `HamOnSiteAtomHarmonicGBN` (366 lines: graphene/hBN, the encapsulated and
  the mixed stacks in one chain).

## Step 3, remaining

- **Which model for which layers.** F2G2 is the more accurate intralayer model for untwisted layers; the
  two-center interlayer model with `KoshinoSR` and `couplingFactor` is the one to use whenever two graphene
  layers are twisted. Today this is spread over `middleTwist`, `forceBilayerF2G2Intralayer` and the stack
  switches. Aim: a table built once, per layer or pair of layers, printed in the output of the run.
- **One loop per model.** The procedures are called per bond from the common loop. Giving each its own loop
  removes the dispatch from the loop and lets the model be chosen once.
- **Parameter types.** The procedures still use the parameters of `HamHopping` by host association. A
  derived type per model family, filled by the parameter procedures, turns them into module procedures with
  short argument lists that can be tested alone.
- **Duplicates.** `HamHopIntralayerEncapsulated3` to `7` are near copies of about 180 lines each; the
  `...GBN`, `...T3GwithBN`, `...BNt2GBN`, `...T2GBN` procedures share most of their text with them.

## Model of each pair of layers chosen from the structure (opt-in)

- `Stack.MeasureTwist .true.`: the orientation of every layer and the twist of every interface are measured
  from the bond directions and printed. Nothing else changes.
- `Stack.ModelFromTwist .true.`: every untwisted AB (Bernal) pair of graphene layers takes the bilayer F2G2
  model of Jung and MacDonald (PRB 89, 035405): site energies, intralayer and interlayer terms, with the dimer
  sites found from the structure. Every other graphene layer keeps the single-layer F2G2 parameters and every
  other interface the model of `TypeOfBL` (the two-center model for twisted layers). AA pairs keep the
  two-center model. One nearest-neighbor hopping (`TB.Hopping`) is used in all layers.
- `Stack.AlignModelEnergies .true.` (needs the switch above): the two parameter sets put the Dirac point at
  different energies, because their second- and fifth-shell hoppings differ (energy of a sublattice at K:
  site energy - 3 t2 + 6 t5 - 3 t6). Without this switch the run prints the energies and a warning when they
  differ by more than 5 meV. With it the site energy of the pair sites without an atom on top is
  `TB.OnSiteC - 3 (t2 - t2A) + 6 (t5 - t5A)`, so the neutrality point of an isolated AB bilayer and the Dirac
  point of an isolated single layer coincide; the dimer sites keep their distance `TB.OnSiteC2 - TB.OnSiteC1`
  (the paper's 0.0211 eV when neither key is given). `TB.OnSiteC` alone then sets the reference energy.
  What it does not do: it is the alignment of isolated layers (no charge transfer, no difference of work
  function between a single layer and a bilayer), and it does not make the velocities of the two fits
  consistent (the bilayer fit was made with a nearest-neighbor hopping of 3.01 eV).
- Tests: survey cases `xyz4_tdbg*` and `xyz3_mixed_from_structure*`; physics checks `twist_measure` and
  `model_alignment` (levels at the folded K points of a decoupled single layer + AB pair, to 2e-6 eV).

## Not planned here

- `diag.F90`: tests first (TAPW with spin-orbit terms, spectral functions,
  Chern numbers have none).
- `AtomsPos` (1192 lines) and `NeighList` (504 lines) can follow step 1 as they
  are; the legacy search routines of `neigh.F90` are candidates for removal
  rather than restructuring.
