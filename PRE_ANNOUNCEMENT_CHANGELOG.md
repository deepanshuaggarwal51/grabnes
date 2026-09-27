# GRABNES pre-announcement collaboration log

This document records repository changes made with Codex assistance while the
public GRABNES repository is being prepared for a broader announcement. The
repository is already publicly accessible, but the work described here should
be treated as pre-release consolidation until collaborators have reviewed and
committed it.

Last updated: 2026-09-27

## Current Git state

- The consolidation and initial examples were committed and pushed to `main`
  as commit `94964340c3322751f352e44337715912ca13505f`.
- The follow-up public update adds the `grabnes_testrun/` harness described
  below for collaborator testing in commits `e57b597` and `47cef60`.
- The next public update completes the four-example suite described below.
- The tag `pre-code-consolidation` points to commit `6133ab7`, the repository
  state before the solver directories were consolidated.
- Deleted files remain recoverable from that tag and from Git history.

## 1. Public examples scaffold

Added an `examples/` entry point for new users:

```text
examples/
├── README.md
├── 01_graphene_bands/
├── 02_graphene_dos/
├── 03_twisted_bilayer_bands/
└── 04_twisted_bilayer_dos/
```

Each directory now contains a minimal `Gendata.in`, location-independent
`run.sh`, Python/Matplotlib `plot.py`, focused README, and known-good reference
data. The launchers automatically find either the canonical executable or the
temporary verified `grabnes_testrun` executable and still accept an explicit
`GRABNES_BIN` override.

The examples are:

- pristine graphene bands along `K - Gamma - M - K'`;
- pristine graphene DOS from exact diagonalization on a `30 x 30` k-grid;
- bands of the 76-atom `(m,n)=(3,2)` commensurate twisted bilayer; and
- DOS of the same twisted bilayer on an `8 x 8` k-grid.

The small 13.17-degree twisted cell is intentional: it demonstrates the full
moire-cell workflow while remaining a laptop-scale test. The DOS examples use
deterministic exact diagonalization rather than the much larger stochastic
Kubo calculations from the 2023 tutorial. All four calculations completed with
`0 errors, 0 warnings`, and a second end-to-end run reproduced every reference
dataset numerically.

## 2. Solver directory consolidation

The repository previously contained three competing solver trees:

- `lanczosKuboCode/`
- `lanczosKuboCode_jiaqi/`
- `lanczosKuboCode_jinwoo/`

Investigation of the private repository history showed that the `jinwoo` tree
contains the active development line through September 2026, including recent
TAPW, Berry-link, Wilson-loop, orbital-moment, and semiclassical-orbit work.
The `jiaqi` source tree was also present byte-for-byte inside
`lanczosKuboCode_jinwoo/SrcJiaqi/`.

The latest `jinwoo` implementation was therefore promoted to the canonical
`lanczosKuboCode/` path. The `_jiaqi` and `_jinwoo` directories were removed.
The embedded duplicate `SrcJiaqi/` tree was also removed.

All repository documentation references were changed from
`lanczosKuboCode_jinwoo` to `lanczosKuboCode`, including:

- `Doxyfile`
- `DOCUMENTATION_SETUP.md`
- `SOC_IMPLEMENTATION_SUMMARY.md`
- `TAPW_CHERN_OPTIMIZATION_GUIDE.md`

The generated Chern-file header in `Src/diag.F90` now identifies the producer
as `GRABNES` instead of referring to the former personal directory name.

## 3. Repository cleanup

Removed material that should not be maintained in the public source tree:

- compiled executables and object, module, and archive files;
- tracked `build/` and `bin/` output;
- generated Sphinx documentation under `docs/_build/`;
- notebook checkpoints and Python bytecode;
- editor swap files and `.DS_Store`;
- backup and temporary Makefiles;
- personal `_prathap` source snapshots;
- the obsolete `diag.F90.org.reduck96` snapshot; and
- other generated sparse-matrix binaries and module files.

The three solver directories previously occupied approximately 81 MB in the
working tree. The consolidated canonical directory is approximately 6.6 MB.
This does not shrink existing Git history; it only cleans the checked-out tree
and future commits.

Expanded `lanczosKuboCode/.gitignore` so that local build products,
documentation output, notebook state, editor files, and backup files are not
accidentally committed again.

## 4. Build and documentation improvements

Replaced the placeholder solver README with current information covering:

- GRABNES's purpose;
- required compilers and numerical libraries;
- configuration and compilation;
- executable location and invocation;
- cleanup; and
- the link to the first public example.

Added `lanczosKuboCode/make.sys.example`, a portable GNU-oriented starting
configuration using `mpif90`, `gfortran`, OpenMP, ARPACK, LAPACK, and BLAS.
Site-specific configurations remain available under `lanczosKuboCode/Sys/`.

Made the nested clean rules tolerate partially created `MIO/` and `math/`
build directories that do not yet contain Makefiles.

Moved the `action` argument declaration in `Src/MIO/MPI/mpitime.F90` outside
the `TIMER` preprocessor guard. Without this change, GNU Fortran rejected the
subroutine when `TIMER` was not defined because `implicit none` left `action`
undeclared.

Trailing whitespace was removed mechanically from the promoted Fortran and
Makefile sources. No intended program logic was changed by that formatting
cleanup.

## 5. Verification performed

The following checks completed successfully:

- exactly one `lanczosKuboCode*` directory remains;
- no repository text references the removed `_jiaqi` or `_jinwoo` paths;
- no compiled objects, modules, archives, binaries, caches, swap files, or
  listed backup patterns remain in the canonical tree;
- `git diff --cached --check` reports no whitespace errors;
- the example shell and Python files pass syntax checks; and
- the graphene reference data match the Twistronics 2023 result exactly.

The consolidated source build was attempted with the local Homebrew GNU
Fortran and Open MPI installation. Compilation reached the generated MPI helper
program after the `mpitime.F90` fix, but the host linker failed with:

```text
ld: library 'crt1.o' not found
```

This was the result of the original build attempt with an obsolete Intel
Homebrew compiler. A later isolated native Apple Silicon build is documented
below.

## 6. Isolated native build and smoke-test harness

Added `grabnes_testrun/` temporarily inside the public checkout so
collaborators can reproduce the current Apple Silicon build while keeping all
compiler products under one ignored test subtree. Its local `.gitignore`
excludes `build/`, `bin/`, and generated example results from Git history.

The harness contains a source snapshot with the portability fixes discovered
during testing, an Apple Silicon `make.sys`, the graphene example, a detailed
`BUILD_REPORT.md`, and a top-level `smoke_test.sh`. Run it with:

```sh
cd grabnes_testrun
./smoke_test.sh
```

The test uses native Homebrew GCC 16.2, Open MPI 5.0.11, ARPACK 3.9.1_1, and
OpenBLAS 0.3.34. A clean checked build completed on this laptop, the program
reported `0 errors, 0 warnings`, and the generated graphene bands matched the
reference file byte-for-byte. Compiler warnings remain in the generated
legacy MPI wrappers and are documented in the build report.

The test-copy compatibility changes include modern logical operators, safe
handling of an unassociated single-process MPI list, corrected memory
accounting during deallocation, support for deeply nested source paths, and
selection of the neighbor routine intended for very small cells. These have
not yet been promoted to the canonical solver and should be reviewed first.

The harness is checkout-location independent: its scripts resolve their own
directories, the MIO generators accept long nested paths, and `make.sys`
discovers Homebrew library prefixes instead of hard-coding one installation
path. The required test `make.sys` is explicitly included, while compiler
products and generated results remain ignored.

## 7. Recommended collaborator review before announcement

1. Review the staged consolidation diff, especially the promoted solver source.
2. Build on a clean Linux environment with MPI, ARPACK, LAPACK, and BLAS.
3. Run `examples/01_graphene_bands/run.sh` and compare the resulting bands with
   `reference/bands.dat`.
4. Decide whether the two sparse-diagonalization test programs should be
   integrated into a maintained test suite or removed.
5. Add the planned graphene DOS and twisted-bilayer examples.
6. Add a repository license, `CITATION.cff`, software DOI, and manuscript
   references before the broader announcement.
7. Replace the current internal-group top-level README with a public-facing
   project landing page.

## Recovery

To inspect the repository before consolidation without changing the current
working tree:

```sh
git show pre-code-consolidation
git ls-tree -r --name-only pre-code-consolidation
```

Collaborators should avoid resetting the current work blindly. Review the
staged changes first, then commit them as one consolidation commit or split
them into examples, cleanup, and solver-promotion commits as appropriate.
