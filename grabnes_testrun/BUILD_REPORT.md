# GRABNES isolated build report

Tested on 2026-09-26 on an Apple Silicon (arm64) MacBook Pro running
macOS 26.5.2. Build products are contained under this test directory and are
excluded from Git by the local `.gitignore`.

## Result

- Native arm64 compilation: **PASS**
- `examples/01_graphene_bands` execution: **PASS**
- Program summary: `0 errors, 0 warnings`
- Output comparison: **PASS**, byte-for-byte identical to
  `reference/bands.dat`
- SHA-256 of both files:
  `312203821daaa76b6f9891988b3a4d59f78e1e7cfc0dc0011a262d394e8c9df8`

Run the complete clean-build smoke test with:

```sh
cd /path/to/grabnes/grabnes_testrun
./smoke_test.sh
```

## Local toolchain

The test uses native Homebrew installations of GCC/GFortran, Open MPI, ARPACK,
and OpenBLAS. `lanczosKuboCode/make.sys` resolves Homebrew library prefixes at
build time and obtains the compiler commands from `PATH`; no Homebrew prefix or
checkout location is hard-coded.

## Compatibility changes made only in this test copy

The source snapshot required several updates for current GFortran. None of
these changes were applied to the public repository during this test.

1. Legacy MPI argument mismatches are accepted with
   `-fallow-argument-mismatch`.
2. Old logical comparisons using `.eq.` were changed to `.eqv.` in
   `Src/diag.F90` and `Src/calc.F90`.
3. An unreachable, self-aliasing `move_alloc` block was removed from
   `Src/diag.F90`.
4. An unassociated one-process MPI receive-list pointer is guarded before
   calling `size()` in `Src/ham.F90`.
5. Memory accounting now obtains an allocation's size before deallocating it
   in `Src/MIO/mem_temp.f90`.
6. The two-atom example selects the small-cell neighbor routine with
   `nonBulkSmall .true.` and disables optional verbose data files. The latter
   avoids the legacy MIO formatted-record limit of 80 characters.
7. OpenMP directives remain disabled for this conservative test build, while
   `libgomp` is linked because the source calls OpenMP query functions.
8. The legacy MIO code generators use larger filename buffers so the harness
   also works from a deeply nested checkout path.

The checked build deliberately uses `-O0 -g -fcheck=all -fbacktrace` so bounds,
pointer, and other runtime errors are detected during the smoke test.

## Remaining modernization note

The generated legacy MPI wrappers compile with argument rank/type warnings.
They did not prevent the one-process example from producing the exact reference
result, but should be modernized before treating a warning-free build as a
release requirement.
