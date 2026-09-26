# GRABNES Fortran code

GRABNES (GRAphene and Boron Nitride Electronic Structure) is a Fortran code
for electronic-structure and quantum-transport calculations in layered
two-dimensional materials.

## Requirements

- GNU Fortran
- MPI with an `mpif90` compiler wrapper
- BLAS and LAPACK
- ARPACK
- GNU Make

## Build

Create a local build configuration from the portable template and compile:

```sh
cp make.sys.example make.sys
make
```

The executable is written to `bin/grabnes`. Site-specific configurations can
instead be copied from `Sys/` and adjusted for the local compiler and library
paths. The local `make.sys`, `build/`, and `bin/` paths are ignored by Git.

## Run

GRABNES accepts the input filename as its first argument:

```sh
bin/grabnes Gendata.in > job.out
```

For a minimal calculation with reference output, see
[`../examples/01_graphene_bands`](../examples/01_graphene_bands/).

## Clean

```sh
make clean
```

The code includes band-structure, density-of-states, Kubo transport, TAPW,
Berry-curvature, Chern-number, and semiclassical-orbit functionality. The
public examples document supported workflows as they are validated.
