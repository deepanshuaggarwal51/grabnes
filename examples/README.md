# GRABNES examples

These examples are ordered from a fast installation check to more complete
moire calculations.

| Example | Purpose | Status |
| --- | --- | --- |
| [`01_graphene_bands`](01_graphene_bands/) | Pristine graphene band structure | Ready |
| `02_graphene_dos` | Pristine graphene density of states | Planned |
| `03_twisted_bilayer_bands` | Twisted-bilayer graphene band structure | Planned |

## Quick start

First compile GRABNES so that the executable is available at
`lanczosKuboCode/bin/grabnes`. Then run:

```sh
cd examples/01_graphene_bands
./run.sh
python3 plot.py
```

Set `GRABNES_BIN` if the executable is located elsewhere:

```sh
GRABNES_BIN=/path/to/grabnes ./run.sh
```

Each completed example contains its input, a short explanation, a plotting
script, and reference data from a known-good calculation.
