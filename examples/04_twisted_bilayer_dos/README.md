# 04: Twisted-bilayer graphene density of states

This example uses exact diagonalization on an `8 x 8` k-point grid to compute
the density of states of the same 76-atom `(m,n) = (3,2)` commensurate cell as
example 03. Gaussian broadening produces a compact, deterministic laptop test.
The Hamiltonian is the one of example 03 (nearest-neighbor intralayer hopping,
two-center interlayer hopping); see its README for what the reference data
represent. They were regenerated together with those of example 03.

## Expected result

![04: Twisted-bilayer graphene density of states](reference/twisted_bilayer_dos.png)

This is the figure `plot.py` draws from the reference data. After `./run.sh`, `python3 plot.py`
writes `twisted_bilayer_dos.png` in this folder from your own run; the two figures should look the same.

## Run

```sh
./run.sh
python3 plot.py
```

The calculation writes `generate.diag.DOS` and `job.out`; the plotting step
writes `twisted_bilayer_dos.png`. If no calculated file exists, the plotter
uses `reference/dos.dat`.

For a smoother production DOS, increase `KGrid` and
`NumberofEnergyPoints`, and reduce `Epsilon` after checking convergence.
