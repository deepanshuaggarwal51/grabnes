#!/bin/sh
set -eu

example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
default_bin="$example_dir/../../lanczosKuboCode/bin/grabnes"
grabnes_bin=${GRABNES_BIN:-$default_bin}

if [ ! -x "$grabnes_bin" ]; then
    printf '%s\n' "GRABNES executable not found or not executable: $grabnes_bin" >&2
    printf '%s\n' "Compile GRABNES first or set GRABNES_BIN=/path/to/grabnes." >&2
    exit 1
fi

cd "$example_dir"
if [ ! -s generate.xyz ] || [ ! -s layerIndex.dat ] || [ ! -s sublatticesSorted.dat ]; then
    python3 make_structure.py 6 5
fi
rm -f generate.spectral* job.out file.13 file.14
"$grabnes_bin" Gendata.in > job.out

if [ ! -s generate.spectral1A ]; then
    printf '%s\n' "GRABNES finished without creating generate.spectral1A." >&2
    exit 1
fi

printf '%s\n' "Created the generate.spectral* files in $example_dir"
printf '%s\n' "Run 'python3 plot.py' to create spectral_function.png."
