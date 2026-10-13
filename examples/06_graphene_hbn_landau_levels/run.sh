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
rm -rf m[0-9]*_B* job.out dos_map.dat
# The recursion method is threaded with OpenMP: OMP_NUM_THREADS sets the number of threads.
"$grabnes_bin" Gendata.in > job.out

if [ ! -s m00064_B26.0867/generate.DOS ]; then
    printf '%s\n' "GRABNES finished without the density of states of the last field (see job.out)." >&2
    exit 1
fi

printf '%s\n' "Created one directory per magnetic field, m<flux quanta>_B<field in T>/, in $example_dir"
printf '%s\n' "Run 'python3 plot.py' to create landau_levels.png."
