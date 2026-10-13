#!/bin/sh
set -eu

example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
default_bin="$example_dir/../../lanczosKuboCode/bin/grabnes"
grabnes_bin=${GRABNES_BIN:-$default_bin}
case $grabnes_bin in /*) ;; *) grabnes_bin="$(pwd)/$grabnes_bin" ;; esac

if [ ! -x "$grabnes_bin" ]; then
    printf '%s\n' "GRABNES executable not found or not executable: $grabnes_bin" >&2
    printf '%s\n' "Compile GRABNES first or set GRABNES_BIN=/path/to/grabnes." >&2
    exit 1
fi

cd "$example_dir"
# One run per twist angle, each in its own directory angle_<degrees>/.
# The recursion method is threaded with OpenMP: OMP_NUM_THREADS sets the number of threads.
for input in Gendata_*.in; do
    angle=${input#Gendata_}; angle=${angle%.in}
    rm -rf "angle_$angle"; mkdir "angle_$angle"
    cp "$input" "angle_$angle/Gendata.in"
    (cd "angle_$angle" && "$grabnes_bin" Gendata.in > job.out)
    if [ ! -s "angle_$angle/generate.DOS" ]; then
        printf '%s\n' "GRABNES finished without creating angle_$angle/generate.DOS." >&2
        exit 1
    fi
    printf '%s\n' "Created $example_dir/angle_$angle/generate.DOS"
done
printf '%s\n' "Run 'python3 plot.py' to create twist_series_dos.png."
