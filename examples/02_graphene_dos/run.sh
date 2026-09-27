#!/bin/sh
set -eu

example_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
default_bin="$example_dir/../../lanczosKuboCode/bin/grabnes"
test_bin="$example_dir/../../grabnes_testrun/lanczosKuboCode/bin/grabnes"
if [ -n "${GRABNES_BIN:-}" ]; then
    grabnes_bin=$GRABNES_BIN
elif [ -x "$default_bin" ]; then
    grabnes_bin=$default_bin
else
    grabnes_bin=$test_bin
fi

if [ ! -x "$grabnes_bin" ]; then
    printf '%s\n' "GRABNES executable not found or not executable: $grabnes_bin" >&2
    printf '%s\n' "Compile GRABNES (or run grabnes_testrun/smoke_test.sh) first." >&2
    printf '%s\n' "You can also set GRABNES_BIN=/path/to/grabnes." >&2
    exit 1
fi

cd "$example_dir"
rm -f generate.diag.DOS job.out file.13 file.14
"$grabnes_bin" Gendata.in > job.out

if [ ! -s generate.diag.DOS ]; then
    printf '%s\n' "GRABNES finished without creating generate.diag.DOS." >&2
    exit 1
fi

printf '%s\n' "Created $example_dir/generate.diag.DOS"
printf '%s\n' "Run 'python3 plot.py' to create graphene_dos.png."
