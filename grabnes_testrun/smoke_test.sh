#!/bin/sh
set -eu

test_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
source_dir="$test_root/lanczosKuboCode"
example_dir="$test_root/examples/01_graphene_bands"

(
    cd "$source_dir"
    make clean
    make
)

rm -f "$example_dir/generate.bands" "$example_dir/job.out"
GRABNES_BIN="$source_dir/bin/grabnes" "$example_dir/run.sh"

cmp "$example_dir/generate.bands" "$example_dir/reference/bands.dat"
printf '%s\n' "PASS: clean build completed and graphene bands match the reference exactly."
