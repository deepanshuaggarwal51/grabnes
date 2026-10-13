#!/bin/sh
# Reproduce the reference data of the longer public examples (06 to 09): each takes minutes and
# example 07 needs 1.5 GB of memory, so they are not part of run_examples.sh.
#
#   GRABNES_BIN=/path/to/grabnes sh tests/regression/run_long_examples.sh [--work-dir DIR] [NAME ...]
#
# Every example is copied to the work directory and run there with its own run.sh; nothing is written
# into examples/. OMP_NUM_THREADS sets the threads of examples 06 and 07. Exit status 1 if a comparison fails.
set -eu

harness_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
examples_dir=$(CDPATH= cd -- "$harness_dir/../../examples" && pwd)
work_dir=""
names=""
while [ $# -gt 0 ]; do
    case $1 in
        --work-dir) work_dir=$2; shift ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        *) names="$names $1" ;;
    esac
    shift
done
[ -n "$names" ] || names="06_graphene_hbn_landau_levels 07_graphene_hbn_twist_series 08_magic_angle_twisted_bilayer 09_alternating_twist_trilayer"
[ -n "$work_dir" ] || work_dir=$(mktemp -d)
mkdir -p "$work_dir"
work_dir=$(CDPATH= cd -- "$work_dir" && pwd)
GRABNES_BIN=${GRABNES_BIN:-$examples_dir/../lanczosKuboCode/bin/grabnes}
case $GRABNES_BIN in /*) ;; *) GRABNES_BIN="$(pwd)/$GRABNES_BIN" ;; esac
export GRABNES_BIN

failures=0
for name in $names; do
    printf '\n%s\n' "== $name"
    run_dir="$work_dir/$name"
    rm -rf "$run_dir"
    cp -r "$examples_dir/$name" "$run_dir"
    if ! (cd "$run_dir" && ./run.sh > run.log 2>&1); then
        printf '%s\n' "  FAIL: run.sh ended with an error (see $run_dir/run.log)"
        failures=$((failures + 1)); continue
    fi
    case $name in
        06_*) (cd "$run_dir" && python3 plot.py --table result.dat > /dev/null)
              set -- --kind table --atol 1e-5 --rtol 0 "$run_dir/result.dat" "$run_dir/reference/dos_map.dat" ;;
        07_*) (cd "$run_dir" && python3 plot.py --table result.dat > /dev/null)
              set -- --kind table --atol 1e-5 --rtol 0 "$run_dir/result.dat" "$run_dir/reference/dos.dat" ;;
        *)    set -- --kind bands --atol 5e-6 --rtol 0 "$run_dir/generate.bands" "$run_dir/reference/bands.dat" ;;
    esac
    if python3 "$harness_dir/compare_output.py" "$@"; then
        printf '%s\n' "  PASS"
    else
        printf '%s\n' "  FAIL: the result does not reproduce the reference"
        failures=$((failures + 1))
    fi
done

printf '\n'
if [ "$failures" -eq 0 ]; then
    printf '%s\n' "PASS: the long examples reproduce their reference data (work directory: $work_dir)"
else
    printf '%s\n' "FAIL: $failures long example(s) failed (work directory: $work_dir)"
    exit 1
fi
