#!/bin/sh
# Draw the figure of every example again from its reference data and store it as
# <example>/reference/<name>.png. Each plot.py is run in a temporary copy of its folder that
# holds only the script and the reference data, so results of a run in the folder are not used.
set -e
here=$(cd "$(dirname "$0")" && pwd)
export MPLBACKEND=Agg
for example in "$here"/0*/; do
    work=$(mktemp -d)
    cp "$example/plot.py" "$work/"
    mkdir "$work/reference"
    cp "$example"/reference/*.dat "$work/reference/"
    (cd "$work" && python3 plot.py > /dev/null)
    cp "$work"/*.png "$example/reference/"
    echo "$(basename "$example"): $(cd "$work" && ls *.png)"
    rm -r "$work"
done
