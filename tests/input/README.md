# Input converter tests

Unit and round-trip tests for `tools/input/parser.py` and `tools/input/gendata.py`
(YAML/JSON to `Gendata.in`). The schema is generated from the source tree, so no
Fortran build is needed.

```sh
uv run python -m pytest tests/input
```

`TestRoundTrip` converts every `examples/**/Gendata*.in` to YAML and back and
compares the parsed values, which keeps the converter honest against the format
the examples use.
