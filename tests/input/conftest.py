"""Fixtures for the converter tests: the schema of this source tree, once per session."""

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'tools', 'input')))

import gendata_schema  # noqa: E402


@pytest.fixture(scope='session')
def schema() -> gendata_schema.Schema:
    """The input schema of this source tree, generated once for the session."""
    return gendata_schema.generate()


@pytest.fixture(scope='session')
def schema_file(tmp_path_factory: pytest.TempPathFactory, schema: gendata_schema.Schema) -> str:
    """The same schema as a stored file, for the --schema path of the CLI."""
    path = tmp_path_factory.mktemp('schema') / 'schema.json'
    gendata_schema.save(schema, str(path))
    return str(path)
