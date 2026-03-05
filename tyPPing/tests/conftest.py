"""Shared fixtures for tyPPing tests."""

from pathlib import Path

import pytest

TYPPING_DIR = Path(__file__).parent.parent
TEST_EXAMPLE_DIR = TYPPING_DIR / "test_example"
INPUT_DATA_DIR = TYPPING_DIR / "tyPPing_input_data"


@pytest.fixture
def test_example_dir():
    return TEST_EXAMPLE_DIR


@pytest.fixture
def input_data_dir():
    return INPUT_DATA_DIR
