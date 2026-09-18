import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def tests_root_dir():
    """Returns the absolute path to the tests directory."""
    return Path(__file__).parent


@pytest.fixture(scope="session")
def test_data_dir(tests_root_dir):
    """Returns the absolute path to the tests directory."""
    return tests_root_dir / "test_data"
