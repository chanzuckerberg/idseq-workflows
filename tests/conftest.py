import pytest
from . import _util


@pytest.fixture(scope="session")
def util():
    return _util
