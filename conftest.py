import pytest
import test_util


@pytest.fixture(scope="session")
def util(tmpdir_factory):
    if not test_util._miniwdl_run_dir:
        test_util._miniwdl_run_dir = str(tmpdir_factory.getbasetemp() / "miniwdl_run")
    return test_util
