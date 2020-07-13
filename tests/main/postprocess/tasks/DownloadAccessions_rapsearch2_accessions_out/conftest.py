
import os
import pytest
import WDL


@pytest.fixture()
def exe(repo_dir, load_task):
    "load the task to be tested"
    return load_task(
        os.path.join(repo_dir, "main/postprocess.wdl"),
        "DownloadAccessions_rapsearch2_accessions_out"
    )
