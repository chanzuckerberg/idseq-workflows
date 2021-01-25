import os
import json
import csv
import tempfile


def test_RunAlignmentBlacklist(
    repo_dir, short_read_mngs_bench3_viral_outputs, miniwdl_inputs_outputs, miniwdl_run, RunFailed_stderr_msg
):
    task_name = "RunAlignment_rapsearch2_out"
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = miniwdl_inputs_outputs(
        os.path.join(short_read_mngs_bench3_viral_outputs["dir"], "call-non_host_alignment", f"call-{task_name}")
    )

    outp = miniwdl_run(
        os.path.join(repo_dir, "short-read-mngs/non_host_alignment.wdl"),
        "--task",
        task_name,
        "-i",
        json.dumps(inputs),
    )

    with open(os.path.join(outp["dir"], outp["outputs"][f"{task_name}.rapsearch2_hitsummary_tab"])) as f:
        taxids = set(row[2] for row in csv.reader(f, delimiter="\t"))

    assert "37124" in taxids, "taxid should be in hitsummary unless filtered out"
    assert "2169701" in taxids, "taxid should be in hitsummary unless filtered out"


    with tempfile.NamedTemporaryFile("w") as blacklist_file:
        blacklist_file.writelines(["37124", "2169701"])
        inputs["taxon_blacklist"] = blacklist_file.name

        outp = miniwdl_run(
            os.path.join(repo_dir, "short-read-mngs/non_host_alignment.wdl"),
            "--task",
            task_name,
            "-i",
            json.dumps(inputs),
        )

        hitsummary = os.path.join(outp["dir"], outp["outputs"][f"{task_name}.rapsearch2_hitsummary_tab"])
        deduped = os.path.join(outp["dir"], outp["outputs"][f"{task_name}.rapsearch2_deduped_m8"])

        with open(hitsummary) as f:
            taxids = set(row[2] for row in csv.reader(f, delimiter="\t"))

        assert "37124" not in taxids, "taxid should be filtered out"
        assert "2169701" not in taxids, "taxid should be filtered out"

        with open(hitsummary) as hf, open(deduped) as df:
            rows = zip(csv.reader(hf, delimiter="\t"), csv.reader(df, delimiter="\t"))
            assert all(hrow[0] == drow[0] for hrow, drow in rows), "hitsummary and deduped output should be aligned"

