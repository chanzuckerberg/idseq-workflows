import os
import json
import csv


def test_RunValidateInput_invalid(
    repo_dir, short_read_mngs_bench3_viral_outputs, miniwdl_inputs_outputs, miniwdl_run, RunFailed_stderr_msg
):
    task_name = "RunAlignment_rapsearch2_out"
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = miniwdl_inputs_outputs(
        os.path.join(short_read_mngs_bench3_viral_outputs["dir"], "call-non_host_alignment", f"call-{task_name}")
    )

    # run the task with the manipulated inputs, expecting an error exit status
    outp = miniwdl_run(
        os.path.join(repo_dir, "short-read-mngs/non_host_alignment.wdl"),
        "--task",
        task_name,
        "-i",
        json.dumps(inputs),
    )

    with open(os.path.join(outp["dir"], f"{task_name}.rapsearch2_m8")) as f:
        taxids = [row[2] for row in csv.reader(f)]
    assert False, ", ".join(taxids)

