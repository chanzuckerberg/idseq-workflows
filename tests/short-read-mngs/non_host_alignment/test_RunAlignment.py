import os
import json


def test_RunValidateInput_invalid(
    repo_dir, short_read_mngs_bench3_viral_outputs, miniwdl_inputs_outputs, miniwdl_run, RunFailed_stderr_msg
):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = miniwdl_inputs_outputs(
        os.path.join(short_read_mngs_bench3_viral_outputs["dir"], "call-non_host_alignment/call-RunAlignment_rapsearch2_out")
    )

    # run the task with the manipulated inputs, expecting an error exit status
    outp = miniwdl_run(
        os.path.join(repo_dir, "short-read-mngs/non_host_alignment.wdl"),
        "--task",
        "RunAlignment_rapsearch2_out",
        "-i",
        json.dumps(inputs),
    )

    assert False, list(outp["outputs"].items())

