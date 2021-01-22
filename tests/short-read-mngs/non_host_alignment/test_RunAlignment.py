import os
import json


def test_RunValidateInput_invalid(
    repo_dir, short_read_mngs_bench3_viral_outputs, miniwdl_inputs_outputs, miniwdl_run, RunFailed_stderr_msg
):
    foo = os.listdir(os.path.join(short_read_mngs_bench3_viral_outputs["dir"], "call-non_host_alignment"))
    assert False, foo
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = miniwdl_inputs_outputs(
        os.path.join(short_read_mngs_bench3_viral_outputs["dir"], "call-non_host_alignment/call-RunAlignment")
    )

    # run the task with the manipulated inputs, expecting an error exit status
    outp = miniwdl_run(
        os.path.join(repo_dir, "short-read-mngs/host_filter.wdl"),
        "--task",
        "RunAlignment",
        "-i",
        json.dumps(inputs),
    )

