import os
import csv
import json
from tempfile import NamedTemporaryFile


def test_RunIDSeqDedup_safe_csv(util, short_read_mngs_bench3_viral_outputs):
    # load the task's inputs from the end-to-end workflow test
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"], "call-host_filter/call-RunIDSeqDedup"
        )
    )

    input_files = []
    try:
        special_char_rows = 0
        ids = set()
        for in_f in inputs["priceseq_fa"]:
            f = NamedTemporaryFile("w")
            for line in open(in_f):
                if line[0] == ">" or line[0] == "@":
                    if special_char_rows < 10:
                        f.write(f"{line[0]}={line[1:]}")
                        special_char_rows += 1
                    else:
                        f.write(line)
                    ids.add(line[1:])
                else:
                    f.write(line)
            f.flush()
            input_files.append(f)

        assert special_char_rows == 10

        inputs["priceseq_fa"] = [f.name for f in input_files]

        outp = util.miniwdl_run(
            util.repo_dir() / "short-read-mngs/host_filter.wdl",
            "--task",
            "RunIDSeqDedup",
            "-i",
            json.dumps(inputs),
        )

        dups = outp["outputs"]["RunIDSeqDedup.duplicate_clusters_csv"]

        found_quotes = 0
        rows = 0
        # check we have an quotes before special characters space to prevent CSV injection
        with open(dups) as f:
            for row in csv.reader(f):
                rows += 1
                for elem in row:
                    if elem[0] == "'":
                        found_quotes += 1
                        continue
                    assert elem[0].isalnum(), f"cell starts with a special character '{elem}'"
        assert rows == len(ids) - 1
        assert found_quotes == 10
    finally:
        for f in input_files:
            f.close()
