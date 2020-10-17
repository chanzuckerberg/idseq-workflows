#!/usr/bin/env python3
import sys
import os
import argparse
import subprocess
import multiprocessing
import concurrent.futures
import threading
import json
from pathlib import Path
from _util import load_benchmarks_yml

BENCHMARKS = load_benchmarks_yml()


def main():
    parser = argparse.ArgumentParser(sys.argv[0], description="run benchmark samples locally")
    parser.add_argument(
        "samples",
        metavar="SAMPLE",
        type=str,
        nargs="+",
        help="any of: " + ", ".join(BENCHMARKS["samples"].keys()),
    )
    parser.add_argument(
        "--docker-image-id", metavar="TAG", type=str, required=True, help="docker image tag"
    )
    parser.add_argument(
        "--dir",
        metavar="NEWDIR",
        type=str,
        required=True,
        help="directory under which to perform runs (mustn't already exist)",
    )
    parser.add_argument(
        "--references",
        metavar="ID",
        type=str,
        default="viral",
        help="reference database set; any of: " + ",".join(BENCHMARKS["references"].keys()),
    )
    parser.add_argument(
        "--verbose", action="store_true", default=False, help="run miniwidl --verbose"
    )

    args = parser.parse_args(sys.argv[1:])
    run_local(**vars(args))


def run_local(samples, docker_image_id, dir, references, verbose):
    """
    run the samples by invoking miniwdl locally on the current WDL revision
    """
    samples = set(samples)
    for sample_i in samples:
        assert sample_i in BENCHMARKS["samples"], f"unknown sample {sample_i}"

    dir = Path(dir)
    assert not dir.exists(), "--dir must't already exist"

    # parallelize runs on a thread pool; miniwdl itself may further restrict concurrency
    failure = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = {
            executor.submit(
                run_local_sample, sample_i, docker_image_id, dir, references, verbose
            ): sample_i
            for sample_i in samples
        }
        for future in concurrent.futures.as_completed(futures):
            exn = future.exception()
            if exn and not failure:
                # on first error, record it and `killall miniwdl` to quickly abort still-outstanding runs
                failure = (futures[future], exn)
                subprocess.run(["killall", "miniwdl"], check=False)

    if failure:
        (failed_sample, exn) = failure
        print(
            f"\nsample {failed_sample} failed; see {dir / failed_sample}/workflow.log\n",
            file=sys.stderr,
        )
        raise exn from None


_localize_lock = threading.Lock()


def run_local_sample(sample, docker_image_id, dir, references, verbose):
    local_driver_wdl = str(
        Path(__file__).parent.parent.parent / "short-read-mngs" / "local_driver.wdl"
    )

    # formulate input by merging selected dicts from BENCHMARK_YML
    wdl_input = {
        **BENCHMARKS["common_inputs"]["default"],
        **BENCHMARKS["references"][references],
        **BENCHMARKS["samples"][sample]["inputs"],
        **{"docker_image_id": docker_image_id},
    }

    # prime the download cache (serialize this to avoid redundant ops)
    with _localize_lock:
        subprocess.run(["miniwdl", "localize", local_driver_wdl, json.dumps(wdl_input)], check=True)

    # run the sample
    cmd = [
        "miniwdl",
        "run",
        local_driver_wdl,
        "--input",
        json.dumps(wdl_input),
        "--dir",
        str(dir / sample) + "/.",
        "--no-color",
    ]
    if verbose:
        cmd.append("--verbose")
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
