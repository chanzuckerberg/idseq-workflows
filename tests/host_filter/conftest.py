import sys
import os
import json
import logging
import pytest
import WDL


COMMON_INPUTS = {
    "docker_image_id": "docker.pkg.github.com/chanzuckerberg/idseq-workflows/idseq-main-public:17c3bd9",
    "dag_branch": "master",
    "s3_wd_uri": "s3://DUMMY_URI/",
}


@pytest.fixture(scope="session")
def miniwdl_run_cfg():
    "miniwdl run session-wide initialization"
    logging.basicConfig(level=logging.INFO + 1)
    with WDL._util.configure_logger():
        logger = logging.getLogger(__name__)
        cfg = WDL.runtime.config.Loader(logger)
        WDL.runtime.task.SwarmContainer.global_init(cfg, logger)
        return cfg


@pytest.fixture(scope="session")
def load_task():
    "function to load the WDL task object"
    return lambda wdl_path, task_name: next(
        task for task in WDL.load(wdl_path).tasks if task.name == task_name
    )


@pytest.fixture(scope="session")
def load_inputs_outputs():
    "function to load inputs.json and expected_outputs.json for a test case"

    def kappa(exe, call_dir):
        with open(os.path.join(call_dir, "inputs.json")) as infile:
            inputs = json.load(infile)
        for key in COMMON_INPUTS:
            if key in inputs:
                inputs[key] = COMMON_INPUTS[key]
        # absolutify paths to make cwd irrelevant
        inputs = WDL.values_from_json(inputs, exe.available_inputs, exe.required_inputs)
        inputs = WDL.Value.rewrite_env_files(
            inputs, lambda fn: os.path.join(call_dir, fn) if "://" not in fn else fn
        )
        inputs = WDL.values_to_json(inputs)
        with open(os.path.join(call_dir, "expected_outputs.json")) as infile:
            expected_outputs = WDL.values_from_json(
                json.load(infile), exe.effective_outputs, namespace=exe.name
            )
        expected_outputs = WDL.Value.rewrite_env_files(
            expected_outputs, lambda fn: os.path.join(call_dir, fn)
        )
        expected_outputs = WDL.values_to_json(expected_outputs, exe.name)
        return (inputs, expected_outputs)

    return kappa


@pytest.fixture
def miniwdl_run(miniwdl_run_cfg, tmpdir):
    "function to run the exe on given inputs (not necessarily the saved ones, e.g. to test error cases)"

    def kappa(exe, inputs, **kwargs):
        inputs = WDL.values_from_json(inputs, exe.available_inputs, exe.required_inputs)
        run_dir, outputs = WDL.runtime.run(
            miniwdl_run_cfg, exe, inputs, run_dir=str(tmpdir), **kwargs
        )
        return (run_dir, WDL.values_to_json(outputs, exe.name))

    return kappa


@pytest.fixture(scope="session")
def compare_outputs():
    "function to perform structured comparison of WDL outputs"

    def kappa(exe, actual, expected, file_digest=None, float_places=None):
        if file_digest is None:
            # default: represent each file just by its basename and size
            file_digest = (
                lambda fn: f"{os.path.basename(fn)}/{os.path.getsize(fn) if os.path.isfile(fn) else None}"
            )
        actual = WDL.values_from_json(actual, exe.effective_outputs, namespace=exe.name)
        actual = WDL.Value.rewrite_env_files(actual, file_digest)
        actual = _round_floats(WDL.values_to_json(actual, exe.name), float_places)
        expected = WDL.values_from_json(expected, exe.effective_outputs, namespace=exe.name)
        expected = WDL.Value.rewrite_env_files(expected, file_digest)
        expected = _round_floats(WDL.values_to_json(expected, exe.name), float_places)
        try:
            _compare_outputs(actual, expected)
        except AssertionError:
            print(json.dumps({"actual": actual, "expected": expected}, indent=2), file=sys.stderr)
            raise

    return kappa


def _round_floats(j, places):
    if places is not None:
        if isinstance(j, list):
            return [_round_floats(elt, places) for elt in j]
        if isinstance(j, dict):
            ans = {}
            for key in j:
                ans[key] = _round_floats(j[key], places)
            return ans
        if isinstance(j, float):
            return round(j, places)
    return j


def _compare_outputs(actual, expected):
    if isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected)
        for i in range(len(expected)):
            _compare_outputs(actual[i], expected[i])
    elif isinstance(expected, dict):
        assert isinstance(actual, dict) and not (set(expected) - set(actual))
        for key in expected:
            _compare_outputs(actual[key], expected[key])
    elif isinstance(expected, (str, int, float)):
        assert type(actual) is type(expected) and actual == expected
    elif expected is None:
        assert actual is None
    else:
        assert False