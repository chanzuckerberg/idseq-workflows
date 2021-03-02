import os
import tempfile
import json
import zipfile
from unittest import TestCase
from subprocess import check_output

import yaml


class TestConsensusGenomes(TestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    with open(os.path.join(os.path.dirname(__file__), "local_test.yml")) as fh:
        common_inputs = yaml.safe_load(fh)

    def run_miniwdl(self, args, task=None, docker_image_id=os.environ["DOCKER_IMAGE_ID"]):
        cmd = ["miniwdl", "run", self.wdl] + args + [f"docker_image_id={docker_image_id}"]
        if task:
            cmd += ["--task", task]
        else:
            cmd += [f"{i}={v}" for i, v in self.common_inputs.items()]
        td = tempfile.TemporaryDirectory(prefix="/mnt/idseq-workflows-test-").name
        cmd += ["--verbose", "--error-json", "--debug", "--dir", td]
        print(cmd)
        res = check_output(cmd)
        return json.loads(res)

    def test_illumina_cg(self):
        fastqs_0 = os.path.join(os.path.dirname(__file__), "sample_sars-cov-2_paired_r1.fastq.gz")
        fastqs_1 = os.path.join(os.path.dirname(__file__), "sample_sars-cov-2_paired_r2.fastq.gz")
        args = ["sample=test_sample", f"fastqs_0={fastqs_0}", f"fastqs_1={fastqs_1}", "technology=Illumina"]
        res = self.run_miniwdl(args)
        print(res)

    def test_zip_outputs(self):
        res = self.run_miniwdl(task="ZipOutputs", args=["prefix=test", "outputFiles=run.wdl"])
        with zipfile.ZipFile(res["outputs"]["ZipOutputs.output_zip"]) as fh:
            self.assertEqual(fh.namelist(), ["run.wdl"])
