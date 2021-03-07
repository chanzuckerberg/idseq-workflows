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
        outputs = res["outputs"]
        with open(outputs["consensus_genome.compute_stats_out_output_stats"]) as fh:
            output_stats = json.load(fh)
        self.assertEqual(output_stats["sample_name"], "test_sample")
        self.assertGreater(output_stats["depth_avg"], 889)
        self.assertEqual(output_stats["total_reads"], 187444)
        self.assertEqual(output_stats["mapped_reads"], 187212)
        self.assertEqual(output_stats["mapped_paired"], 187151)
        self.assertEqual(output_stats["ercc_mapped_reads"], 0)
        self.assertEqual(output_stats["ref_snps"], 7)
        self.assertEqual(output_stats["ref_mnps"], 0)
        for output_name, output in outputs.items():
            if not isinstance(output, list):
                output = [output]
            for filename in output:
                with open(filename, "rb") as fh:
                    self.assertGreater(os.path.getsize(filename), 0)

    def test_zip_outputs(self):
        res = self.run_miniwdl(task="ZipOutputs", args=["prefix=test", "outputFiles=run.wdl"])
        with zipfile.ZipFile(res["outputs"]["ZipOutputs.output_zip"]) as fh:
            self.assertEqual(fh.namelist(), ["run.wdl"])
