import os
import sys
# import tempfile
# import json

from test_util import WDLTestCase

class TestPhylotree(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    samples_dir = os.path.join(os.path.dirname(__file__), "full_zika_test_data")

    samples = {}
    for i, name in enumerate(os.listdir(samples_dir)):
        first_dot = name.find(".")
        sample_name = name[:first_dot]
        sample = samples.get(sample_name, {
            "workflow_run_id": i,
            "sample_name": sample_name,
        })
        if "contig_summary" in name:
            sample["combined_contig_summary"] = os.path.join(samples_dir, name)
        else:
            sample["contig_fasta"] = os.path.join(samples_dir, name)
        samples[sample_name] = sample
        
    common_inputs = {
        # "samples": [
        #     {
        #         "sample_name": "test_sample1",
        #         "workflow_run_id": 12345,
        #         "contig_fasta": os.path.join(samples_dir, "test_sample1", "contigs.fasta"),
        #         "combined_contig_summary": os.path.join(samples_dir, "test_sample1", "combined_contig_summary.json"),
        #     },
        #     {
        #         "sample_name": "test_sample2",
        #         "workflow_run_id": 23456,
        #         "contig_fasta": os.path.join(samples_dir, "test_sample2", "contigs.fasta"),
        #         "combined_contig_summary": os.path.join(samples_dir, "test_sample2", "combined_contig_summary.json"),
        #     },
        # ],
        "samples": list(samples.values()),
        "reference": {"taxon_id": 64320},
        "additional_references": [
            {"accession_id": "NC_012532.1"},
            {"accession_id": "NC_035889.1"}
        ],
        "superkingdom_name": "viral"
    }

    def test_phylotree(self):
        res = self.run_miniwdl()
        self.assertEqual(res['outputs'], {})
