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
        "samples": list(samples.values()),
        "reference_taxon_id": 64320,
        "additional_reference_accession_ids": [
            "NC_012532.1",
            "NC_035889.1",
        ],
        "superkingdom_name": "viral"
    }

    def test_phylotree(self):
        res = self.run_miniwdl()
        a = res["outputs"]["phylotree.foo"]
        with open(a) as f:
            for line in f:
                sys.stderr.write(line)
        self.assertEqual(res['outputs'], {})
