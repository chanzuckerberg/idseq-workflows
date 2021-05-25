import os
import sys

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

    accession_ids = ["NC_012532.1", "NC_035889.1"]
        
    common_inputs = {
        "samples": list(samples.values()),
        "reference_taxon_id": 64320,
        "additional_reference_accession_ids": accession_ids,
        "superkingdom_name": "viruses"
    }

    def test_phylotree(self):
        res = self.run_miniwdl()
        outputs = res["outputs"]

        assert "phylotree.phylotree_newick" in outputs
        assert "phylotree.ska_distances" in outputs
        assert "phylotree.clustermap_png" in outputs
        assert "phylotree.clustermap_svg" in outputs
        assert "phylotree.variants" in outputs
        assert len(outputs) == 5
            
        with open(outputs["phylotree.phylotree_newick"]) as f:
            tree_text = f.readlines()[0]
            for sample in self.common_inputs["samples"]:
                assert str(sample["workflow_run_id"]) in tree_text
            for accession_id in self.accession_ids:
                assert accession_id in tree_text
