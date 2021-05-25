import os
import json
from csv import DictReader

from Bio import SeqIO

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
        workflow_run_ids = [s["workflow_run_id"] for s in self.common_inputs["samples"]]
        sample_names = [s["sample_name"] for s in self.common_inputs["samples"]]

        res = self.run_miniwdl()
        outputs = res["outputs"]

        assert "phylotree.phylotree_newick" in outputs
        assert "phylotree.ska_distances" in outputs
        assert "phylotree.clustermap_png" in outputs
        assert "phylotree.clustermap_svg" in outputs
        assert "phylotree.variants" in outputs
        assert "phylotree.ncbi_metadata_json" in outputs
        assert len(outputs) == 6
            
        with open(outputs["phylotree.phylotree_newick"]) as f:
            tree_text = f.readlines()[0]
            for workflow_run_id in workflow_run_ids:
                assert str(workflow_run_id) in tree_text
            for accession_id in self.accession_ids:
                assert accession_id in tree_text
        
        identifiers = sorted(sample_names + self.accession_ids)
        with open(outputs["phylotree.ska_distances"]) as f:
            pairs = [sorted([r["Sample 1"], r["Sample 2"]]) for r in DictReader(f, delimiter="\t")]
            expected = [[a, b] for a in identifiers for b in identifiers if a < b]
            assert sorted(pairs) == expected

        with open(outputs["phylotree.variants"]) as f:
            assert identifiers == sorted([r.id for r in SeqIO.parse(f, "fasta")])

        with open(outputs["phylotree.ncbi_metadata_json"]) as f:
            m = json.load(f)
            assert sorted(list(m.keys())) == self.accession_ids
