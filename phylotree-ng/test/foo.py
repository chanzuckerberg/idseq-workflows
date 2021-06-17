import os
import json

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
        sample["combined_contig_summary"] = os.path.join("s3://idseq-samples-development/tmorse-test/1/inputs/", name)
    else:
        sample["contig_fasta"] = os.path.join("s3://idseq-samples-development/tmorse-test/1/inputs/", name)
    samples[sample_name] = sample

accession_ids = ["NC_012532.1", "NC_035889.1"]

common_inputs = {
    "samples": list(samples.values()),
    "reference_taxon_id": 64320,
    "additional_reference_accession_ids": accession_ids,
    "superkingdom_name": "viruses"
}

print(json.dumps(common_inputs, indent=4))
