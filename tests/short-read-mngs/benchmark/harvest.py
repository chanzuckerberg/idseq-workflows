#!/usr/bin/env python3
import sys
import argparse
import json
from pathlib import Path
from _util import load_benchmarks_yml

BENCHMARKS = load_benchmarks_yml()


def main():
    parser = argparse.ArgumentParser(
        sys.argv[0], description="harvest benchmark runs into importable JSON"
    )

    parser.add_argument(
        "outputs",
        metavar="sample=path/to/rundir",
        type=str,
        nargs="+",
        help="sample identifier & path to completed run directory",
    )
    parser.add_argument(
        "-o",
        metavar="FILENAME",
        dest="outfilename",
        type=str,
        required=True,
        help="output JSON file",
    )

    args = parser.parse_args(sys.argv[1:])
    harvest(**vars(args))


def harvest(outputs, outfilename):
    queue = []

    # process command line args
    i = 0
    while i < len(outputs):
        if outputs[i].endswith("="):
            assert i + 1 < len(outputs), f"missing path after {outputs[i]}"
            sample = outputs[i][:-1]
            rundir = Path(outputs[i + 1])
            i += 2
        else:
            parts = outputs[i].split("=")
            assert len(parts) == 2, f"invalid SAMPLE_ID=RUNDIR pair: {outputs[i]}"
            sample = parts[0]
            rundir = Path(parts[1])
            i += 1
        assert sample in BENCHMARKS["samples"], f"unknown sample {sample}"
        assert (rundir / "outputs.json").is_file(), f"couldn't find outputs.json in {rundir}"
        queue.append((sample, rundir))

    # harvest each supplied sample
    rslt = {}
    for sample, rundir in queue:
        assert sample not in rslt, f"repeated sample {sample}"
        with open(rundir / "outputs.json") as infile:
            outputs_json = json.load(infile)
        rslt[sample] = harvest_sample(sample, outputs_json)

    with open(outfilename, mode="w") as outfile:
        print(json.dumps(rslt, indent=2), file=outfile)


def harvest_sample(sample, outputs_json):
    ans = {}

    # collect read counts at various pipeline steps
    paired = "fastqs_1" in BENCHMARKS["samples"][sample]["inputs"]
    ans["input_read_pairs" if paired else "input_reads"] = read_output_jsonfile(
        outputs_json, "host_filter.input_read_count"
    )["fastqs"]
    for step in [
        "validate_input",
        "star",
        "trimmomatic",
        "priceseq",
        "idseq_dedup",
        "lzw",
        "bowtie2",
        "subsampled",
        "gsnap_filter",
    ]:
        ans[step + ("_out_read_pairs" if paired else "_out_reads")] = read_output_jsonfile(
            outputs_json, "host_filter." + step + "_out_count"
        )[step + "_out"]

    # TODO: count reads mapped to ERCC controls (target name ERCC-*)
    # TODO: count reads with only poor E-value hits to either NR or NT (quasi-unmapped)
    ans["initial_unmapped_reads"] = read_output_jsonfile(
        outputs_json, "non_host_alignment.annotated_out_count"
    )["unidentified_fasta"]
    ans["refined_unmapped_reads"] = read_output_jsonfile(
        outputs_json, "postprocess.refined_annotated_out_count"
    )["unidentified_fasta"]

    # collect NR/NT taxon counts
    ans = {"read_counts": ans, "taxon_counts": {}}
    ans["taxon_counts"]["NT"] = harvest_sample_taxon_counts(sample, outputs_json, "NT")
    ans["taxon_counts"]["NR"] = harvest_sample_taxon_counts(sample, outputs_json, "NR")
    return ans


def harvest_sample_taxon_counts(sample, outputs_json, dbtype):
    assert dbtype in ("NR", "NT")

    # read in the taxon counts & contig summary JSON files
    taxon_counts = read_output_jsonfile(
        outputs_json,
        "postprocess.refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json",
    )["pipeline_output"]["taxon_counts_attributes"]
    contig_summary = read_output_jsonfile(
        outputs_json, "postprocess.contig_summary_out_assembly_combined_contig_summary_json"
    )
    contig_summary = {(elt["count_type"], elt["taxid"]): elt for elt in contig_summary}

    # for each species in the taxon counts (excluding genus/family for now)
    ans = {}
    for rslt in taxon_counts:
        if rslt["count_type"] == dbtype and rslt["tax_level"] == 1 and int(rslt["tax_id"]) >= 0:
            # lookup corresponding contig summary
            contig_summary_key = (dbtype, rslt["tax_id"])
            if contig_summary_key in contig_summary:
                rslt_contigs = contig_summary[contig_summary_key]["contig_counts"]
            else:
                rslt_contigs = None
            assert rslt["tax_id"] not in ans

            # combine info
            ans[rslt["tax_id"]] = {
                "tax_level": rslt["tax_level"],
                "reads": rslt["nonunique_count"],
                "reads_dedup": rslt["unique_count"],
                "contigs": sum(1 for k in rslt_contigs if k != "*") if rslt_contigs else 0,
                "contigs_reads": sum(rslt_contigs[k] for k in rslt_contigs if k != "*")
                if rslt_contigs
                else 0,
                "avg_aln_len": rslt["alignment_length"],
                "avg_pct_id": rslt["percent_identity"],
                "avg_log10_E": rslt["e_value"],
            }
            # TODO: get taxon name?

    # sort by abundance, then E
    return {
        k: v
        for k, v in sorted(
            ans.items(), key=lambda kv: (0 - kv[1]["reads_dedup"], kv[1]["avg_log10_E"])
        )
    }


def read_output_jsonfile(outputs_json, key):
    """
    From the WDL outputs dict, read the contents of an output JSON file
    """
    filename = outputs_json["idseq_short_read_mngs." + key]
    with open(filename) as infile:
        return json.load(infile)


if __name__ == "__main__":
    main()
