#!/usr/bin/env python3
import sys
import argparse
import json
from pathlib import Path
from _util import load_benchmarks_yml
from taxadb.taxid import TaxID

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
        "--taxadb", metavar="FILENAME", type=str, help="taxadb SQLite file, if available"
    )

    args = parser.parse_args(sys.argv[1:])
    harvest(**vars(args))


def harvest(outputs, taxadb):
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

    if taxadb:
        taxadb = TaxID(dbtype="sqlite", dbname=taxadb)

    # harvest each supplied sample
    rslt = {}
    for sample, rundir in queue:
        assert sample not in rslt, f"repeated sample {sample}"
        with open(rundir / "outputs.json") as infile:
            outputs_json = json.load(infile)
        rslt[sample] = harvest_sample(sample, outputs_json, taxadb)
        with open(rundir / "inputs.json") as infile:
            rslt[sample]["inputs"] = json.load(infile)

    print(json.dumps(rslt, indent=2))


def harvest_sample(sample, outputs_json, taxadb):
    ans = {}

    # collect read counts at various pipeline steps
    ans["paired"] = "fastqs_1" in BENCHMARKS["samples"][sample]["inputs"]
    ans["input_reads"] = read_output_jsonfile(outputs_json, "host_filter.input_read_count")[
        "fastqs"
    ]
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
        ans[step + "_reads"] = read_output_jsonfile(
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

    contig_summary = read_output_jsonfile(
        outputs_json, "postprocess.contig_summary_out_assembly_combined_contig_summary_json"
    )
    contig_summary = {(elt["count_type"], elt["taxid"]): elt for elt in contig_summary}
    contigs_mapped = set()
    for elt in contig_summary.values():
        for key in elt["contig_counts"]:
            if key != "*":
                contigs_mapped.add(key)

    contig_lengths = load_contig_lengths(outputs_json)
    ans.update(contigs_stats(contig_lengths))
    for key, value in contigs_stats(contig_lengths, contigs_mapped).items():
        ans["mapped_" + key] = value

    # collect NR/NT taxon counts
    ans = {"counts": ans, "taxa": {}}
    ans["taxa"]["NT"] = harvest_sample_taxon_counts(
        sample, outputs_json, contig_summary, contig_lengths, "NT", taxadb
    )
    ans["taxa"]["NR"] = harvest_sample_taxon_counts(
        sample, outputs_json, contig_summary, contig_lengths, "NR", taxadb
    )

    return ans


def harvest_sample_taxon_counts(
    sample, outputs_json, contig_summary, contig_lengths, dbtype, taxadb
):
    assert dbtype in ("NR", "NT")

    # read in the taxon counts & contig summary JSON files
    taxon_counts = read_output_jsonfile(
        outputs_json,
        "postprocess.refined_taxon_count_out_assembly_refined_taxon_counts_with_dcr_json",
    )["pipeline_output"]["taxon_counts_attributes"]

    # for each species in the taxon counts (excluding genus/family for now)
    ans = {}
    for rslt in taxon_counts:
        if rslt["count_type"] == dbtype and rslt["tax_level"] == 1 and int(rslt["tax_id"]) >= 0:
            # lookup corresponding contig summary
            contig_summary_key = (dbtype, rslt["tax_id"])
            if contig_summary_key in contig_summary:
                rslt_contigs = contig_summary[contig_summary_key]["contig_counts"]
            else:
                rslt_contigs = {}
            assert rslt["tax_id"] not in ans

            # combine info
            info = {
                "tax_level": rslt["tax_level"],
                "reads": rslt["nonunique_count"],
                "reads_dedup": rslt["unique_count"],
                "avg_aln_len": rslt["alignment_length"],
            }
            info.update(contigs_stats(contig_lengths, (k for k in rslt_contigs if k != "*")))
            info["contigs_reads"] = sum(rslt_contigs[k] for k in rslt_contigs if k != "*")
            if taxadb:
                info["tax_name"] = taxadb.sci_name(int(rslt["tax_id"]))
            ans[rslt["tax_id"]] = info

    # sort by abundance
    return {
        k: v
        for k, v in sorted(
            ans.items(), key=lambda kv: (kv[1]["reads_dedup"], kv[1]["avg_aln_len"]), reverse=True
        )
    }


def read_output_jsonfile(outputs_json, key):
    """
    From the WDL outputs dict, read the contents of an output JSON file
    """
    filename = outputs_json["idseq_short_read_mngs." + key]
    with open(filename) as infile:
        return json.load(infile)


def load_contig_lengths(outputs_json):
    """
    Generate dict contig id -> length
    """
    fasta = outputs_json["idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"]
    lengths = {}
    with open(fasta) as lines:
        cur = None
        for line in lines:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur[0]] = cur[1]
                cur = (line[1:].strip(), 0)
            else:
                cur = (cur[0], cur[1] + len(line.strip()))
        if cur is not None:
            lengths[cur[0]] = cur[1]
    return lengths


def contigs_stats(contig_lengths, ids=None):
    lengths = []
    for id in ids if ids else contig_lengths.keys():
        lengths.append(contig_lengths[id])
    lengths.sort(reverse=True)
    total_nt = sum(lengths)
    cumlen = 0
    N50 = 0
    for i, length_i in enumerate(lengths):
        if cumlen + length_i >= total_nt / 2:
            N50 = length_i
            break
        cumlen += length_i
    return {"contigs": len(lengths), "contigs_nt": total_nt, "contigs_N50": N50}


if __name__ == "__main__":
    main()
