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
        "-o",
        metavar="FILENAME",
        dest="outfilename",
        type=str,
        required=True,
        help="output JSON file",
    )
    parser.add_argument(
        "--taxadb", metavar="FILENAME", type=str, help="taxadb SQLite file, if available"
    )

    args = parser.parse_args(sys.argv[1:])
    harvest(**vars(args))


def harvest(outputs, outfilename, taxadb):
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

    with open(outfilename, mode="w") as outfile:
        print(json.dumps(rslt, indent=2), file=outfile)


def harvest_sample(sample, outputs_json, taxadb):
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

    contig_summary = read_output_jsonfile(
        outputs_json, "postprocess.contig_summary_out_assembly_combined_contig_summary_json"
    )
    contig_summary = {(elt["count_type"], elt["taxid"]): elt for elt in contig_summary}
    contigs_mapped = set()
    for elt in contig_summary.values():
        for key in elt["contig_counts"]:
            if key != "*":
                contigs_mapped.add(key)

    ans.update(contigs_stats(outputs_json))
    for key, value in contigs_stats(outputs_json, contigs_mapped).items():
        ans["mapped_" + key] = value

    # collect NR/NT taxon counts
    ans = {"read_counts": ans, "taxon_counts": {}}
    ans["taxon_counts"]["NT"] = harvest_sample_taxon_counts(
        sample, outputs_json, contig_summary, "NT", taxadb
    )
    ans["taxon_counts"]["NR"] = harvest_sample_taxon_counts(
        sample, outputs_json, contig_summary, "NR", taxadb
    )
    return ans


def harvest_sample_taxon_counts(sample, outputs_json, contig_summary, dbtype, taxadb):
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
                rslt_contigs = None
            assert rslt["tax_id"] not in ans

            # combine info
            info = {
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
            if taxadb:
                info["tax_name"] = taxadb.sci_name(int(rslt["tax_id"]))
            ans[rslt["tax_id"]] = info

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


def contigs_stats(outputs_json, ids=None):
    fasta = outputs_json["idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"]
    lengths = []
    with open(fasta) as lines:
        cur = None
        for line in lines:
            if line.startswith(">"):
                if cur is not None and (not ids or line[1:].strip() in ids):
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur is not None:
            lengths.append(cur)
    lengths.sort(reverse=True)
    total_nt = sum(lengths)
    cur = 0
    for i, length_i in enumerate(lengths):
        if cur + length_i > total_nt / 2:
            N50 = length_i
            break
        cur += length_i
    assert N50
    return {"contigs_count": len(lengths), "contigs_total_nt": total_nt, "contigs_N50": N50}


if __name__ == "__main__":
    main()
