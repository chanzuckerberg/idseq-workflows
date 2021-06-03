import os
import json
import csv
import tempfile


def test_RunAssembly_defaults(util, short_read_mngs_bench3_viral_outputs):
    """
    On default settings, the assembly_out_contigs_fasta and assembly_out_contigs_all_fasta files
    should be identicaly (because all contigs pass the default min_contig_length filter)
    """
    task_name = "RunAlignment_gsnap_out"
    assembly_contigs_fasta = short_read_mngs_bench3_viral_outputs["outputs"][
        "idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"
    ]
    assembly_contigs_all_fasta = short_read_mngs_bench3_viral_outputs["outputs"][
        "idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_all_fasta"
    ]
    assert os.path.getsize(assembly_contigs_fasta) == os.path.getsize(assembly_contigs_all_fasta)


def test_RunAssembly_filtered(util, short_read_mngs_bench3_viral_outputs):
    """
    Re-run the assembly task with a higher min_contig_length and observe that contigs are removed.
    """
    # load the task's inputs from the end-to-end workflow test
    task_name = "RunAssembly"
    inputs, _ = util.miniwdl_inputs_outputs(
        os.path.join(
            short_read_mngs_bench3_viral_outputs["dir"],
            "call-postprocess",
            f"call-{task_name}",
        )
    )

    # override min_contig_length
    inputs["min_contig_length"] = 150

    # rerun
    outp = util.miniwdl_run(
        util.repo_dir() / "short-read-mngs/postprocess.wdl",
        "--task",
        task_name,
        "-i",
        json.dumps(inputs),
    )

    # count fasta
    assembly_contigs_fasta = os.path.join(
        outp["dir"], outp["outputs"][f"{task_name}.assembly_contigs_fasta"]
    )
    assembly_contigs_all_fasta = os.path.join(
        outp["dir"], outp["outputs"][f"{task_name}.assembly_contigs_all_fasta"]
    )
    assert count_fasta(assembly_contigs_fasta) and count_fasta(
        assembly_contigs_fasta
    ) < count_fasta(assembly_contigs_all_fasta)


def count_fasta(fn):
    ans = 0
    with open(fn) as infile:
        while True:
            line = infile.readline()
            if not line:
                return ans
            if line[0] == ">":
                ans += 1
