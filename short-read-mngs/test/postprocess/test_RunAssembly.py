import os
import json


def test_RunAssembly_defaults(util, short_read_mngs_bench3_viral_outputs):
    """
    On default settings, the assembly_out_contigs_fasta and assembly_out_contigs_all_fasta files
    should be identical (because all contigs pass the default min_contig_length filter)
    """
    assembly_contigs_fasta = short_read_mngs_bench3_viral_outputs["outputs"][
        "idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_fasta"
    ]
    assembly_contigs_all_fasta = short_read_mngs_bench3_viral_outputs["outputs"][
        "idseq_short_read_mngs.postprocess.assembly_out_assembly_contigs_all_fasta"
    ]
    assert fasta_headers(assembly_contigs_fasta) == fasta_headers(assembly_contigs_all_fasta)
    # quick&dirty: the contents should be identical modulo newlines
    with open(assembly_contigs_fasta) as in1, open(assembly_contigs_all_fasta) as in2:
        assert in1.read().replace("\n", "") == in2.read().replace("\n", "")


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
    contigs_headers = fasta_headers(assembly_contigs_fasta)
    contigs_all_headers = fasta_headers(assembly_contigs_all_fasta)
    assert contigs_headers and (set(contigs_all_headers) - set(contigs_headers))


def fasta_headers(fn):
    ans = []
    with open(fn) as infile:
        while True:
            line = infile.readline()
            if not line:
                return ans
            if line[0] == ">":
                ans.append(line)
