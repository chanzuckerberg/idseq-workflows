import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.exceptions import InsufficientReadsError
from idseq_dag.util.idseq_dedup_clusters import parse_clusters_file
from idseq_dag.util.count import save_duplicate_cluster_sizes


class PipelineStepRunIDSeqDedup(PipelineStep):  # Deliberately not PipelineCountingStep
    # TODO: add description
    """"""

    def validate_input_files(self):
        if not count.files_have_min_reads(self.input_files_local[0], 2):
            raise InsufficientReadsError("Insufficient reads before idseq-dedup")

    def run(self):
        input_fas = self.input_files_local[0]
        output_files = self.output_files_local()
        assert len(output_files) == len(input_fas) + 2, f"Context: {input_fas} -> {output_files}."
        output_fas = output_files[:len(input_fas)]
        duplicate_cluster_sizes_path = output_files[-1]
        assert duplicate_cluster_sizes_path.endswith(".tsv"), str(output_files)
        duplicate_clusters_path = output_files[-2]
        assert duplicate_clusters_path.endswith(".csv"), str(output_files)

        # See docstring above for explanation of these options.
        idseq_dedup_params = [
            '-i', input_fas[0], '-o', output_fas[0],
            '-l', '70',
            '-c', duplicate_clusters_path,
        ]
        if len(input_fas) == 2:
            idseq_dedup_params += ['-i', input_fas[1], '-o', output_fas[1]]
        command.execute(
            command_patterns.SingleCommand(
                cmd='idseq-dedup',
                args=idseq_dedup_params
            )
        )

        # Emit cluster sizes.  One line per cluster.  Format "<cluster_size> <cluster_read_id>".
        # This info is loaded in multiple subsequent steps using m8.load_duplicate_cluster_sizes,
        # and used to convert unique read counts to original read counts, and also to compute
        # per-taxon DCRs emitted alongside taxon_counts.
        clusters_dict = parse_clusters_file(duplicate_clusters_path)
        save_duplicate_cluster_sizes(duplicate_cluster_sizes_path, clusters_dict)

    def count_reads(self):
        self.should_count_reads = True
        # Here we intentionally count unique reads.
        self.counts_dict[self.name] = count.reads_in_group(
            self.output_files_local()[:-2])  # last two outputs are not fastas

