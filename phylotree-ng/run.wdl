version 1.1

struct SampleInfo {
    String sample_name
    Int workflow_run_id
    File contig_fasta # assembly_out_assembly_contigs_fasta from mngs
    File combined_contig_summary # contig_summary_out_assembly_combined_contig_summary_json from mngs
}

workflow phylotree {
    input {
        Array[SampleInfo] samples
        Int reference_taxon_id

        # TODO: pass this to the relevant tasks and adjust SKA parameters as appropriate
        # (kSNP3 used a variable kmer length for each superkingdom: Viruses: 13, Bacteria: 19, Eukaryota: 19)
        # These kmer names can't be ported directly, because kSNP requires them to be odd while SKA requires a multiple of 3
        # Try flanking sequence length (-k) 12 for viruses, 18 for bacteria/eukaryotes for SKA
        # (SKA uses split kmers so the total resulting SKA kmer length here is (12*2 + 1) for 1 wobble base in the middle)
        String superkingdom_name # viruses, bacteria, or eukaryota

        # allow the user to pass specific reference taxids/accessions to include with the tree
        Array[Int] additional_reference_taxon_ids = []
        Array[String] additional_reference_accession_ids = []

        String cut_height = .16
        String ska_align_p = .9
        String docker_image_id
        String outgroup = ""
    }

    call GetSampleContigFastas {
        input:
        reference = reference,
        samples = samples,
        docker_image_id = docker_image_id
    }

    call GetReferenceFastas {
        input:
        references = additional_references,
        docker_image_id = docker_image_id
    }

    call RunSKA {
        input:
        sample_and_reference_fastas = flatten([GetSampleContigFastas.sample_contig_fastas, GetReferenceFastas.reference_fastas]),
        docker_image_id = docker_image_id
    }

    call ComputeClusters {
        input:
        ska_distances = RunSKA.distances,
        cut_height = cut_height,
        docker_image_id = docker_image_id
    }

    call GenerateClusterPhylos {
        input:
        clusters_directory = ComputeClusters.clusters_directory,
        ska_hashes = RunSKA.ska_hashes,
        ska_align_p = ska_align_p,
        docker_image_id = docker_image_id
    }

    call PlotClusterPhylos {
        input:
        ska_results = GenerateClusterPhylos.ska_results,
        outgroup = outgroup,
        docker_image_id = docker_image_id
    }

    output {
        File ska_hashes = RunSKA.ska_hashes
        File ska_distances = RunSKA.distances
        File stats_json = ComputeClusters.stats_json
        File dendrogram_png = ComputeClusters.dendrogram_png
        File clustermap_png = ComputeClusters.clustermap_png
        File clusters_directory = ComputeClusters.clusters_directory
        File dummy_subcluster_output = GenerateClusterPhylos.dummy_subcluster_output
        File ska_results = GenerateClusterPhylos.ska_results
        File dummy_plotphylos = PlotClusterPhylos.dummy_plotphylos
        File phylo_plot_outputs = PlotClusterPhylos.phylo_plot_outputs

        # TODO: These are the output names from the old phylotree. Check which of these we still need to emit
        # File phylo_tree_newick = GeneratePhyloTree.phylo_tree_newick
        # File ncbi_metadata_json = GeneratePhyloTree.ncbi_metadata_json
    }
}

task GetSampleContigFastas {
    # Given a list of samples, their workflow run IDs, and a reference taxon ID, retrieve contigs assigned to that taxon
    # ID or its descendants within each sample. Emit one fasta file per sample, to be graphed in a phylotree along with
    # references.

    # For each input sample, scan its contig summary to identify the contigs that need to be pulled from contigs.fasta.
    # TODO: For now, we do exact matching of the given reference taxon ID to the taxon IDs in the contig summary.
    #       This implies that the reference is selected at the species level, and that the contig summary is rolled up
    #       to species level as well. In the future we want to remove this assumption and select all contigs that map
    #       at all ranks under the given reference taxon.
    input {
        ReferenceInfo reference
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<
    python3 /bin/get_sample_contig_fastas.py --reference-taxid ~{reference_taxon_id} --samples "~{write_json(samples)}"
    >>>

    output {
        Array[File] sample_contig_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task GetReferenceFastas {
    # Given a list of reference taxon IDs or sequence accession IDs, retrieve the GenBank nt reference sequences
    # associated with them. Emit one fasta file per reference, to be graphed in a phylotree along with samples
    input {
        Array[ReferenceInfo] references
        String docker_image_id
    }

    command <<<
    for accession_id in $(jq -r .[].accession_id "~{write_json(references)}"); do
        taxoniq get-from-s3 --accession-id $accession_id > $accession_id.fasta
    done
    >>>

    output {
        Array[File] reference_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task RunSKA {
    input {
        Array[File] sample_and_reference_fastas
        String docker_image_id
    }

    command <<<
    for i in ~{sep=' ' sample_and_reference_fastas}; do
        ska fasta -o $(basename ${i%%.*}) $i
    done

    mkdir ska_hashes
    mv *.skf ska_hashes

    ska distance -o ska ska_hashes/*.skf

    tar -czf ska_hashes.tar.gz ska_hashes
    >>>

    output {
        File distances = "ska.distances.tsv"
        File ska_hashes = "ska_hashes.tar.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task ComputeClusters {
    # TODO: throw warning or error if sequences are too divergent
    input {
        File ska_distances
        String cut_height
        String docker_image_id
    }

    command <<<
    mkdir cluster_files
    python3 /bin/compute_clusters.py --ska-distances ~{ska_distances} --cut-height ~{cut_height}
    >&2 find .
    tar -czf clusters.tar.gz cluster_files
    >>>

    output {
        File stats_json = "stats.json"
        File clusters_directory = "clusters.tar.gz"
        File dendrogram_png = "dendrogram.png"
        File clustermap_png = "clustermap.png"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateClusterPhylos {
    input {
        File clusters_directory
        File ska_hashes
        String ska_align_p
        String docker_image_id
    }

    command <<<
    tar -xzvf "~{clusters_directory}"
    tar -xzvf "~{ska_hashes}"

    ls . >> output.txt
    echo "unzipped files" >> output.txt
    ls cluster_files/* >> output.txt
    ls ska_hashes/* >> output.txt

    CLUSTER_COUNTER=1
    mkdir ska_outputs

    for i in `ls cluster_files/*`
    do
       rm -r temp_cluster # remove if already existed
       mkdir temp_cluster

       for j in `cat $i`
       do
           cp ska_hashes/$j.skf temp_cluster
       done

       ska distance -o ska temp_cluster/*.skf
       ska merge -o ska.merged temp_cluster/*.skf
       ska distance -o ska ska_hashes/*.skf
       ska merge -o ska.merged ska_hashes/*.skf
       ska align -p "~{ska_align_p}" -o ska -v ska.merged.skf
       mv ska_variants.aln ska.variants.aln
       iqtree -s ska.variants.aln

       mkdir ska_outputs/cluster_$CLUSTER_COUNTER
       mv ska.* ska_outputs/cluster_$CLUSTER_COUNTER

       (( CLUSTER_COUNTER++ ))
    done

    tar -czf ska_outputs.tar.gz ska_outputs
    >>>

    output {
        File dummy_subcluster_output = "output.txt"
        File ska_results = "ska_outputs.tar.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task PlotClusterPhylos {
    input {
        File ska_results
        String outgroup
        String docker_image_id
    }

    command <<<
    tar -xzvf "~{ska_results}"

    ls > output.txt
    ls ska_outputs >> output.txt

    mkdir phylo_plot_outputs

    python3 /bin/plot_cluster_phylos.py --outgroup ~{outgroup}

    tar -czf phylo_plot_outputs.tar.gz phylo_plot_outputs
    >>>

    output {
        File dummy_plotphylos = "output.txt"
        File phylo_plot_outputs = "phylo_plot_outputs.tar.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task FetchSequenceByAccessionId {
    input {
        String accession_id
        String docker_image_id
    }

    command <<<
        taxoniq get-from-s3 --accession-id "~{accession_id}" > sequence.fa
        if [[ $? == 4 ]]; then
            export error=AccessionIdNotFound cause="Accession ID ~{accession_id} not found in the index"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 4
        fi
        exit $?
    >>>

    output {
        File sequence_fa = "sequence.fa"
    }

    runtime {
        docker: docker_image_id
    }
}
