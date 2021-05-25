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
        reference_taxon_id = reference_taxon_id,
        samples = samples,
        docker_image_id = docker_image_id
    }

    call GetReferenceAccessionFastas {
        input:
        accession_ids = additional_reference_accession_ids,
        docker_image_id = docker_image_id
    }

    call GetReferenceTaxonFastas {
        input:
        taxon_ids = additional_reference_taxon_ids,
        docker_image_id = docker_image_id
    }

    call RunSKA {
        input:
        sample_and_reference_fastas = flatten([GetSampleContigFastas.sample_contig_fastas, GetReferenceAccessionFastas.reference_fastas, GetReferenceTaxonFastas.reference_fastas]),
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

    output {
        File ska_hashes = RunSKA.ska_hashes
        File ska_distances = RunSKA.distances
        File stats_json = ComputeClusters.stats_json
        File clustermap_png = ComputeClusters.clustermap_png
        File clustermap_svg = ComputeClusters.clustermap_svg
        File treefile = GenerateClusterPhylos.treefile
        File distances = GenerateClusterPhylos.distances
        # TODO: These are the output names from the old phylotree. Check which of these we still need to emit
        # File phylo_tree_newick = GeneratePhyloTree.phylo_tree_newick
        # File ncbi_metadata_json = GeneratePhyloTree.ncbi_metadata_json
    }
}

task GetReferenceAccessionFastas {
    input {
        Array[String] accession_ids
        String docker_image_id
    }

    command <<<
        for accession_id in ~{sep=' ' accession_ids}; do
            taxoniq get-from-s3 --accession-id $accession_id > $accession_id.fasta
            if [[ $? == 4 ]]; then
                export error=AccessionIdNotFound cause="Accession ID $accession_id not found in the index"
                jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
                exit 4
            fi
        done
    >>>

    output {
        Array[File] reference_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task GetReferenceTaxonFastas {
    input {
        Array[Int] taxon_ids
        String docker_image_id
    }

    command <<<
        for taxon_id in ~{sep=' ' taxon_ids}; do
            taxoniq get-from-s3 --taxon-id $taxon_id > $taxon_id.fasta
            if [[ $? == 4 ]]; then
                export error=TaxonIdNotFound cause="Taxon ID $taxon_id not found in the index"
                jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
                exit 4
            fi
        done
    >>>

    output {
        Array[File] reference_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
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
        Int reference_taxon_id
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
        File clustermap_png = "clustermap.png"
        File clustermap_svg = "clustermap.svg"
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

    N_CLUSTERS=$(ls cluster_files | wc -l)
    
    if [[ $N_CLUSTERS != 1 ]]; then
        export error=MultipleClusters cause="SKA found multiple clusters ($N_CLUSTERS)"
        jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
        exit 1
    fi

    CLUSTER_FILE=$(ls cluster_files | head -n 1)

    mkdir cluster

    for hash in `cat $CLUSTER_FILE`
    do
        cp ska_hashes/$hash.skf cluster
    done

    mkdir ska_outputs

    ska distance -o ska cluster/*.skf
    ska merge -o ska.merged cluster/*.skf
    ska distance -o ska ska_hashes/*.skf
    ska merge -o ska.merged ska_hashes/*.skf
    ska align -p "~{ska_align_p}" -o ska -v ska.merged.skf
    mv ska_variants.aln ska.variants.aln
    iqtree -s ska.variants.aln
    mv ska.variants.aln.treefile tree.nwk
    >>>

    output {
        File treefile = "tree.nwk"
        File distances = "ska.distances.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}

