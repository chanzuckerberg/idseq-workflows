# IDseq PhyloTree-NG workflow

version 1.1

task PrepareTaxonFasta {
    input {
        String docker_image_id
        String s3_wd_uri
        Array[SampleInfo] samples
    }
    command<<<
    set -euxo pipefail
    pass
    >>>
    output {
        Array[File] taxon_fastas = glob("taxon_fastas/*.fasta")
    }
    runtime {
        docker: docker_image_id
    }
}

task GeneratePhyloTree {
    input {
        String docker_image_id
        String s3_wd_uri
        String superkingdom_name
        Array[SampleInfo] samples
        Array[File] taxon_fastas
        Int taxid
        Array[Int] reference_taxids
        File nt_loc_db
        String nt_db
    }
    command<<<
    set -euxo pipefail
    echo input-files '[["~{sep('", "', taxon_fastas)}"], ["~{write_json(samples)}"]]'
    echo output-files '["phylo_tree.newick", "ncbi_metadata.json"]'
    echo additional-attributes '{"superkingdom_name": "~{superkingdom_name}", "taxid": ~{taxid}, "reference_taxids": [~{sep(", ", prefix("", reference_taxids))}]}'
    >>>
    output {
        File phylo_tree_newick = "phylo_tree.newick"
        File ncbi_metadata_json = "ncbi_metadata.json"
        Directory ska_outputs = "ska_outputs/"
    }
    runtime {
        docker: docker_image_id
    }
}

workflow phylotree-ng {
    input {
        String docker_image_id
        String s3_wd_uri
        String superkingdom_name
        Array[SampleInfo] samples
        Int taxid
        Array[Int] reference_taxids
    }

    call PrepareTaxonFasta {
        input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri,
        samples = samples,
    }

    call GeneratePhyloTree {
        input:
        docker_image_id = docker_image_id,
        s3_wd_uri = s3_wd_uri,
        superkingdom_name = superkingdom_name,
        samples = samples,
        taxid = taxid,
        reference_taxids = reference_taxids,
        taxon_fastas = PrepareTaxonFasta.taxon_fastas,
    }

    output {
        Array[File] taxon_fastas = PrepareTaxonFasta.taxon_fastas
        File phylo_tree_newick = GeneratePhyloTree.phylo_tree_newick
        File ncbi_metadata_json = GeneratePhyloTree.ncbi_metadata_json
        Directory ksnp3_outputs = GeneratePhyloTree.ksnp3_outputs
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
