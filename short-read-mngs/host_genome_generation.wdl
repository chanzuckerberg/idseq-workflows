version 1.0

task RunGenerateHostGenome {
  input {
    String docker_image_id
    String s3_wd_uri
    String host_name
    File host_fasta
    File ercc_fasta
    File ercc_gtf
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_genome_generation \
    --step-module idseq_dag.steps.generate_host_genome \
    --step-class PipelineStepGenerateHostGenome \
    --step-name generate_host_genome \
    --input-files '[["~{host_fasta}"]]' \
    --output-files '["fasta_with_ercc.fa", "gtf_with_ercc.fa", "~{host_name}_STAR_genome.tar", "~{host_name}_bowtie2_genome.tar"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-files '{ "ercc_fasta": "~{ercc_fasta}", "ercc_gtf": "~{ercc_gtf}" }' \
    --additional-attributes '{ "host_name": "~{host_name}", "max_star_part_size": null }'
  >>>
  output {
    File fasta_with_ercc = "fasta_with_ercc.fa"
    File gtf_with_ercc = "gtf_with_ercc.fa"
    File star_genome_tar = "~{host_name}_STAR_genome.tar"
    File bowtie2_genome_tar = "~{host_name}_bowtie2_genome.tar"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow idseq_host_filter {
  input {
    String docker_image_id
    String s3_wd_uri
    String host_name
    File host_fasta
    File ercc_fasta
    File ercc_gtf
  }

  call RunGenerateHostGenome {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      host_name = host_name
      host_fasta = host_fasta,
      ercc_fasta = ercc_fasta,
      ercc_gtf = ercc_gtf,
  }

  output {
    File fasta_with_ercc = RunGenerateHostGenome.fasta_with_ercc
    File gtf_with_ercc = RunGenerateHostGenome.gtf_with_ercc
    File star_genome_tar = RunGenerateHostGenome.star_genome_tar
    File bowtie2_genome_tar = RunGenerateHostGenome.bowtie2_genome_tar
  }
}
