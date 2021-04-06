# Running Consensus-Genome locally

For Consensus Genome workflow we can follow a similar workflow to the `short-read-mngs` presented in wiki: [Running-WDL-workflows-locally](https://github.com/chanzuckerberg/idseq-workflows/wiki/Running-WDL-workflows-locally).

We first build a local Docker container image with the consensus genome workflow:

```bash
docker build -t idseq-consensus-genome idseq-workflows/consensus-genome
```

TIPS: For more detailed setup information
 - for miniwdl - https://github.com/openwdl/learn-wdl/blob/master/6_miniwdl_course/0_setup.md
 - for this example - https://github.com/openwdl/learn-wdl/blob/master/6_miniwdl_course/5c_cloud_spec_consensus-genome.md 

## Run 

We then use our local sample configuration file that points to IDseq's public references for smaller runs:

```bash
miniwdl run --verbose consensus-genome/run.wdl \
    docker_image_id=idseq-consensus-genome \
    fastqs_0=idseq-workflows/tests/consensus-genome/sample_sars-cov-2_paired_r1.fastq.gz \
    fastqs_1=idseq-workflows/tests/consensus-genome/sample_sars-cov-2_paired_r2.fastq.gz \
    sample=sample_sars-cov-2_paired \
    -i idseq-workflows/tests/consensus-genome/local_test.yml
```

Where:

* `docker_image_id=` should be set to the docker image tag you used when building the image (in our example, `idseq-consensus-genome`)
* `consensus-genome/run.wdl` is the WDL for the consensus genome sequencing workflow.
* `fastqs_0` and `fastqs_1` are the pair of FASTQ files. The ones referred to are small files to run locally.
* `sample` is the name to use where referencing the sample in the output files.
* `local_test.yml` supplies boilerplate workflow inputs, such as the S3 paths for the reference databases. For local run purposes, we use lighter references:
  * The human database for host removal only contains chromosome 1.
  * The kraken db used locally only has coronavirus sequences.

## Reference files
Filename | Provenance
---------|-----------
s3://idseq-public-references/consensus-genome/vadr-models-corona-1.1.3-1.tar.gz | Downloaded from https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/coronaviridae/1.2-1/vadr-models-corona-1.2-1.tar.gz on 2021-03-05
s3://idseq-public-references/consensus-genome/artic-primer-schemes.tar.gz | `primer_schemes` directory of https://github.com/artic-network/artic-ncov2019/commit/7e359dae37d894b40ae7e35c3582f14244ef4d36
`test/MT007544.fastq.gz` | Copied from https://github.com/artic-network/fieldbioinformatics/blob/master/test-data/MT007544/MT007544.fastq on 2021-03-06
## More Information

For more information, including a screencast of this example, see the `learn-miniwdl` open source course
- Screencast at - https://www.youtube.com/watch?v=bnXOoPm_F2I
- Miniwdl Course at - https://github.com/openwdl/learn-wdl/tree/master/6_miniwdl_course
- WDL course at - https://github.com/openwdl/learn-wdl
