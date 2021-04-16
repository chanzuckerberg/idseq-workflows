# clustered-phylotree
WDL workflow for split-kmer phylogeny with pre-filtering by hierarchical clustering

This repo contains the workflow used for prototyping a phylogenetic tree workflow that uses [SKA](https://github.com/simonrharris/SKA) split k-mers to create phylogenetic trees for an input set of samples. The pipeline works for complete genomes as well as raw sequences (currently supporting only .fasta inputs).  


**To run the clustered-phylotree pipeline locally...**

First, clone the repo

```
git clone git@github.com:katrinakalantar/clustered-phylotree.git
```

Then, build the docker image:

```
docker build -t clustphylo clustered-phylotree/
```

We then use the docker image to run the pipeline as follows:

```
miniwdl run --verbose clustered-phylotree/run.wdl docker_image_id=clustphylo data_directory=full_entero_data.tar.gz cut_height=.14 ska_align_p=.9
```

note: the pipeline requires as input the `data_directory`, which is a directory containing .fasta files (one file per sample) which is then tar zipped using the following command:

```
tar -czf full_entero_data.tar.gz full_entero_data
```

The optional parameters `cut_height` and `ska_align_p` allow you to specify the dendrogram cut height for pre-clustering and the ska alignment proportion parameters, respectively. The default options are shown in the command above.



The `/analysis_support` directory contains scripts that have been used to support the experimentation and validation of the phylotree pipeline.

