## Introduction
Nanopore_assembly pipeline contains processes from Nanopore reads QC to assembly polishing.

short_reads consensus pipeline contains processes from pair-end reads QC to consensus sequences based on a reference provided by the user.

## Usage
To run the Nextflow pipelines with Docker (placing in the working directory the docker nextflow configuration file)
```
nextflow run Nanopore_assembly_docker.nf
nextflow run short_reads_consensus_docker.nf
```

To run the Nextflow pipelines with conda (placing in the working directory the conda nextflow configuration file)
```
nextflow run short_reads_consensus_conda.nf
```
