# scRNA-seq Differential Expression Analysis Benchmarking

## Overview
This repository provides a benchmarking framework for evaluating differential expression analysis methods applied to single-cell transcriptomics. It assesses the practical performance of diverse statistical approaches—including DESeq2, edgeR, limma, nebula and bootstrap-based methods. The goal is to offer empirical guidance for method selection in scRNA-seq studies.

## Nextflow Pipeline
This workflow accompanies our upcoming manuscript: **"Benchmarking at nominal FDR thresholds refines method selection in single-cell differential expression analysis"**. The provided Nextflow pipeline is designed to easily reproduce the analytical workflows, core results, and figures generated during our evaluations.

To ensure complete reproducibility, all software dependencies are encapsulated in Docker containers. Nextflow automatically retrieves these images from Docker Hub via Apptainer during execution, requiring no manual environment setup.

## Getting Started

### Prerequisites
To run the pipeline, the following tools are required:
- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)

### Running the Pipeline
First, clone the repository to your local machine:

```bash
git clone https://github.com/zw0124/sc-RNAseq-DEA-benchmarking.git
cd sc-RNAseq-DEA-benchmarking
```

Then, launch the pipeline from the project root:

```bash
nextflow run .
```

## Outputs
By default, the pipeline outputs results and generated figures to the `./output` directory. You can override this behavior by altering the `nextflow.config` file or passing the `--output` flag during execution.




