# Mapping the Spatial transcriptomic signature of the hippocampus during memory consolidation.

## Abstract

Memory consolidation involves discrete patterns of transcriptional events in the hippocampus. Despite the emergence of single-cell transcriptomic profiling techniques, defining learning-responsive gene expression across subregions of the hippocampus has remained challenging. Here, we utilized unbiased spatial sequencing to elucidate transcriptome-wide changes in gene expression in the hippocampus following learning, enabling us to define molecular signatures unique to each hippocampal subregion. We find that each subregion of the hippocampus exhibits distinct yet overlapping transcriptomic signatures. Although the CA1 region exhibited increased expression of genes related to transcriptional regulation, the DG showed upregulation of genes associated with protein folding. We demonstrate the functional relevance of subregion-specific gene expression by genetic manipulation of a transcription factor selectively in the CA1 hippocampal subregion, leading to long-term memory deficits. This work demonstrates the power of using spatial molecular approaches to reveal transcriptional events during memory consolidation.



This repository describes our SPLIT-seq and scATAC-seq analyses between HC and SOR conditions.

The file "split-pipe code" contains the code used to run the split-pipe pipelien to pre-process the raw fastq files.

The file "SPLIT-SEQ R data analysis.R" contains the code in R language for downstream analyses, from normalizaing/scaling the data to the identification of speicifc clusters of celltype.

The file "scATAC-seq analysis.R" contains the code in R language for downstream analyses of the single-cell ATAC-seq data.

## SPLIT-seq experimental setup

![image](https://github.com/YannVRB/SPLIT-seq-SOR-data-analysis/assets/69206510/18f9b9a1-cfbf-4a10-b1a0-00124eaf8dc7)


![image](https://github.com/YannVRB/SPLIT-seq-SOR-data-analysis/assets/69206510/02b16c05-d7d1-4498-945f-1d20c05a806b)
