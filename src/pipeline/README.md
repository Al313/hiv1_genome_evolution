# Sequencing Data Processing Pipeline

All steps for processing raw sequencing reads were implemented using a **Snakemake** pipeline.  
The current version of the pipeline consists of **10 modular components**, located in the `./module` directory.

## Pipeline Overview

The directed acyclic graph (DAG) of the workflow is shown below:

![Pipeline DAG](rulegraph.png)

---

## Workflow Summary

### 1. Quality Control
Raw sequencing reads were assessed using **FASTQC** and **MULTIQC**.  
Samples that did not meet quality standards were re-sequenced.

---

### 2. Read Trimming
Reads were processed using **Trimmomatic** with the following criteria:
- Removal of low-quality bases (**Phred score ≤ 20**, error probability ≥ 0.01)
- Removal of short reads (**length ≤ 36 bp**)

---

### 3. Read Mapping
Trimmed reads were aligned to the ancestral **HIV-1 NL4-3 consensus sequence** using **BWA-MEM** (default parameters).

- Reads with **mapping quality ≤ 35** (error probability ≥ 0.00032) were discarded

---

### 4. Coverage Assessment
- Average sequencing depth was calculated for each of the **5 amplicons**
- Amplicons with coverage below **1,000×** were re-sequenced
- The full viral genome achieved a minimum coverage of **200 high-quality reads**, enabling full-length reconstruction

---

### 5. Variant Calling
Variants were identified using **custom Python scripts** (functionally comparable to tools such as `smalt-align` and `mpileup`).

- Reference genome: **HIV-1 NL4-3 plasmid** (NCBI accession: *AF324493.2*)
- Filtering thresholds:
  - **Minor Allele Frequency (MAF) ≥ 0.01**
  - **Sequencing depth ≥ 200 reads**

---

### 6. Variant Classification

Variants were categorized based on allele frequency:

| Category        | Frequency Range        |
|----------------|----------------------|
| Fixed          | ≥ 0.99               |
| Majority       | ≥ 0.50               |
| Minority       | 0.01 ≤ freq < 0.50   |

---

### 7. Functional Annotation
Variants were annotated using **SnpEff (v5.1)** following official guidelines.

Annotations include:
- Predicted gene impact
- Amino acid substitutions
- Functional effects on proteins

The genomic feature set for the **HIV-1 NL4-3 strain** was obtained from NCBI (*AF324493.2*) and manually refined based on the reconstructed ancestral sequence.

- Includes:
  - **9 protein-coding genes**
  - **5 non-coding regions**

---

### 8. Reversion Analysis
Reversion sites were identified through gene-by-gene alignment of:
- The virus stock consensus sequence
- A 2004 cross-subtype consensus from the Los Alamos HIV Sequence Database

A mutation was classified as a **reversion** if it returned to the database consensus nucleotide in later passages.

---

### 9. Phylogenetic Reconstruction
The pipeline reconstructs phylogenetic relationships among samples using a **maximum likelihood approach**, based on the consensus sequence of each population.

---

## Summary

This pipeline provides a fully automated workflow for:
- Quality control
- Read processing
- Genome reconstruction
- Variant detection and annotation
- Evolutionary analysis

All steps are reproducible and modular, facilitating extension and adaptation.
