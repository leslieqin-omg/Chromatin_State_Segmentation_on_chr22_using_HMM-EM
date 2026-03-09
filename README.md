# Chromatin-State-Segmentation-on-chr22-using-HMM-EM-Baum-Welch-
# HMM-based Chromatin State Segmentation

## Overview
This repository contains the code and methodology for unsupervised chromatin state segmentation on human chromosome 22 (hg38). Using continuous ChIP-seq bigWig tracks, the pipeline employs a multivariate Gaussian Hidden Markov Model (HMM) to convert noisy, bin-level epigenetic signals into interpretable, discrete regulatory annotations. 

The core model is trained on three highly informative histone marks (H3K4me3, H3K27ac, and H3K36me3) to discover recurrent chromatin states, with built-in extensions for robust model selection and biological validation against known genomic landmarks.


## Repository Structure
* `data/` - Place your raw `.BigWig` files here (e.g., `H3K4me3.BigWig`).
* `scripts/` - R Markdown and script files containing the analysis pipeline.
* `results/` - Output tables (CSV), decoded state tracks (BED), and generated PDF/SVG plots.
* `docs/` - Project manuscript and supplementary diagnostic reports.

## Prerequisites
All analyses are performed in **R**. The following packages are required:

**Core Modeling & Genomics:**
* `depmixS4` (Hidden Markov Models)
* `GenomicRanges`, `rtracklayer`, `GenomeInfoDb` (Genomic intervals and bigWig I/O)
* `TxDb.Hsapiens.UCSC.hg38.knownGene`, `GenomicFeatures` (Transcriptome annotations)

**Visualization & Data Manipulation:**
* `ggplot2`, `pheatmap`, `corrplot`, `gridExtra`
* `reshape2`

## Pipeline Workflow

### 1. Data Preprocessing
* **Binning:** Chromosome 22 is partitioned into non-overlapping 200-bp bins using `tileGenome`.
* **Signal Extraction:** Average coverage is calculated for each bin.
* **Normalization:** To stabilize variance and prevent extreme peaks from dominating the Gaussian estimations, the raw signal undergoes a $log(1+x)$ transformation followed by per-mark z-score standardization.

### 2. HMM Training & Decoding
* **Model:** A multivariate Gaussian HMM with diagonal covariance.
* **Training:** Parameters are estimated via the Expectation-Maximization (EM) algorithm (Baum-Welch).
* **Decoding:** The Viterbi algorithm is used to extract the most likely contiguous chromatin state path.
* *Note:* The pipeline includes a multi-restart wrapper to avoid local optima during EM training.

### 3. Biological Validation
To interpret the functional relevance of the unsupervised states, the pipeline computes spatial fold-enrichment against orthogonal datasets:
* **Transcription Start Sites (TSS)**
* **CpG Islands** (Directly fetched from UCSC GoldenPath)
* **Gene Bodies** (Excluding $TSS \pm 2kb$ to isolate transcription-elongation signatures)

### 4. Model Diagnostics & Evaluation
The codebase includes comprehensive diagnostic suites for model selection and stability:
* **Information Criteria:** AIC / BIC evaluation across $K=2$ to $K=8$ states.
* **Persistence Analysis:** Empirical transition matrices and state segment-length distributions.
* **Residual Analysis:** Evaluation of conditional normality assumptions.
* **Correlation:** Pearson correlation matrices of input epigenetic marks.

## Identified Chromatin States (K=4 Model)
Based on the emission parameters and biological enrichment, the standard 4-state model resolves into:
1. **Promoter-like State:** High H3K4me3/H3K27ac; strongly enriched at TSS and CpG islands.
2. **Gene-body-like State:** High H3K36me3; preferentially enriched within transcribed gene bodies.
3. **High-activity Composite State:** Elevated signal across multiple marks, representing broadly active or dense regulatory hubs.
4. **Background-like State:** Uniformly low signal representing quiescent or unannotated chromatin.

## Usage
1. Clone this repository and ensure the required `data/*.BigWig` files are in place.
2. Open `final_project.Rmd` in RStudio.
3. Execute the chunks sequentially. The script will automatically generate standardized plots (e.g., `Figure_EmissionHeatmap.pdf`, `TSS_enrich_validation.pdf`) and export the final segmentation as `K562_HMM_Segmentation.bed` for direct UCSC Genome Browser visualization.

## Author
**Haifeng Qin**
