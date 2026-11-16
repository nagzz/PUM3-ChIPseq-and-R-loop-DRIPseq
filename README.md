# PUM3 ChIP-seq and R-loop DRIP-seq analysis

This repository currently contains an R workflow for summarizing ChIP-seq peak signal from a BigWig track, annotating peaks to nearby genes, and generating peak- and gene-level summaries. The repository name also refers to R-loop DRIP-seq, but the available script focuses on BED peak regions and BigWig signal; additional DRIP-seq analysis steps are not currently documented here.

> **Research-use note:** This is a study-oriented analysis script, not a validated clinical or turnkey pipeline. It expects local input files and a particular genome annotation. Review the code, confirm coordinate conventions and reference-genome compatibility, and adapt file paths before use. The underlying sequencing data are not included.

## Repository contents

| File | Purpose |
| --- | --- |
| [`peak_size_annotation.R`](peak_size_annotation.R) | Imports BED peak regions and a BigWig signal track, assigns each peak its maximum overlapping signal, summarizes and plots scores, and annotates peaks using ChIPseeker with the hg19 known-gene annotation. It also creates gene-level summaries and example filtered subsets. |

## What the script does

`peak_size_annotation.R` is organized as a step-by-step RStudio workflow:

1. Loads genomic ranges and BigWig import tools.
2. Imports a BED file of peak regions and a BigWig signal file.
3. Finds overlaps and assigns each peak the maximum BigWig score over its overlapping signal intervals (or zero when there is no overlap).
4. Summarizes peak scores and plots their distribution.
5. Annotates peaks relative to genes and transcription start sites using ChIPseeker and the UCSC hg19 known-gene TxDb.
6. Creates annotation plots, example promoter/gene-of-interest subsets, a score-by-feature plot, and a gene-level summary.

## Requirements

The script uses R and Bioconductor packages:

- `rtracklayer`
- `GenomicRanges`
- `ChIPseeker`
- `TxDb.Hsapiens.UCSC.hg19.knownGene`
- `org.Hs.eg.db`

Install packages as needed from an R session; for example:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "rtracklayer",
  "GenomicRanges",
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "org.Hs.eg.db"
))
```

## Inputs and configuration

The script's example input filenames are set near the top of the file:

- `bless80_hg19.bed` — BED peak regions
- `ChIP_neg_average.bw` — BigWig signal track

Edit these variables and the working directory to point to your own files. The annotation section uses hg19; use a matching genome build for both the peak coordinates and the BigWig track. If your data use a different organism or assembly, replace the TxDb and organism annotation package accordingly.

To run the script from a terminal after updating its input paths:

```sh
Rscript peak_size_annotation.R
```

The file contains interactive plot calls, so it can also be opened and run in RStudio. The script writes outputs to the current working directory.

## Outputs

The script writes or displays the following results:

- `peaks_with_scores.bed` — peak coordinates, names, maximum signal scores, and strand.
- `bless_peaks_with_scores_and_genes_ChIP_neg.txt` — peak table with gene and genomic-feature annotations.
- `peaks_with_genes.bed` — BED-style peak table with gene names in the feature label.
- `gene_level_summary.txt` — per-gene peak counts and score summaries.
- Histograms and annotation plots displayed by the active R graphics device.

Check output filenames and the current working directory before running. Some comments near the end of the script mention a different full-annotation filename than the one actually written; the code writes `bless_peaks_with_scores_and_genes_ChIP_neg.txt`.

## Notes for adapting the analysis

- The script calculates the maximum signal among BigWig intervals overlapping a peak; this is not the same as the mean signal across the full peak region.
- Peak annotation uses a ±3 kb TSS window and hg19 gene models. These choices should match the biological question and reference assembly.
- Example genes of interest (`TP53`, `MYC`, and `BRCA1`) are placeholders that can be changed for the study.
- Review BED coordinate conventions, chromosome naming (for example, `chr1` versus `1`), and genome-build consistency before interpreting overlaps or annotations.
- The repository does not currently declare a license. Add one if you want to specify reuse and redistribution terms.
