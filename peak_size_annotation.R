# ============================================================================
# STEP-BY-STEP GUIDE: Extract Peak Scores from BigWig Files in RStudio
# ============================================================================

# STEP 1: Install Required Packages
# Run this once - it will take a few minutes
# ============================================================================
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("rtracklayer")
BiocManager::install("GenomicRanges")


# STEP 2: Load Required Libraries
# Run this every time you start a new R session
# ============================================================================
library(rtracklayer)
library(GenomicRanges)


# STEP 3: Set Your Working Directory
# Change this to where your files are located
# ============================================================================
setwd("./")  # CHANGE THIS PATH!

# To check your current directory:
getwd()

# To see what files are in this directory:
list.files()


# STEP 4: Define Your Input Files
# Update these paths to match your actual file names
# ============================================================================
bed_file <- "bless80_hg19.bed"           # Your BED file with peak regions
bigwig_file <- "ChIP_neg_average.bw"        # Your BigWig file with signal data
output_file <- "peaks_with_scores.bed"  # Output file name


# STEP 5: Read the BED File
# This loads your peak regions
# ============================================================================
cat("Step 5: Reading BED file...\n")
peaks <- import.bed(bed_file)

# Check what you loaded:
print(peaks)
cat("Number of peaks:", length(peaks), "\n")


# STEP 6: Read the BigWig File
# This loads your signal data
# ============================================================================
cat("\nStep 6: Reading BigWig file...\n")
bw <- import.bw(bigwig_file)

# Check the bigWig data:
print(head(bw))
cat("Number of entries:", length(bw), "\n")


# STEP 7: Extract Scores for Each Peak
# This calculates the maximum score in each peak region
# ============================================================================
cat("\nStep 7: Extracting scores...\n")

# Find overlaps between peaks and bigWig data
overlaps <- findOverlaps(peaks, bw)

# Initialize score vector
scores <- numeric(length(peaks))

# Calculate max score for each peak
for (i in seq_along(peaks)) {
  # Get bigWig entries that overlap with this peak
  idx <- subjectHits(overlaps)[queryHits(overlaps) == i]
  
  if (length(idx) > 0) {
    # Extract scores and take maximum
    peak_scores <- score(bw[idx])
    scores[i] <- max(peak_scores, na.rm = TRUE)
  } else {
    scores[i] <- 0
  }
  
  # Progress indicator (every 1000 peaks)
  if (i %% 1000 == 0) {
    cat("Processed", i, "peaks...\n")
  }
}

cat("Score extraction complete!\n")


# STEP 8: Create Results Data Frame
# Organize the results in a clean format
# ============================================================================
cat("\nStep 8: Creating results table...\n")

results <- data.frame(
  chr = as.character(seqnames(peaks)),
  start = start(peaks),
  end = end(peaks),
  name = if ("name" %in% names(mcols(peaks))) {
    mcols(peaks)$name
  } else {
    paste0("peak_", seq_along(peaks))
  },
  score = scores,
  strand = if ("strand" %in% names(mcols(peaks))) {
    as.character(strand(peaks))
  } else {
    rep(".", length(peaks))
  }
)

# View first few rows
head(results)


# STEP 9: Examine the Results
# Get summary statistics
# ============================================================================
cat("\nStep 9: Summary statistics:\n")
cat("Total peaks:", nrow(results), "\n")
cat("Score range:", min(scores), "to", max(scores), "\n")
cat("Mean score:", mean(scores), "\n")
cat("Median score:", median(scores), "\n")

# Summary statistics
summary(results$score)


# STEP 10: Visualize Score Distribution
# Create a histogram of scores
# ============================================================================
cat("\nStep 10: Creating visualization...\n")

hist(results$score, 
     main = "Distribution of Peak Scores",
     xlab = "Score",
     ylab = "Number of Peaks",
     col = "lightblue",
     breaks = 50,
     border = "white")


# STEP 11: Save Results to File
# Write the output BED file
# ============================================================================
cat("\nStep 11: Saving results...\n")

write.table(results, 
            output_file, 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE, 
            col.names = FALSE)

cat("Results saved to:", output_file, "\n")


# STEP 12: View Top Peaks
# Find the peaks with highest scores
# ============================================================================
cat("\nStep 12: Top 10 peaks by score:\n")

top_peaks <- results[order(results$score, decreasing = TRUE), ][1:10, ]
print(top_peaks)


# ============================================================================
# STEP 13: ANNOTATE PEAKS WITH GENE NAMES
# ============================================================================

cat("\n\nStep 13: Annotating peaks with gene names...\n")

# Install annotation packages (run once)
# ----------------------------------------------------------------------------
if (!require("ChIPseeker", quietly = TRUE)) {
  BiocManager::install("ChIPseeker")
}
if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
}
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

# Load annotation libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# NOTE: Change these if using a different organism or genome version:
# For hg19: TxDb.Hsapiens.UCSC.hg19.knownGene
# For mm10 (mouse): TxDb.Mmusculus.UCSC.mm10.knownGene and org.Mm.eg.db
# For mm39 (mouse): TxDb.Mmusculus.UCSC.mm39.knownGene and org.Mm.eg.db

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


# STEP 14: Perform Peak Annotation
# ----------------------------------------------------------------------------
cat("\nStep 14: Performing peak annotation...\n")

# Convert results back to GRanges for annotation
# Set all strands to "*" (unstranded) if not properly formatted
if (!"strand" %in% colnames(results) || all(results$strand == ".")) {
  results$strand <- "*"
}

# Convert results back to GRanges for annotation
peaks_gr <- GRanges(
  seqnames = results$chr,
  ranges = IRanges(start = results$start, end = results$end),
  strand = results$strand
)

# Add the score to the GRanges object
mcols(peaks_gr)$score <- results$score

# Annotate peaks
peakAnno <- annotatePeak(peaks_gr, 
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")

# Convert to data frame
anno_df <- as.data.frame(peakAnno)


# STEP 15: Extract Gene Names and Add to Results
# ----------------------------------------------------------------------------
cat("\nStep 15: Adding gene names to results...\n")

# Add gene names to results
results$gene_name <- anno_df$SYMBOL
results$gene_id <- anno_df$geneId
results$annotation <- anno_df$annotation
results$distance_to_TSS <- anno_df$distanceToTSS
results$gene_strand <- anno_df$geneStrand

# View results with gene names
cat("\nResults with gene annotations:\n")
head(results[, c("chr", "start", "end", "score", "gene_name", "gene_strand", "annotation", "distance_to_TSS")])


# STEP 16: Summary of Annotations
# ----------------------------------------------------------------------------
cat("\nStep 16: Annotation summary:\n")

# Count peaks by genomic feature
cat("\nPeaks by genomic region:\n")
print(table(anno_df$annotation))

# Count unique genes
cat("\nNumber of unique genes:", length(unique(results$gene_name[!is.na(results$gene_name)])), "\n")


# STEP 17: Visualize Peak Annotations
# ----------------------------------------------------------------------------
cat("\nStep 17: Creating annotation visualizations...\n")

# Plot annotation distribution
plotAnnoPie(peakAnno)

# Plot distance to TSS
plotDistToTSS(peakAnno,
              title = "Distribution of Peaks Relative to TSS")


# STEP 18: Save Annotated Results
# ----------------------------------------------------------------------------
cat("\nStep 18: Saving annotated results...\n")

# Save full annotated results
write.table(results,
            "bless_peaks_with_scores_and_genes_ChIP_neg.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# Save simplified version (BED format with gene name)
simple_output <- data.frame(
  chr = results$chr,
  start = results$start,
  end = results$end,
  name = paste0(results$gene_name, "_", results$name),
  score = results$score,
  strand = results$strand
)

write.table(simple_output,
            "peaks_with_genes.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("Annotated results saved!\n")


# STEP 19: Find Top Scoring Peaks with Gene Names
# ----------------------------------------------------------------------------
cat("\nStep 19: Top 20 peaks with gene annotations:\n")

top_peaks_genes <- results[order(results$score, decreasing = TRUE), 
                           c("chr", "start", "end", "score", "gene_name", 
                             "annotation", "distance_to_TSS")][1:20, ]
print(top_peaks_genes)


# ============================================================================
# OPTIONAL: Additional Analyses
# ============================================================================

# Filter peaks by gene name
# ----------------------------------------------------------------------------
cat("\n\nOPTIONAL: Filtering peaks by gene name...\n")

# Example: Find peaks near specific genes
genes_of_interest <- c("TP53", "MYC", "BRCA1")  # CHANGE THESE!
peaks_near_genes <- results[results$gene_name %in% genes_of_interest, ]

if (nrow(peaks_near_genes) > 0) {
  cat("Peaks near genes of interest:\n")
  print(peaks_near_genes[, c("chr", "start", "end", "score", "gene_name", "annotation")])
} else {
  cat("No peaks found near specified genes.\n")
}


# Filter promoter peaks only
# ----------------------------------------------------------------------------
cat("\nFiltering promoter peaks...\n")
promoter_peaks <- results[grepl("Promoter", results$annotation), ]
cat("Number of promoter peaks:", nrow(promoter_peaks), "\n")
cat("Mean score in promoters:", mean(promoter_peaks$score), "\n")


# Compare scores across genomic regions
# ----------------------------------------------------------------------------
cat("\nScore distribution by genomic region:\n")

# Simplify annotation categories
results$simple_annotation <- sapply(results$annotation, function(x) {
  if (grepl("Promoter", x)) return("Promoter")
  if (grepl("Exon|5' UTR|3' UTR", x)) return("Exon")
  if (grepl("Intron", x)) return("Intron")
  if (grepl("Downstream", x)) return("Downstream")
  return("Intergenic")
})

# Boxplot by region
boxplot(score ~ simple_annotation, data = results,
        main = "Peak Scores by Genomic Region",
        xlab = "Genomic Region",
        ylab = "Score",
        col = rainbow(5),
        las = 2)


# Create summary table by gene
# ----------------------------------------------------------------------------
cat("\nCreating gene-level summary...\n")

gene_summary <- aggregate(score ~ gene_name, data = results, 
                          FUN = function(x) c(count = length(x),
                                              max = max(x),
                                              mean = mean(x)))
gene_summary <- do.call(data.frame, gene_summary)
colnames(gene_summary) <- c("gene", "peak_count", "max_score", "mean_score")

# Sort by peak count
gene_summary <- gene_summary[order(gene_summary$peak_count, decreasing = TRUE), ]

cat("\nTop 10 genes by number of peaks:\n")
print(head(gene_summary, 10))


# Save gene summary
write.table(gene_summary,
            "gene_level_summary.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nOutput files created:\n")
cat("  1. peaks_with_scores.bed - Original output with scores\n")
cat("  2. peaks_with_scores_and_genes.txt - Full annotation table\n")
cat("  3. peaks_with_genes.bed - BED file with gene names\n")
cat("  4. gene_level_summary.txt - Summary statistics by gene\n")

