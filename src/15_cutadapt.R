# This script the analyzes cutadapt-processed reads (DADA2 to Supplementary Figure 3)

# Establish the file path to the directory containing the demultiplexed samples
cutadapt_path <- "/home/ccarr/Documents/lab/implant_microbiota_study/cutadapt_reads/trimmed_reads"
list.files(cutadapt_path) # Check that all the files are there

# Forward and reverse fastq filenames have the format:
# SAMPLENAME.1.fastq.gz and SAMPLENAME.2.fastq.gz
# Use this fact to sort them into two groups
cutadapt_fnFs <- sort(list.files(cutadapt_path, pattern = ".1.fastq.gz", full.names = TRUE))
cutadapt_fnRs <- sort(list.files(cutadapt_path, pattern = ".2.fastq.gz", full.names = TRUE))

# Extract sample names (i.e., exclude the forward/reverse identifier),
# assuming filenames have format SAMPLENAME-Rn.fastq and SAMPLENAME
# does not include any "-"
cutadapt_sample_names <- sapply(strsplit(basename(cutadapt_fnFs), ".1.fastq.gz"), `[`, 1)
any(duplicated(cutadapt_sample_names)) # FALSE, so we can proceed

# Grab four samples, the reads of which will be examined in terms of
# their quality profiles
cutadapt_ids <- round(runif(4, 1, length(cutadapt_sample_names)))

# Output the quality profiles for the forward and reverse reads of
# those samples
pdf("/home/ccarr/Documents/lab/implant_microbiota_study/figures/cutadapt_final_qualprofiles.pdf")
plotQualityProfile(cutadapt_fnFs[cutadapt_ids])
plotQualityProfile(cutadapt_fnRs[cutadapt_ids])
dev.off()

# Create a directory and names for the filtered files
cutadapt_filtFs <- file.path("/home/ccarr/Documents/lab/implant_microbiota_study/cutadapt_filtered_reads",
                             paste0(cutadapt_sample_names, "-F-filt.fastq"))
cutadapt_filtRs <- file.path("/home/ccarr/Documents/lab/implant_microbiota_study/cutadapt_filtered_reads",
                             paste0(cutadapt_sample_names, "-R-filt.fastq"))
names(cutadapt_filtFs) <- cutadapt_sample_names
names(cutadapt_filtRs) <- cutadapt_sample_names

# Actually perform filtering on the basis of quality
cutadapt_out <- filterAndTrim(cutadapt_fnFs, cutadapt_filtFs, cutadapt_fnRs, cutadapt_filtRs, truncLen = c(175,155), truncQ = 2,
                              maxN = 0, maxEE = c(2,2), rm.phix = TRUE, compress = TRUE,
                              verbose = TRUE, multithread = TRUE)

# Learn the error rates
cutadapt_errF <- learnErrors(cutadapt_filtFs, multithread = TRUE)
cutadapt_errR <- learnErrors(cutadapt_filtRs, multithread = TRUE)

# Plot the error rates
pdf("/home/ccarr/Documents/lab/implant_microbiota_study/figures/cutadapt_error_plot_final_F.pdf")
plotErrors(cutadapt_errF, nominalQ = TRUE)
dev.off()

pdf("/home/ccarr/Documents/lab/implant_microbiota_study/figures/cutadapt_error_plot_final_R.pdf")
plotErrors(cutadapt_errR, nominalQ = TRUE)
dev.off()

# Use the filtered files and error rates to perform
# sample inference (without pooling)
cutadapt_dadaFs <- dada(cutadapt_filtFs, err = cutadapt_errF, pool = FALSE, multithread = TRUE)
cutadapt_dadaRs <- dada(cutadapt_filtRs, err = cutadapt_errR, pool = FALSE, multithread = TRUE)

# Merge the forward and reverse paired reads
cutadapt_mergers <- mergePairs(cutadapt_dadaFs, cutadapt_filtFs, cutadapt_dadaRs, cutadapt_filtRs, verbose = TRUE)

# Build the counts table
cutadapt_counts_raw <- makeSequenceTable(cutadapt_mergers)

# Check that all of the sequence lengths are within the expected range
table(nchar(getSequences(cutadapt_counts_raw)))
# Output:
# 175  183  185  186  187  188  195  197  203  207  210  218  221  222  223  227  250  251  252  253  254  261  270  273  278  280  303
# 5    1    2   22    7    1    1    1    3    5    2    1    5    1    2    1    2    1  107 1670   31    1    3    4    1    1    1

# Most reads are acceptable lengths, but the others must be removed
cutadapt_counts_trimmed <- cutadapt_counts_raw[, nchar(colnames(cutadapt_counts_raw)) %in% seq(250, 256)]

# Filter chimeras
cutadapt_counts_nochim <- removeBimeraDenovo(cutadapt_counts_trimmed, method = "consensus",
                                             multithread = TRUE, verbose = TRUE)
# Output: Identified 698 bimeras out of 1811 input sequences.

sum(cutadapt_counts_nochim)/sum(cutadapt_counts_trimmed)
# Output: [1] 0.9961283

# Assign taxonomy
cutadapt_tax_nochim <- assignTaxonomy(cutadapt_counts_nochim,
                                      "/home/ccarr/Documents/lab/tax_data/silva_nr99_v138_train_set.fa.gz",
                                      multithread = TRUE)

cutadapt_tax_nochim <- addSpecies(cutadapt_tax_nochim,
                                  "/home/ccarr/Documents/lab/tax_data/silva_species_assignment_v138.fa.gz")

# Filter by taxonomy
cutadapt_tax_filtered <- as.data.frame(cutadapt_tax_nochim) %>%
  filter(!is.na(Kingdom)) %>%
  filter(Kingdom != "Eukaryota") %>%
  filter(Family != "Mitochondria") %>%
  filter(Order != "Chloroplast")

cutadapt_counts_filtered <- cutadapt_counts_nochim[, rownames(cutadapt_tax_filtered)]

# Assign readable names
cutadapt_tax_filtered <- cutadapt_tax_filtered %>%
  mutate(Sequence = rownames(cutadapt_tax_filtered))
rownames(cutadapt_tax_filtered) <- paste0("SV_", 1:nrow(cutadapt_tax_filtered))

any(colnames(cutadapt_counts_filtered) != cutadapt_tax_filtered$Sequence) # FALSE, so the ASVs are in the same order

colnames(cutadapt_counts_filtered) <- paste0("SV_", 1:ncol(cutadapt_counts_filtered))

# Construct a table to summarize the removal of reads throughout the pipeline
getN <- function(x) {sum(getUniques(x))}
cutadapt_track <- cbind(cutadapt_out, sapply(cutadapt_dadaFs, getN), sapply(cutadapt_dadaRs, getN), sapply(cutadapt_mergers, getN), rowSums(cutadapt_counts_trimmed), rowSums(cutadapt_counts_nochim), rowSums(cutadapt_counts_filtered))
colnames(cutadapt_track) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Trimmed", "Non-Chimeric", "Tax-Filtered")
rownames(cutadapt_track) <- cutadapt_sample_names

# Output the counts, tax, and tracking tables
write.table(cutadapt_counts_filtered,
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/cutadapt_counts.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(cutadapt_tax_filtered,
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/cutadapt_tax.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(cutadapt_track,
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/cutadapt_track.txt",
            sep = "\t", col.names = NA, quote = F)

# Load metadata
cutadapt_meta <- read.table("/home/ccarr/Documents/lab/implant_microbiota_study/data/cutadapt_meta.txt",
                           header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)


# Filter
cutadapt_counts <- cutadapt_counts_filtered[rownames(cutadapt_meta), ]

cutadapt_counts <- t(cutadapt_counts)

cutadapt_props <- apply(cutadapt_counts, 2, function(x) {x/sum(x)}) # Generate a relative abundance table
cutadapt_filt_by_props <- cutadapt_counts[apply(cutadapt_props, 1, max) >= 0.01, ] # Remove ASVs accounting for less than 1% of reads in every sample

cutadapt_filtered_counts_0 <- cutadapt_filt_by_props[rowSums(cutadapt_filt_by_props) >= 250, ] # ASVs

cutadapt_filtered_tax_0 <- cutadapt_tax_filtered[rownames(cutadapt_filtered_counts_0), ]

cutadapt_filtered_meta_0 <- cutadapt_meta[colnames(cutadapt_filtered_counts_0), ]

# Convert to relative abundance data
cutadapt_counts_by_genus_0 <- aggregate(cutadapt_filtered_counts_0, by = list(cutadapt_filtered_tax_0$Genus), FUN = sum)
rownames(cutadapt_counts_by_genus_0) <- cutadapt_counts_by_genus_0$Group.1
cutadapt_counts_by_genus_0$Group.1 <- NULL

cutadapt_props_by_genus_0 <- apply(cutadapt_counts_by_genus_0, 2, function(x) {x/sum(x)})

# Plot
cutadapt_esch_shig_heatmap <- location_heatmap(cutadapt_props_by_genus_0, "Escherichia-Shigella", cutadapt_filtered_meta_0$Plate, cutadapt_filtered_meta_0$Row, cutadapt_filtered_meta_0$Column)

cutadapt_esch_shig_heatmap <- cutadapt_esch_shig_heatmap +
  theme(legend.position = "bottom")

cutadapt_staph_heatmap <- location_heatmap(cutadapt_props_by_genus_0, "Staphylococcus", cutadapt_filtered_meta_0$Plate, cutadapt_filtered_meta_0$Row, cutadapt_filtered_meta_0$Column)

supp_cutadapt_fig <- ggarrange(cutadapt_staph_heatmap, cutadapt_esch_shig_heatmap,
                               nrow = 1,
                               labels = c("A", "B"),
                               common.legend = TRUE,
                               legend = "bottom")

tiff("/home/ccarr/Documents/lab/implant_microbiota_study/figures/supp_cutadapt_fig.tiff", units = "in", width = 7.5, height = 3.4, res = 300)
supp_cutadapt_fig
dev.off()
