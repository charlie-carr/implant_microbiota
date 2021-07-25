# DADA2 pipeline for reads resulting from the custom demultiplexing script.
# The code is very minimally modified from the tutorial at 
# https://benjjneb.github.io/dada2/tutorial.html

# Load packages
library(dada2)
library(dplyr)

# Establish the file path to the directory containing the demultiplexed samples
path <- "/home/ccarr/Documents/lab/implant_microbiota_study/demultiplex_reads"
list.files(path) # Check that all the files are there 

# Forward and reverse fastq filenames have the format:  
# SAMPLENAME-R1.fastq and SAMPLENAME-R2.fastq
# Use this fact to sort them into two groups
fnFs <- sort(list.files(path, pattern = "-R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "-R2.fastq", full.names = TRUE)) 

# Extract sample names (i.e., exclude the forward/reverse identifier), 
# assuming filenames have format SAMPLENAME-Rn.fastq and SAMPLENAME
# does not include any "-"
sample_names <- sapply(strsplit(basename(fnFs), "-R1"), `[`, 1)
any(duplicated(sample_names)) # FALSE, so we can proceed

# Grab four samples, the reads of which will be examined in terms of 
# their quality profiles
ids <- round(runif(4,1,length(sample_names)))

# Output the quality profiles for the forward and reverse reads of 
# those samples
pdf("/home/ccarr/Documents/lab/implant_microbiota_study/figures/final_qualprofiles.pdf")
plotQualityProfile(fnFs[ids])
plotQualityProfile(fnRs[ids])
dev.off()

# Create a directory and names for the filtered files
filtFs <- file.path("/home/ccarr/Documents/lab/implant_microbiota_study/filtered_reads", 
                    paste0(sample_names, "-F-filt.fastq"))
filtRs <- file.path("/home/ccarr/Documents/lab/implant_microbiota_study/filtered_reads", 
                    paste0(sample_names, "-R-filt.fastq"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Perform quality filtering 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(175, 155), maxEE = c(2, 2), 
                     compress = TRUE, multithread = TRUE, verbose = TRUE)

# Learn the error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot the error rates
pdf("/home/ccarr/Documents/lab/implant_microbiota_study/figures/error_plot_final_F.pdf")
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf("/home/ccarr/Documents/lab/implant_microbiota_study/figures/error_plot_final_R.pdf")
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Use the filtered files and error rates to perform
# sample inference (without pooling)
dadaFs <- dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)

# Merge the forward and reverse paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Build the counts table
counts_raw <- makeSequenceTable(mergers)

# Check that all of the sequence lengths are within the expected range
table(nchar(getSequences(counts_raw)))
# Output:
# 175 177 178 195 209 210 211 215 239 240 241 242 249 266 268 291 
# 29   3   1   4   4   1   2   1   1  46 607  19   1   1   1   1 

# Most ASVs are acceptable lengths, but the others must be removed
counts_trimmed <- counts_raw[, nchar(colnames(counts_raw)) %in% seq(240, 246)]

# Filter chimeras
counts_nochim <- removeBimeraDenovo(counts_trimmed, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)
# Output:
# Identified 78 bimeras out of 672 input sequences.

sum(counts_nochim)/sum(counts_trimmed)
# Output:
# 0.9964818

# Assign taxonomy
tax_nochim <- assignTaxonomy(counts_nochim,
                             "/home/ccarr/Documents/lab/tax_data/silva_nr99_v138_train_set.fa.gz", 
                             multithread = TRUE)

tax_nochim <- addSpecies(tax_nochim, 
                         "/home/ccarr/Documents/lab/tax_data/silva_species_assignment_v138.fa.gz")

# Filter by taxonomy
tax_filtered <- as.data.frame(tax_nochim) %>%
  filter(!is.na(Kingdom)) %>%
  filter(Kingdom != "Eukaryota") %>%
  filter(Family != "Mitochondria") %>%
  filter(Order != "Chloroplast")
  
counts_filtered <- counts_nochim[, rownames(tax_filtered)]

# Assign readable names
tax_filtered <- tax_filtered %>%
  mutate(Sequence = rownames(tax_filtered))
rownames(tax_filtered) <- paste0("SV_", 1:nrow(tax_filtered))

any(colnames(counts_filtered) != tax_filtered$Sequence) # FALSE, so the ASVs are in the same order

colnames(counts_filtered) <- paste0("SV_", 1:ncol(counts_filtered))

# Construct a table to summarize the removal of reads throughout the pipeline
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(counts_trimmed), rowSums(counts_nochim), rowSums(counts_filtered))
colnames(track) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Trimmed", "Non-Chimeric", "Tax-Filtered")
rownames(track) <- sample_names

# Output the counts, tax, and tracking tables
write.table(counts_filtered, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/counts.txt", 
            sep = "\t", col.names = NA, quote = F)

write.table(tax_filtered, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/tax.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(track, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/track.txt",
            sep = "\t", col.names = NA, quote = F)
