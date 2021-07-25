# This script loads the packages, external functions, and data required for the analysis

# Load packages
library(ggplot2) # Plotting
library(ggpubr) 

library(tibble) # Data wrangling
library(tidyr)
library(forcats)

library(zCompositions) # PCA

library(vegan) # Diversity calculations

library(ALDEx2) # ALDEx2

library(phyloseq) # decontam
library(decontam) 

# Load external functions
source('/home/ccarr/Documents/lab/sourcetracker-master/src/SourceTracker.r') # SourceTracker

# Load data
counts <- read.table("/home/ccarr/Documents/lab/implant_microbiota_study/data/counts.txt", 
                     header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
counts <- t(counts)
counts <- counts[, order(colnames(counts))]

tax <- read.table("/home/ccarr/Documents/lab/implant_microbiota_study/data/tax.txt", 
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

meta <- read.table("/home/ccarr/Documents/lab/implant_microbiota_study/data/meta.txt", 
                   header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
meta <- meta[order(rownames(meta)), ]
meta[meta == " PowerSoil"] <- "PowerSoil"

st_meta <- read.table("/home/ccarr/Documents/lab/implant_microbiota_study/data/meta_st.txt", 
                      header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
st_meta <- st_meta[order(rownames(st_meta)), ]

conc_data <- read.table("/home/ccarr/Documents/lab/implant_microbiota_study/data/dna_concs.txt", 
                        header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
