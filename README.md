# Deciphering the low abundance microbiota of presumed aseptic hip and knee implants

This repository contains scripts to reproduce the analyses in "Deciphering the low abundance microbiota of presumed aseptic hip and knee implants".

## Contents

`01_custom_dada2.R`: DADA2 pipeline for reads resulting from the custom demultiplexing script

`02_load_requirements.R`: loads the packages, external functions, and data required for the analysis

`03_filtering.R`: filters out low-abundance and rare ASVs + low-depth samples

`04_heatmap_and_aitchison_dists.R`: generates Figure 2

`05_sourcetracker.R`: uses SourceTracker to select patient samples

`06_decontam.R`: executes the decontam pipeline

`07_data_reorg.R`: removes contaminants and re-organizes data before further analysis

`08_envfit_and_pca.R`: generates PCA biplots and performs envfit testing required for Figures 3 and 4

`09_shannon_diversity.R`: computes Shannon diversity and performs Shannon diversity comparisons required for Figures 3 and 4

`10_aldex2.R`: performs and plots ALDEx2 testing required for Figures 3 and 4

`11_assemble_figs_3_4.R`: assembles Figures 3 and 4

`12_barplot.R`: generates Figure 5 

`13_select_for_pcr.R`: pseudorandomly selects samples to examine with species-specific PCR

`14_run_cutadapt.sh`: performs demultiplexing and primer/adapter removal with cutadapt

`15_cutadapt.R`: analyzes cutadapt-processed reads (DADA2 to Supplementary Figure 3)

## Data Availability

To use these scripts, please download demultiplexed reads from the Sequence Read Archive (PRJNA726194) and/or data from Zenodo (10.5281/zenodo.5136551).

## Citation

Carr C, Wilcox H, Burton JP, Menon S, Al KF, O’Gorman D, et al. (2021) Deciphering the low abundance microbiota of presumed aseptic hip and knee implants. PLoS ONE 16(9): e0257471. https://doi.org/10.1371/journal.pone.0257471
