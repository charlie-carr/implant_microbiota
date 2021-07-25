# This script filters out low-abundance and rare ASVs + low-depth samples

# Find initial samples, ASVs, and reads
sum(counts) # Output: 8177228
dim(counts) # Output: 570 307

# Generate a relative abundance table to remove ASVs accounting for less than 1% of reads in every sample
props <- apply(counts, 2, function(x) {x/sum(x)})
filtered_by_props <- counts[apply(props, 1, max) >= 0.01, ]

# Re-check dimensions
sum(filtered_by_props) # Output: 8060299
dim(filtered_by_props) # Output: 346 307

any(colSums(filtered_by_props) == 0) # Output: FALSE

# Filter ASVs and samples according to some additional filtering parameters. To pick these parameters, generate PCA biplots
# for the data sets that result from different filtering conditions
test_filt_params <- function(counts_table, sv_param, sample_param) {

  # Perform filtering
  filtered_counts_0 <- counts_table[rowSums(counts_table) >= sv_param, ]
  filtered_counts <- filtered_counts_0[, colSums(filtered_counts_0) >= sample_param]

  # CZM, CLR, PCA
  czm <- cmultRepl(t(filtered_counts), label = 0, method = "CZM")
  clr <- t(apply(czm, 1, function(x){log(x) - mean(log(x))}))
  pca <- prcomp(clr)

  # Select coordinates of samples and features
  sample_positions <- data.frame(pca[["x"]])
  sv_positions <- data.frame(pca[["rotation"]])

  # Produce the PCA biplot
  ggplot(sample_positions, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.7) + # Plot samples
    geom_segment(data = sv_positions, aes(x = 0, y = 0, xend = 15 * PC1, yend = 15 * PC2),
                 arrow = arrow(length = unit(1/2, 'picas')),
                 color = "#DC0000") + # Plot features
    geom_text(data = sv_positions, aes(x = 15 * PC1, y = 15 * PC2), nudge_x = 0.5, nudge_y = 0.25, check_overlap = TRUE,
              label = rownames(sv_positions),
              color = "#DC0000") + # Label features. To show them all, check_overlap = FALSE
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0 , linetype = "dashed") +
    labs(title = paste0("ASVs \u2265 ", sv_param, "; Samples \u2265 ", sample_param),
         x = paste("PC1: ", round(pca$sdev[1]^2/sum(pca$sdev^2),3)),
         y = paste("PC2: ", round(pca$sdev[2]^2/sum(pca$sdev^2),3)),
         caption = paste0("ASVs: ", nrow(filtered_counts), "; Samples: ", ncol(filtered_counts), "; Reads: ", sum(filtered_counts))) +
    theme_classic()

}

# Store all PCA biplots for inspection
for (i in seq(50, 650, 200)) {
  for (j in seq(50, 650, 200)) {
    assign(paste0("plot_", i, "_", j), test_filt_params(filtered_by_props, i, j))
  }
}

# Build the overall figure
varied_filtering <- ggarrange(plot_50_50, plot_50_250, plot_50_450, plot_50_650,
                              plot_250_50, plot_250_250, plot_250_450, plot_250_650,
                              plot_450_50, plot_450_250, plot_450_450, plot_450_650,
                              plot_650_50, plot_650_250, plot_650_450, plot_650_650,
                              labels = c("A", "B", "C", "D",
                                         "E", "F", "G", "H",
                                         "I", "J", "K", "L",
                                         "M", "N", "O", "P"),
                              ncol = 4, nrow = 4)

tiff("/home/ccarr/Documents/lab/implant_microbiota_study/figures/varied_filtering.tiff", units = "in", width = 15, height = 10, res = 150)
varied_filtering
dev.off()

# Tidy environment
rm(i, j)

rm(plot_50_50,  plot_50_250,  plot_50_450,  plot_50_650,
   plot_250_50, plot_250_250, plot_250_450, plot_250_650,
   plot_450_50, plot_450_250, plot_450_450, plot_450_650,
   plot_650_50, plot_650_250, plot_650_450, plot_650_650)

# Perform filtering
filtered_counts_0 <- filtered_by_props[rowSums(filtered_by_props) >= 250, ] # ASVs
filtered_counts <- filtered_counts_0[, colSums(filtered_counts_0) >= 50] # Samples

filtered_tax_0 <- tax[rownames(filtered_counts_0), ]
filtered_tax <- tax[rownames(filtered_counts), ]

filtered_meta_0 <- meta[colnames(filtered_counts_0), ]
filtered_meta <- meta[colnames(filtered_counts), ]
