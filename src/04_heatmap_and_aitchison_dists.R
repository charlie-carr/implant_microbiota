# This script generates Figure 2

# Function to build heatmaps
location_heatmap <- function(counts, taxon, plates, rows, cols) {

  # Create a long data frame of abundance, plate, row, and column information
  abun <- as.data.frame(counts[grep(taxon, rownames(counts)), ])
  abun$Sample <- colnames(counts)
  colnames(abun) <- c("Reads", "Sample")

  plate_data <- as.data.frame(plates) %>%
    as_tibble() %>%
    rowid_to_column(var = "Sample")
  plate_data$Sample <- colnames(counts)

  row_data <- as.data.frame(rows) %>%
    as_tibble() %>%
    rowid_to_column(var = "Sample")
  row_data$Sample <- colnames(counts)

  col_data <- as.data.frame(cols) %>%
    as_tibble() %>%
    rowid_to_column(var = "Sample")
  col_data$Sample <- colnames(counts)

  grouped <- merge(abun, plate_data, by = "Sample", all = FALSE)
  grouped <- merge(grouped, row_data, by = "Sample", all = FALSE)
  grouped <- merge(grouped, col_data, by = "Sample", all = FALSE)

  # Set grey colour for the wells that were not used to sequence samples
  # for this study
  na_df <- data.frame(Sample = "Filler", Reads = NA,
                      plates = c(3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
                      rows = c("7", "8", "2", "3", "6", "7", "6", "2", "3", "4", "5", "7", "8"),
                      cols = c(12, 12, 1, 1, 2, 2, 3, 4, 4, 4, 4, 4, 4))

  grouped <- rbind(grouped, na_df)

  # Reorder factors for an organized plot
  grouped$rows <- factor(grouped$rows)
  grouped$rows <- fct_relevel(grouped$rows,
                              "8", "7", "6", "5",
                              "4", "3", "2", "1")
  levels(grouped$rows) <- c("H", "G", "F", "E",
                            "D", "C", "B", "A")

  grouped$cols <- factor(grouped$cols)
  grouped$cols <- fct_relevel(grouped$cols,
                              "1", "2", "3", "4", "5", "6",
                              "7", "8", "9", "10", "11", "12")

  # Add an extra label to the cells that contained the positive controls
  if (taxon == "Escherichia-Shigella") {
    grouped <- grouped %>%
      mutate(annotation = ifelse(Sample %in% c("Spike_1_1D3", "Spike_1_2E6", "Spike_1_3C4"), "C", ""))
  } else {
    grouped <- grouped %>%
      mutate(annotation = ifelse(Sample %in% c("Spike_2_1F8", "Spike_2_2C8", "Spike_2_3D7"), "C", ""))
  }

  # Create a vector to rename facets
  facet_labs <- paste("Plate", 1:4)
  names(facet_labs) <- 1:4

  # Plot
  ggplot(grouped, aes(x = cols, y = rows, fill = Reads)) +
    geom_tile() +
    scale_fill_distiller(palette = "YlGnBu", name = "Proportional Abundance    ", na.value = "grey92") +
    theme_classic() +
    facet_wrap(~plates, scales = "free_x", labeller = labeller(plates = facet_labs)) +
    labs(title = NULL, x = "Column", y = "Row") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5, hjust = 1),
          panel.grid = element_blank(), legend.position = "none") +
    geom_text(aes(label = annotation, fontface = "bold"))

}

# Collapse by genera. To ensure that patterns are clear, use counts data
# filtered only by ASVs (not samples) without removal of the decontam-identified
# contaminants
counts_by_genus_0 <- aggregate(filtered_counts_0, by = list(filtered_tax_0$Genus), FUN = sum)
rownames(counts_by_genus_0) <- counts_by_genus_0$Group.1
counts_by_genus_0$Group.1 <- NULL

# Convert to relative abundance data
props_by_genus_0 <- apply(counts_by_genus_0, 2, function(x) {x/sum(x)})

# Store plots for both spikes
esch_shig_heatmap <- location_heatmap(props_by_genus_0, "Escherichia-Shigella", filtered_meta_0$Plate, filtered_meta_0$Row, filtered_meta_0$Column)

esch_shig_heatmap <- esch_shig_heatmap +
  theme(legend.position = "bottom")

staph_heatmap <- location_heatmap(props_by_genus_0, "Staphylococcus", filtered_meta_0$Plate, filtered_meta_0$Row, filtered_meta_0$Column)

# Compute the Aitchison distances between samples and spikes
filtered_counts_0_czm <- cmultRepl(t(filtered_counts_0), label = 0, method = "CZM")
filtered_counts_0_clr <- apply(filtered_counts_0_czm, 1, function(x){log(x) - mean(log(x))})

aitch_dist <- as.matrix(dist(t(filtered_counts_0_clr)))

# Compute proportional abundance data frame
filtered_props_0 <- apply(filtered_counts_0, 2, function(x){x/sum(x)})

# Select all samples with at least 50% of their reads accounted for by the
# positive control ASV, then find the required Aitchison distances
dist_finder <- function(spike_name, ctrl_type, plate_no) {

  # Select samples that satisfy the abundance cutoff
  if (ctrl_type == "E") {
    abund_limit <- colnames(filtered_props_0)[filtered_props_0["SV_2", ] >= 0.5]
  } else {
    abund_limit <- colnames(filtered_props_0)[filtered_props_0["SV_1", ] >= 0.5]
  }

  # Select samples from the same plate and sort them according to whether they passed the previous threshold
  contam <- Reduce(intersect, list(abund_limit, rownames(filtered_meta_0)[filtered_meta_0$Plate == plate_no]))
  clean <- setdiff(rownames(filtered_meta_0)[filtered_meta_0$Plate == plate_no], contam)

  # Find the distances between the positive control and all samples
  contam_dists <- as.data.frame(aitch_dist[grep(spike_name, rownames(aitch_dist)), contam[!contam %in% c(spike_name)]])
  clean_dists <- as.data.frame(aitch_dist[grep(spike_name, rownames(aitch_dist)), clean[!clean %in% c(spike_name)]])

  # Perform and output a Mann-Whitney U test
  test <- wilcox.test(aitch_dist[grep(spike_name, rownames(aitch_dist)), contam[!contam %in% c(spike_name)]],
                      aitch_dist[grep(spike_name, rownames(aitch_dist)), clean[!clean %in% c(spike_name)]])

  print(test)
  print(length(contam[!contam %in% c(spike_name)]))
  print(length(clean[!clean %in% c(spike_name)]))

  # Assemble all distances into a data frame for easy plotting
  contam_dists <- mutate(contam_dists, "Artifact-Affected", paste("Plate", plate_no))
  colnames(contam_dists) <- c("Distance", "Type", "Plate")
  clean_dists <- mutate(clean_dists, "Artifact-Free", paste("Plate", plate_no))
  colnames(clean_dists) <- c("Distance", "Type", "Plate")

  all_dists <- as.data.frame(rbind(contam_dists, clean_dists))

}

# Run the distance function for all spikes
plate_1_spike_1 <- dist_finder("Spike_1_1D3", "E", 1)
plate_2_spike_1 <- dist_finder("Spike_1_2E6", "E", 2)
plate_3_spike_1 <- dist_finder("Spike_1_3C4", "E", 3)

plate_1_spike_2 <- dist_finder("Spike_2_1F8", "S", 1)
plate_2_spike_2 <- dist_finder("Spike_2_2C8", "S", 2)
plate_3_spike_2 <- dist_finder("Spike_2_3D7", "S", 3)

# Combine the distances together by spike
spike_1_dists <- rbind(plate_1_spike_1, plate_2_spike_1, plate_3_spike_1)

spike_2_dists <- rbind(plate_1_spike_2, plate_2_spike_2, plate_3_spike_2)

# Set the colour scheme
clean_colour <- rgb(81/255, 90/255, 154/255)
contam_colour <- rgb(139/255, 195/255, 205/255)

# Plot
spike_1_dists_plot <- ggplot(spike_1_dists, aes(x = Plate, y = Distance, fill = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.25), alpha = 0.8) +
  scale_fill_manual(values = rep(c(contam_colour, clean_colour), 3), name = "") +
  ylim(0, 32.5) +
  geom_signif(annotation = "***",
              y_position = 31.5, xmin = c(0.8125, 1.8125, 2.8125), xmax = c(1.187, 2.187, 3.187)) +
  labs(title = NULL, x = NULL, y = "Aitchison Distance to Positive Control") +
  theme_classic() +
  theme(legend.position = "bottom")

spike_2_dists_plot <- ggplot(spike_2_dists, aes(x = Plate, y = Distance, fill = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.25), alpha = 0.8) +
  scale_fill_manual(values = rep(c(contam_colour, clean_colour), 3)) +
  ylim(0, 32.5) +
  geom_signif(annotation = "***",
              y_position = 31.5, xmin = c(0.8125, 1.8125, 2.8125), xmax = c(1.187, 2.187, 3.187)) +
  labs(title = NULL, x = NULL, y = "Aitchison Distance to Positive Control") +
  theme_classic() +
  theme(legend.position = "none")

# Build the whole figure
staph_half <- ggarrange(staph_heatmap, spike_2_dists_plot,
                        ncol = 2,
                        labels = c("A", "B"))

staph_half <- annotate_figure(staph_half,
                              left = text_grob("Staphylococcus", size = 16, rot = 90, face = "bold.italic"))

esch_half <- ggarrange(esch_shig_heatmap, spike_1_dists_plot,
                       ncol = 2,
                       labels = c("C", "D"))

esch_half <- annotate_figure(esch_half,
                             left = text_grob("Escherichia-Shigella", size = 16, rot = 90, face = "bold.italic"))

positive_ctrl_artifact <- ggarrange(staph_half, esch_half,
                                    nrow = 2,
                                    heights = c(0.88, 1),
                                    common.legend = TRUE)

tiff("/home/ccarr/Documents/lab/implant_microbiota_study/figures/positive_ctrl_artifact.tiff", units = "in", width = 11, height = 10, res = 300)
positive_ctrl_artifact
dev.off()

# Also output the Aitchison distance matrix and distances used in comparisons
write.table(aitch_dist,
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/aitch_dist.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(spike_1_dists,
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/spike_1_dists.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(spike_2_dists,
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/spike_2_dists.txt",
            sep = "\t", col.names = NA, quote = F)
