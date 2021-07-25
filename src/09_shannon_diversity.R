# This script computes Shannon diversity and performs Shannon diversity 
# comparisons required for Figures 3 and 4

# Compute Shannon diversity
extraction_shannon_diversity <- as.data.frame(vegan::diversity(pt_counts, index = "shannon", MARGIN = 2))
colnames(extraction_shannon_diversity) <- "Shannon_Diversity"

# Organize data for plotting
extraction_shannon_diversity <- merge(extraction_shannon_diversity, pt_meta[, "Extraction", drop = FALSE], by = 0)
extraction_shannon_diversity$Row.names <- NULL

# Plot Shannon diversity
extraction_shannon_diversity_plot <- ggplot(extraction_shannon_diversity, aes(x = Extraction, y = Shannon_Diversity, fill = Extraction)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.25), alpha = 0.8) +
  geom_signif(comparisons = list(c(" CTAB", "PowerSoil")), map_signif_level = TRUE) + 
  scale_fill_manual(values = c("#b3b3b3", "#66c2a5"), name = "") +
  labs(title = NULL, x = NULL, y = "Shannon Diversity") + 
  theme_classic() +
  theme(legend.position = "bottom")

# Export Shannon diversity data
write.table(extraction_shannon_diversity, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/extraction_shannon_diversity.txt", 
            sep = "\t", col.names = NA, quote = F)

# Compute Shannon diversity
pt_oac_shannon_diversity <- as.data.frame(vegan::diversity(pt_oac_counts, index = "shannon", MARGIN = 2))
colnames(pt_oac_shannon_diversity) <- "Shannon_Diversity"

# Organize data for plotting
pt_oac_shannon_diversity <- merge(pt_oac_shannon_diversity, pt_oac_meta[, c("Type", "Extraction")], by = 0)
pt_oac_shannon_diversity$Row.names <- NULL

# Compute p-values
wilcox.test(pt_oac_shannon_diversity[which(pt_oac_shannon_diversity$Type == "Patient" & pt_oac_shannon_diversity$Extraction == " CTAB"), "Shannon_Diversity"], 
            pt_oac_shannon_diversity[which(pt_oac_shannon_diversity$Type == "Open-air control" & pt_oac_shannon_diversity$Extraction == " CTAB"), "Shannon_Diversity"])

wilcox.test(pt_oac_shannon_diversity[which(pt_oac_shannon_diversity$Type == "Patient" & pt_oac_shannon_diversity$Extraction == "PowerSoil"), "Shannon_Diversity"], 
            pt_oac_shannon_diversity[which(pt_oac_shannon_diversity$Type == "Open-air control" & pt_oac_shannon_diversity$Extraction == "PowerSoil"), "Shannon_Diversity"])

# Plot Shannon diversity
pt_oac_shannon_diversity_plot <- ggplot(pt_oac_shannon_diversity, aes(x = Extraction, y = Shannon_Diversity, fill = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.25), alpha = 0.8) +
  geom_signif(annotations = c("**", "***"),
              y_position = 2.9, xmin = c(0.8125, 1.8125), xmax = c(1.187, 2.187)) +
  scale_fill_manual(values = c("#86adc8", "#cc453a"), name = "") +
  labs(title = NULL, x = NULL, y = "Shannon Diversity") + 
  theme_classic() +
  theme(legend.position = "bottom")

# Export Shannon diversity data
write.table(pt_oac_shannon_diversity, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/pt_oac_shannon_diversity.txt", 
            sep = "\t", col.names = NA, quote = F)
