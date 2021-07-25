# This script performs and plots ALDEx2 testing required for Figures 3 and 4

# Function to run ALDEx2. This function is based on code at 
# https://github.com/ggloor/CoDa_microbiome_tutorial
run_aldex2 <- function(counts, groups, taxon, tax) {
  
  da_x <- aldex.clr(counts, groups, mc.samples = 1000, verbose = TRUE)
  da_tt <- aldex.ttest(da_x, groups, paired.test = FALSE)
  da_effect <- aldex.effect(da_x, groups, include.sample.summary = FALSE, verbose = TRUE)
  da_all <- data.frame(da_tt, da_effect, stringsAsFactors = FALSE)

  # Add the tax information to the table
  da_all <- data.frame(da_all, tax[rownames(da_all), taxon])
  colnames(da_all)[ncol(da_all)] <- taxon
  
  # Add significance information to the table
  sig <- da_all$wi.eBH < 0.05
  eff <- abs(da_all$effect) > 1
  pos <- da_all$effect > 0
  da_all <- data.frame(da_all, sig, eff, pos)
  
  da_all
  
}

# Run and plot ALDEx2 to compare extraction methodologies
extraction_aldex2 <- run_aldex2(pt_counts, pt_meta$Extraction, "Genus", pt_tax)

extraction_aldex2_plot <- ggplot(extraction_aldex2, aes(x = Genus, y = effect, fill = pos, shape = sig)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey92", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey92", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey92", size = 0.5) +
  geom_point(alpha = 0.8) +
  scale_fill_manual(values =  c("#b3b3b3", "#66c2a5"), labels = c("More abundant in CTAB-extracted samples", "More abundant in PowerSoil-extracted samples"), guide = guide_legend(override.aes = list(shape = 21))) +
  scale_shape_manual(values = c(21, 23), labels = c(expression(italic("p ") >= " 0.05"), expression(italic("p") < " 0.05"))) +
  scale_x_discrete(limits = rev) +
  labs(y = "Median Effect Size of Abundance Difference") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"))

# Export ALDEx2 results
write.table(extraction_aldex2, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/extraction_aldex2.txt", 
            sep = "\t", col.names = NA, quote = F)

# Split patient and open-air control data by DNA extraction methodology to run
# ALDEx2
ctab_pt_oac_meta <- pt_oac_meta[pt_oac_meta$Extraction == " CTAB", ]
ctab_pt_oac_counts <- pt_oac_counts[, rownames(ctab_pt_oac_meta)]

ps_pt_oac_meta <- pt_oac_meta[pt_oac_meta$Extraction == "PowerSoil", ]
ps_pt_oac_counts <- pt_oac_counts[, rownames(ps_pt_oac_meta)]

# Run ALDEx2 to compare patient samples and open-air controls
ctab_pt_oac_aldex2 <- run_aldex2(ctab_pt_oac_counts, ctab_pt_oac_meta$Type, "Genus", pt_oac_tax)
ps_pt_oac_aldex2 <- run_aldex2(ps_pt_oac_counts, ps_pt_oac_meta$Type, "Genus", pt_oac_tax)

# Combine ALDEx2 results for CTAB and PowerSoil groups
ctab_pt_oac_aldex2$Extraction = "CTAB"
ps_pt_oac_aldex2$Extraction = "PowerSoil"

pt_oac_aldex2 <- rbind(ctab_pt_oac_aldex2, ps_pt_oac_aldex2)

# Plot results of ALDEx2 comparison of patient samples and open-air controls
pt_oac_aldex2_plot <- ggplot(pt_oac_aldex2, aes(x = Genus, y = effect, fill = pos, shape = sig)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey92", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey92", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey92", size = 0.5) +
  geom_point(alpha = 0.8) +
  scale_fill_manual(values =  c("#86adc8", "#cc453a"), labels = c("More abundant in open-air controls", "More abundant in patient samples"), guide = guide_legend(override.aes = list(shape = 21))) +
  scale_shape_manual(values = c(21, 23), labels = c(expression(italic("p ") >= " 0.05"), expression(italic("p") < " 0.05"))) +
  scale_x_discrete(limits = rev) +
  labs(y = "Median Effect Size of Abundance Difference") +
  coord_flip() +
  facet_wrap(~Extraction) +
  theme_minimal() +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        axis.text.y = element_text(face = "italic"))

# Export ALDEx2 results
write.table(pt_oac_aldex2, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/pt_oac_aldex2.txt", 
            sep = "\t", col.names = NA, quote = F)
