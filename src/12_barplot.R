# This script generates Figure 5 

prop_bar_plot <- function(count_data, tax, level, tax_name, extraction) {

  # Aggregate the counts data at the selected level
  counts_by_tax <- aggregate(count_data, by = list(tax[, level]), FUN = sum)
  rownames(counts_by_tax) <- counts_by_tax$Group.1
  counts_by_tax$Group.1 <- NULL

  # Convert counts to proportions
  props_by_tax <- data.frame(apply(counts_by_tax, 2, function(x){x/sum(x)}))

  # Convert the data into "long" format for ggplot
  props_by_tax$Taxon <- rownames(props_by_tax)

  long_props_by_tax <- props_by_tax %>%
    pivot_longer(!Taxon, names_to = "Sample", values_to = "Rel_Abun")

  # Add extraction information to label the plot
  long_props_by_tax$Extraction <- extraction

  # Plot
  ggplot(long_props_by_tax, aes(x = Sample, y = Rel_Abun, fill = Taxon)) +
    geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
    scale_fill_manual(values = colours, name = tax_name) +
    labs(y = "Relative Abundance", x = NULL) +
    theme_minimal() +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(legend.text = element_text(face = "italic"), legend.key.size = unit(0.5, "cm"),
          legend.key.width = unit(0.75,"cm")) +
    guides(fill = guide_legend(ncol = 1)) +
    facet_wrap(~Extraction)

}

colours <- c("steelblue3", "skyblue1", "#1C4466", "indianred1", "mediumpurple1", "#FFED6F", "ivory2", "olivedrab3", "aquamarine3", "pink", "mediumorchid3", "#C0C0C0", "#3399FF", "#339966", "#999933", "#666699", "mediumvioletred", "#666666", "#993300", "#006666", "#FFCC66", "#9999CC", "#663366", "#999966", "#D7002C")

# Create the plots for CTAB and PowerSoil samples separately
ctab_bar_plot <- prop_bar_plot(pt_counts[, rownames(pt_meta)[which(pt_meta$Extraction == " CTAB")]], pt_tax, 6, "Genus", "CTAB")
ps_bar_plot <- prop_bar_plot(pt_counts[, rownames(pt_meta)[which(pt_meta$Extraction == "PowerSoil")]], pt_tax, 6, "Genus", "PowerSoil")

# Remove y-axis labels from the second part of the plot
ps_bar_plot <- ps_bar_plot +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  labs(y = NULL)

# Assemble the full plot
full_bar_plot <- ggarrange(ctab_bar_plot, ps_bar_plot, common.legend = TRUE, legend = "right", widths = c(1, 2.975))

tiff("/home/ccarr/Documents/lab/implant_microbiota_study/figures/prop_barplot.tiff", units = "in", width = 12, height = 6, res = 300)
full_bar_plot
dev.off()
