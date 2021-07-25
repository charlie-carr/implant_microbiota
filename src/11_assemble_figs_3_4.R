# This script assembles Figures 3 and 4

# Build and output Figure 3
fig_3_top <- ggarrange(extraction_shannon_diversity_plot, extraction_envfit_and_pca[[2]],
                       ncol = 2,
                       widths = c(0.5, 1),
                       labels = c("A", "B"))

fig_3 <- ggarrange(fig_3_top, extraction_aldex2_plot,
                   nrow = 2,
                   labels = c("", "C"))

tiff("/home/ccarr/Documents/lab/implant_microbiota_study/figures/extraction.tiff", units = "in", width = 11, height = 10, res = 300)
fig_3
dev.off()

# Build and output Figure 4
fig_4_top <- ggarrange(pt_oac_shannon_diversity_plot, pt_oac_envfit_and_pca[[2]],
                       ncol = 2,
                       widths = c(0.5, 1),
                       labels = c("A", "B"))

fig_4 <- ggarrange(fig_4_top, pt_oac_aldex2_plot,
                   nrow = 2,
                   labels = c("", "C"))

tiff("/home/ccarr/Documents/lab/implant_microbiota_study/figures/pt_oac.tiff", units = "in", width = 11, height = 10, res = 300)
fig_4
dev.off()
