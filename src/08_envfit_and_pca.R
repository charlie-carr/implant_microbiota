# This script generates PCA biplots and performs envfit testing required for 
# Figures 3 and 4

# Function to run envfit and make PCA biplot to test and visualize associations
# between metadata and ordination
envfit_and_pca_biplot <- function(count_data, metadata, make_biplot, colour_scheme = NULL) {

  czm <- cmultRepl(t(count_data), label = 0, method = "CZM")
  clr <- t(apply(czm, 1, function(x){log(x) - mean(log(x))}))
  pca <- prcomp(clr)

  # Run envfit
  envfit_result <- envfit(pca, data.frame(metadata))

  if (make_biplot) {
    
    # Make PCA biplot
    sample_positions <- data.frame(pca[["x"]], metadata)
    sv_positions <- data.frame(pca[["rotation"]])

    plot <- ggplot(sample_positions, aes(x = PC1, y = PC2)) +
      stat_ellipse(level = 0.95, aes(fill = metadata, color = metadata), geom = "polygon", alpha = 0.05, linetype = "dashed") +
      geom_point(aes(fill = metadata), pch = 21, alpha = 0.8) +
      scale_color_manual(values = colour_scheme, name = "") +
      scale_fill_manual(values = colour_scheme, name = "") +
      geom_segment(data = sv_positions, aes(x = 0, y = 0, xend = 15 * PC1, yend = 15 * PC2),
                   arrow = arrow(length = unit(1/2, 'picas')),
                   color = "black",
                   size = 0.1,
                   alpha = 0.5) +
      geom_text(data = sv_positions, aes(x = 15 * PC1, y = 15 * PC2), check_overlap = TRUE,
                label = rownames(sv_positions),
                color = "black",
                size = 2.5,
                fontface = "italic") +
      labs(x = paste("PC1: ", round(pca$sdev[1]^2/sum(pca$sdev^2),3)),
           y = paste("PC2: ", round(pca$sdev[2]^2/sum(pca$sdev^2),3))) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size = 1, fill = NA),
            legend.position = "bottom")

    list(envfit_result, plot)

  } else {

    list(envfit_result)

  }

}

# Perform CTAB v PowerSoil and OAC v patient sample comparisons, and test
# associations with clinical data
extraction_envfit_and_pca <- envfit_and_pca_biplot(pt_counts_by_genus, pt_meta[, "Extraction"], TRUE, c("#b3b3b3", "#66c2a5"))

pt_oac_envfit_and_pca <- envfit_and_pca_biplot(pt_oac_counts_by_genus, pt_oac_meta[, "Type"], TRUE, c("#86adc8", "#cc453a"))

clinical_envfit <- envfit_and_pca_biplot(pt_counts_by_genus, clinical_meta, FALSE)

# Export results
saveRDS(extraction_envfit_and_pca, "~/Documents/lab/implant_microbiota_study/data/extraction_envfit_and_pca.RData")

saveRDS(pt_oac_envfit_and_pca, "~/Documents/lab/implant_microbiota_study/data/pt_oac_envfit_and_pca.RData")

saveRDS(clinical_envfit, "~/Documents/lab/implant_microbiota_study/data/clinical_envfit.RData")
