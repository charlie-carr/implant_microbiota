# This script removes contaminants and re-organizes data before further analysis

# Create counts and metadata tables with only SourceTracker-selected patient samples
pt_counts <- filtered_counts[, rownames(st_data_95)]
pt_meta <- filtered_meta[rownames(st_data_95), ]

# Combine SourceTracker-selected patient data with open-air control data before 
# additional filtering
pt_oac_names <- c(rownames(st_data_95), rownames(filtered_meta)[which(filtered_meta$Type == "Open-air control")])
pt_oac_counts <- filtered_counts[, pt_oac_names]
pt_oac_meta <- filtered_meta[pt_oac_names, ]

# Update counts and tax tables to remove contaminants identified by decontam
pt_counts <- pt_counts[!rownames(pt_counts) %in% c("SV_2", "SV_9"), ]
pt_tax <- filtered_tax[rownames(pt_counts), ]

pt_oac_counts <- pt_oac_counts[!rownames(pt_oac_counts) %in% c("SV_2", "SV_9"), ]
pt_oac_tax <- filtered_tax[rownames(pt_oac_counts), ]

# Collapse ASVs by genus
pt_counts_by_genus <- aggregate(pt_counts, by = list(pt_tax$Genus), FUN = sum)
rownames(pt_counts_by_genus) <- pt_counts_by_genus$Group.1
pt_counts_by_genus$Group.1 <- NULL

pt_oac_counts_by_genus <- aggregate(pt_oac_counts, by = list(pt_oac_tax$Genus), FUN = sum)
rownames(pt_oac_counts_by_genus) <- pt_oac_counts_by_genus$Group.1
pt_oac_counts_by_genus$Group.1 <- NULL
