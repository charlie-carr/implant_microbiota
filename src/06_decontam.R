# This scripts executes the decontam pipeline

# Build the phyloseq object with filtered data
count_mat <- as.matrix(filtered_counts)
tax_mat <- as.matrix(filtered_tax)

count_phylo <- otu_table(count_mat, taxa_are_rows = TRUE)
tax_phylo <- tax_table(tax_mat)
meta_phylo <- sample_data(filtered_meta)

physeq <- phyloseq(count_phylo, tax_phylo, meta_phylo)

# Search for contaminants by prevalence 
sample_data(physeq)$neg_ctrls <- sample_data(physeq)$Type %in% c("Blank", "CTAB negative control", "Open-air control")
contam_by_prev_0.2 <- isContaminant(physeq, method = "prevalence", neg = "neg_ctrls", threshold = 0.2)
rownames(contam_by_prev_0.2)[which(contam_by_prev_0.2$contaminant)] 
# Output: "SV_2"

# Export decontam by prevalence results
write.table(contam_by_prev_0.2, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/contam_by_prev_0.2.txt", 
            sep = "\t", col.names = NA, quote = F)

# Search for contaminants by DNA concentration
conc_meta_phylo <- sample_data(conc_data)

conc_physeq <- phyloseq(count_phylo, tax_phylo, conc_meta_phylo)

contam_by_freq_0.2 <- isContaminant(conc_physeq, method = "frequency", conc = "DNA_Concentration", threshold = 0.2)
rownames(contam_by_freq_0.2)[which(contam_by_freq_0.2$contaminant)] 
# Output: "SV_9"

# Export decontam by DNA concentration results
write.table(contam_by_freq_0.2, 
            file = "/home/ccarr/Documents/lab/implant_microbiota_study/data/contam_by_freq_0.2.txt", 
            sep = "\t", col.names = NA, quote = F)
