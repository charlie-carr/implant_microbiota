# This script pseudorandomly selects samples to examine with
# species-specific PCR

# Group the counts by genus
counts_by_genus <- aggregate(filtered_counts, by = list(filtered_tax$Genus), FUN = sum)
rownames(counts_by_genus) <- counts_by_genus$Group.1
counts_by_genus$Group.1 <- NULL

# Aggregate these counts to determine reads/genus
total_counts_by_genus <- as.data.frame(rowSums(counts_by_genus))
colnames(total_counts_by_genus) <- "Total Reads"
# Output:
#                                                    Total Reads
# Achromobacter                                            17606
# Acinetobacter                                             3667
# Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium         463
# Bacillus                                                  6604
# Brachybacterium                                           1252
# Bradyrhizobium                                             461
# Brevundimonas                                              746
# Cloacibacterium                                           1250
# Corynebacterium                                            480
# Cutibacterium                                             3976
# Delftia                                                   3234
# Enhydrobacter                                              777
# Enterococcus                                              2897
# Escherichia-Shigella                                   3681925
# Hydrotalea                                                 479
# Jeotgalicoccus                                             383
# Klebsiella                                                 644
# Lactobacillus                                             2461
# Ottowia                                                   3415
# Phyllobacterium                                           1960
# Pseudomonas                                               8752
# Sporolactobacillus                                        6960
# Staphylococcus                                         4290874
# Stenotrophomonas                                          8605
# Terrilactibacillus                                         909
# Thauera                                                    275
# Xanthomonas                                                389

# With the target taxa selected, select samples to test
sample_select <- function(counts_table, taxon, n_zero, n_nonzero) {
  
  # Select the counts for the genus of interest
  abbrev_counts <- counts_table[grep(taxon, rownames(counts_table)), ]
  
  # Sort the samples by number of reads (zero or nonzero)
  zeros <- colnames(counts_table)[apply(abbrev_counts, 1, function(x) {x == 0})]
  nonzeros <- colnames(counts_table)[apply(abbrev_counts, 1, function(x) {x != 0})]
  
  # Grab the appropriate numbers of each type of sample
  zero_inds <- sample(1:length(zeros), n_zero, replace = FALSE) 
  nonzero_inds <- sample(1:length(nonzeros), n_nonzero, replace = FALSE) 
  
  # Print results
  print(paste("Samples with zero", taxon, "reads:"))
  print(zeros[zero_inds])
  print(paste("Samples with more than zero", taxon, "reads:")) 
  print(nonzeros[nonzero_inds])
  
  # Check these to verify the output:
  View(abbrev_counts[, zeros[zero_inds]])
  View(abbrev_counts[, nonzeros[nonzero_inds]])
  
}

sample_select(counts_by_genus, "Staphylococcus", 13, 13)
# "Samples with zero Staphylococcus reads:"
# "MN2_3"     "MN30_4"    "MN36_4"    "MN15_OA11" "MN29_6"    "MN10_7"    "MN12_8"    "MN11_8"   
# "4AC"       "MN16_OA7"  "MN36_1"    "MN34_6"    "MN8_3"    

# "Samples with more than zero Staphylococcus reads:"
# "MN34_1"    "MN8_2"     "MN26_2"    "MN10_10"   "MN19_5"    "3_3"       "MN1_9"     "MN32_4"   
# "MN24_OA11" "MN4_4"     "MN33_4"    "MN12_5"    "MN9_6" 

sample_select(counts_by_genus, "Cutibacterium", 26, 0) 
# "Samples with zero Cutibacterium reads:"
# "MN18_OA5"       "82919_OA"       "MN33_OA7"       "MN7_5"          "MN9_OA11"       "MN29_OA8"      
# "MN11_4"         "MN8_7"          "MN15_7"         "MN9_2"          "MN27_4"         "MN14_3"        
# "MN20_3"         "PCR_blank_2F10" "MN10_1"         "4_1"        "MN15_9"         "MN13_3"        
# "MN1_2"          "MN23_2"         "6AC"         "DNA_blank_3B9"  "MN14_1"         "MN27_9" 
# "MN25_1"         "4_6" 

# Samples with more than zero Cutibacterium reads (selected manually): MN18_1 and MN18_3
