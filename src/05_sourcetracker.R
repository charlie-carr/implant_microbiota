# This script uses SourceTracker to select patient samples. Modified from the 
# tutorial at https://github.com/danknights/sourcetracker/blob/master/example.r

# Transpose the filtered counts table produced by filtering.R
st_counts <- t(filtered_counts)

# Subset the st_meta according to the samples that remain after filtering
st_meta <- st_meta[rownames(st_counts), ]

# Select samples for training and applying the algorithm
train.ix <- which(st_meta$SourceSink == "source")
test.ix <- which(st_meta$SourceSink == "sink")

# Store the environments
enviros <- st_meta$Env

# Select alpha values by cross-validation
tune.results <- tune.st(st_counts[train.ix, ], enviros[train.ix])
alpha1 <- tune.results$best.alpha1 # 0.001
alpha2 <- tune.results$best.alpha2 # 0.01

# Train and run SourceTracker
st <- sourcetracker(st_counts[train.ix, ], enviros[train.ix])

st_results <- predict(st, st_counts[test.ix, ], alpha1 = alpha1, alpha2 = alpha2)

# Select samples that pass the 95% threshold
st_data <- as.data.frame(st_results$proportions)

st_data_95 <- st_data[which(st_data$Unknown >= 0.95), ]

