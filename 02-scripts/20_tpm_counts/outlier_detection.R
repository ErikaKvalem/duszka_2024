library(bigutilsr)
library(ggplot2)
#theme_set(bigstatsr::theme_bigstatsr(0.8))

## INTESTINE
counts = read.csv("/data/scratch/kvalem/projects/2024/duszka_2024/02-scripts/tables/20_tpm_counts/counts_matrix_metadata_intestine.csv")
counts <- counts[, -(1:2)]
# List of column names to remove
columns_to_remove <- c('SampleID', 'mouseID', 'group', 'subgroup', 'source', 'sequencingID','sample_type')

# Remove columns from the dataframe
X <- counts[, !(names(counts) %in% columns_to_remove)]

#PCA
pca <- prcomp(X, retx = TRUE,  scale = TRUE,tol=0.4)
predicted <-predict(pca,X)

predicted <- cbind(predicted, sequencingID = counts$sequencingID)
ggplot(data.frame(predicted))+aes(x=PC1,y=PC2,color=counts$sample_type)+geom_point(aes(size=0.5))+stat_ellipse()+stat_ellipse(level=0.8)


merged_data <- merge(predicted, counts, by = "sequencingID")
#PLOT PCA
ggplot(merged_data) +
  aes(x = PC1, y = PC2, color = sample_type) +
  geom_point(size = 0.5) +
  stat_ellipse() +
  stat_ellipse(level = 0.8) +
  geom_text(aes(label = sequencingID), hjust = 0, vjust = 0)

# Calculate outliers
data_transposed <- t(X)
pca_result <- prcomp(data_transposed, scale. = TRUE)

#
## Calculate distances from centroid
centroid <- colMeans(pca_result$x)
distances <- apply(pca_result$x, 2, function(x) sqrt(sum((x - centroid)^2)))

# Determine outliers based on threshold (more than 3 standard deviations)
outlier_threshold <- mean(distances) + 3 * sd(distances)
outliers <- which(distances > outlier_threshold)

# Print names of outlier samples
print(names(data_transposed[outliers, ]))

## LIVER
counts = read.csv("/data/scratch/kvalem/projects/2024/duszka_2024/02-scripts/tables/20_tpm_counts/counts_matrix_metadata_liver.csv")
counts <- counts[, -(1:2)]
# List of column names to remove
columns_to_remove <- c('SampleID', 'mouseID', 'group', 'subgroup', 'source', 'sequencingID','sample_type')

# Remove columns from the dataframe
X <- counts[, !(names(counts) %in% columns_to_remove)]

#PCA
pca <- prcomp(X, retx = TRUE,  scale = TRUE,tol=0.4)
predicted <-predict(pca,X)

predicted <- cbind(predicted, sequencingID = counts$sequencingID)
ggplot(data.frame(predicted))+aes(x=PC1,y=PC2,color=counts$sample_type)+geom_point(aes(size=0.5))+stat_ellipse()+stat_ellipse(level=0.8)


merged_data <- merge(predicted, counts, by = "sequencingID")
#PLOT PCA
ggplot(merged_data) +
  aes(x = PC1, y = PC2, color = sample_type) +
  geom_point(size = 0.5) +
  stat_ellipse() +
  stat_ellipse(level = 0.8) +
  geom_text(aes(label = sequencingID), hjust = 0, vjust = 0)

# Calculate outliers
data_transposed <- t(X)
pca_result <- prcomp(data_transposed, scale. = TRUE)

#
## Calculate distances from centroid
centroid <- colMeans(pca_result$x)
distances <- apply(pca_result$x, 2, function(x) sqrt(sum((x - centroid)^2)))

# Determine outliers based on threshold (more than 3 standard deviations)
outlier_threshold <- mean(distances) + 3 * sd(distances)
outliers <- which(distances > outlier_threshold)

# Print names of outlier samples
print(names(data_transposed[outliers, ]))