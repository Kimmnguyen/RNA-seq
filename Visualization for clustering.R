#set_up & assigned data
library("tidyverse")
library("dplyr")
library(tidyr)
library("FNN")
library("ggplot2")
library("tibble")
library(DESeq2)

test_result= as.data.frame(res)
test_result$gene= rownames(test_result)
sample_info= as.data.frame(colData)
sample_info$hours= gsub(".*?(\\d+)h.*", "\\1",sample_info$groups)
sample_info$type= gsub("_.*","",sample_info$conditions)
trans_cts=as.data.frame( counts(dds,normalized=T))
# Filter significant genes with |log2FoldChange| > 1 and padj < 0.05
significant_genes <- subset(res, padj < 0.1)

trans_cts$gene= row.names(trans_cts)
#set candidate genes
candidate_genes <- test_result %>% 
  
  filter(padj< 0.1) %>%   # filter table
  pull(gene) %>%
  unique()

candidate_genes= as.data.frame(candidate_genes)


#### Summarise counts#### 
library(tidyverse)

trans_cts_mean <- trans_cts %>% 
  pivot_longer(cols = D0h.rep1:R4h.rep2, 
               names_to = "groups", 
               values_to = "cts") %>%
  # Join with sample_info
  full_join(sample_info, by = "groups") %>% 
  
  # Scale the cts column for each gene
  group_by(gene) %>% 
  mutate(cts_scaled = (cts - mean(cts)) / sd(cts)) %>% 
  # Group by gene and conditions
  group_by(gene, hours,type) %>%
  # Calculate the mean (scaled) cts and count of replicates
  summarise(mean_cts_scaled = mean(cts_scaled, na.rm = TRUE),  # Handle NA values
            nrep = n(), .groups = 'drop') %>% 
  ungroup()


# Print the summarized data
print(trans_cts_mean)


####Hirechal clustering####
# Create a matrix
hclust_matrix <- trans_cts %>% 
  select(-gene) %>% 
  as.matrix()

# assign rownames
rownames(hclust_matrix) <- trans_cts$gene

hclust_matrix <- hclust_matrix[candidate_genes, ]
# Scale datase
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

#caculate the distance between each gene
gene_dist <- dist(hclust_matrix)

#perform clustering
gene_hclust <- hclust(gene_dist, method = "complete")

#by trend
gene_hclust <- hclust(gene_dist, method = "ward.D2")
# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 6, col = "brown", lwd = 2)
cluster_assignments=cutree(gene_hclust, k = 6)

cutree(gene_hclust, k = 6)


# Load required libraries
library("dplyr")
library("tibble")

gene_cluster <- cutree(gene_hclust, k = 6) %>% 
  # turn the named vector into a tibble
  enframe() %>%
  # rename some of the columns
rename(gene = name, cluster = value)

head(gene_cluster)



#Visualize per cluster
trans_cts_cluster <- trans_cts %>% 
  inner_join(gene_cluster, by = "gene")

head(trans_cts_cluster)
row.names(trans_cts_cluster)= trans_cts_cluster$gene

# Load necessary library
library(ggplot2)

# Convert hours to numeric if it isn't already
trans_cts_cluster$hours <- as.numeric(trans_cts_cluster$hours)

# Plot the expression patterns by hours, type, and cluster
ggplot(trans_cts_cluster, aes(x = hours, y = mean_cts_scaled, group = gene)) +
  geom_line(alpha = 0.3) +  # Light lines for individual genes
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +  # Median trend line for each facet
  facet_grid( rows= vars(type),cols = vars(cluster)) +  # Facet by type and cluster
  theme_minimal() +  # Use a clean theme
  labs(
    x = "Hours",
    y = "Mean Scaled Counts",
    title = "Expression Pattern per Cluster",
    subtitle = "Trends across time for each cluster and condition"
  ) + 
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Bold facet labels
    axis.text = element_text(size = 10),  # Adjust axis text size
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
Heatmap(hclust_matrix, show_row_names = FALSE)

####K-means clustering####
wcss <- numeric()  # Initialize a vector to store WCSS
for (k in 1:10) {
  kmean_result <- kmeans(hclust_matrix, centers = k)
  wcss[k] <- kmean_result$tot.withinss
}
# Plot WCSS vs. number of clusters
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "WCSS", main = "Elbow Method")
#caculate the RMSE
# Example: Perform K-means clustering
set.seed(123)  # Ensure reproducibility
k <- 6 # Number of clusters
kmeans_result <- kmeans(hclust_matrix, centers = k, nstart = 25)

# Get cluster centers
centers <- kmeans_result$centers

# Assign each point to its cluster center
predicted_centers <- centers[kmeans_result$cluster, ]

# Calculate the squared error for each point
squared_errors <- rowSums((hclust_matrix - predicted_centers) ^ 2)

# Compute RMSE
rmse <- sqrt(mean(squared_errors))

# Print RMSE
print(paste("RMSE:", rmse))

#sill
library(cluster)
silhouette_scores <- silhouette(kmeans_result$cluster, dist(hclust_matrix))
avg_silhouette <- mean(silhouette_scores[, 3])  # Average silhouette score
print(paste("Average Silhouette Score:", avg_silhouette))


####finding best representer for each cluster####
#Distance to centroid
# Calculate cluster centroids
centroids <- trans_cts_cluster %>%
  group_by(cluster, hours) %>%
  summarise(mean_cts_scaled = mean(mean_cts_scaled, na.rm = TRUE), .groups = "drop")

# Calculate Euclidean distance for each gene from the cluster's centroid
distances <- trans_cts_cluster %>%
  inner_join(centroids, by = c("cluster", "hours"), suffix = c("", "_centroid")) %>%
  mutate(distance = (mean_cts_scaled - mean_cts_scaled_centroid)^2) %>%
  group_by(gene, cluster) %>%
  summarise(total_distance = sqrt(sum(distance, na.rm = TRUE)), .groups = "drop")

# Find the gene with the minimum distance in each cluster
best_representative <- distances %>%
  group_by(cluster) %>%
  slice_min(total_distance, n = 1)  # Select the gene with the smallest distance

print(best_representative)


# Visualize expression of the representative genes
trans_cts_cluster %>%
  inner_join(best_representative, by = c("gene", "cluster")) %>%
  ggplot(aes(x = as.numeric(hours), y = mean_cts_scaled, 
             group = gene, colour = factor(cluster))) +
  geom_line(size = 1.2) +
  facet_grid(type ~ cluster, scales = "free_y") +  # Facet by type (Drought/Recovery) and cluster
  theme_minimal() +
  scale_x_continuous(breaks = unique(as.numeric(trans_cts_cluster$hours))) +  # Ensure correct x-axis ticks
  labs(
    title = "Expression Patterns of Representative Genes",
    x = "Hours", y = "Mean Scaled Counts", colour = "Cluster"
  ) +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    legend.position = "right"
  )




# Install and load the necessary package
if (!require(dbscan)) install.packages("dbscan")
library(dbscan)

# Step 1: Standardize the data (optional, but recommended)
scaled_matrix <- scale(hclust_matrix)

# Step 2: Run DBSCAN
set.seed(123)  # Ensure reproducibility
dbscan_result <- dbscan(scaled_matrix, eps = 1, minPts = 5)  # Adjust eps and minPts

# Step 3: Check the cluster assignments
print(dbscan_result)

# Step 4: Add cluster information to the data
cluster_labels <- as.factor(dbscan_result$cluster)
hclust_matrix_with_labels <- data.frame(hclust_matrix, cluster = cluster_labels)

# Step 5: Visualize the clusters (if data is 2D or 3D)
library(ggplot2)

# Convert data to data frame for plotting
plot_data <- as.data.frame(scaled_matrix)
plot_data$cluster <- factor(dbscan_result$cluster)
head(plot_data)

pca_result <- prcomp(scaled_matrix)  # Perform PCA
pca_data <- as.data.frame(pca_result$x[, 1:2])  # Extract the first two principal components
colnames(pca_data) <- c("PC1", "PC2")  # Rename columns for clarity
pca_data$cluster <- factor(dbscan_result$cluster)

# Plot using PCA components
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "DBSCAN Clustering (PCA)", x = "PC1", y = "PC2")

# Plot the clusters (for the first 2 dimensions)
ggplot(plot_data, aes(V1, V2, color = cluster)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "DBSCAN Clustering", x = "Dimension 1", y = "Dimension 2")



#####decided K-means clustering vs DBSCAN####
library(cluster)

# K-means clustering
kmeans_result <- kmeans(scaled_matrix, centers = 5)
sil_kmeans <- silhouette(kmeans_result$cluster, dist(scaled_matrix))

# DBSCAN clustering
dbscan_result <- dbscan(scaled_matrix, eps = 1, minPts = 5)
sil_dbscan <- silhouette(dbscan_result$cluster, dist(scaled_matrix))

# Compute average silhouette scores
avg_sil_kmeans <- mean(sil_kmeans[, 3])
avg_sil_dbscan <- mean(sil_dbscan[, 3])

cat("Average Silhouette Score for K-means:", avg_sil_kmeans, "\n")
cat("Average Silhouette Score for DBSCAN:", avg_sil_dbscan, "\n")

