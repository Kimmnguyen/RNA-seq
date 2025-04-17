####Reclustering####
# Calculate row means only for non-zero values
trans_cts_mean <- data.frame(
  D0h = apply(trans_cts_cluster[, c("D0h.rep1", "D0h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  D1h = apply(trans_cts_cluster[, c("D1h.rep1", "D1h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  D2h = apply(trans_cts_cluster[, c("D2h.rep1", "D2h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  D4h = apply(trans_cts_cluster[, c("D4h.rep1", "D4h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  R1h = apply(trans_cts_cluster[, c("R1h.rep1", "R1h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  R2h = apply(trans_cts_cluster[, c("R2h.rep1", "R2h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  R4h = apply(trans_cts_cluster[, c("R4h.rep1", "R4h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  gene= trans_cts_cluster$gene,
  cluster= trans_cts_cluster$cluster
  )
rownames(trans_cts_mean)= trans_cts$gene

#####For sig gene####
trans_cts_mean <- data.frame(
  D0h = apply(trans_cts_cluster_.1[, c("D0h.rep1", "D0h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  D1h = apply(trans_cts_cluster_.1[, c("D1h.rep1", "D1h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  D2h = apply(trans_cts_cluster_.1[, c("D2h.rep1", "D2h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  D4h = apply(trans_cts_cluster_.1[, c("D4h.rep1", "D4h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  R1h = apply(trans_cts_cluster_.1[, c("R1h.rep1", "R1h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  R2h = apply(trans_cts_cluster_.1[, c("R2h.rep1", "R2h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  R4h = apply(trans_cts_cluster_.1[, c("R4h.rep1", "R4h.rep2")], 1, function(x) ifelse(all(x == 0), 0, mean(x[x != 0]))),
  gene= trans_cts_cluster_.1$gene,
  cluster= trans_cts_cluster_.1$cluster
  
)
rownames(trans_cts_mean)= trans_cts_cluster_.1$gene
# Create a matrix
hclust_matrix <- log2FC_df_complete  %>% 
  as.matrix()

candidate_genes <- rownames(log2FC_df_complete)

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
gene_hclust <- hclust(gene_dist, method = "ward.D2")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 4, col = "brown", lwd = 2)
cutree(gene_hclust, k = 4)

wcss <- numeric()  # Initialize a vector to store WCSS
for (k in 1:10) {
  kmean_result <- kmeans(hclust_matrix, centers = k)
  wcss[k] <- kmean_result$tot.withinss
}
# Plot WCSS vs. number of clusters
plot(1:10, wcss, type = "b", xlab = "Number of Clusters (k)", ylab = "WCSS", main = "Elbow Method")


# Load required libraries
library(dplyr)
library(tibble)

gene_cluster <- cutree(gene_hclust, k = 6) %>% 
  # turn the named vector into a tibble
  enframe() %>%
  # rename some of the columns
  rename(gene = name, cluster = value)

head(gene_cluster)
#add gene column for merging
#trans_cts_significant= log2FC_df[rownames(log2FC_df) %in% log2FC_df$gene, ]
trans_cts_significant <- trans_cts[rownames(trans_cts) %in% trans_cts$gene, ]
trans_cts_mean= trans_cts_cluster
trans_cts_cluster <-trans_cts_mean %>% 
  inner_join(gene_cluster, by = "gene")
rownames(trans_cts_cluster)=trans_cts_mean$gene


trans_cts_cluster_matching= formatted_summary[formatted_summary$sequence_name %in% trans_cts_cluster$gene ,]
trans_cts_cluster_matching_sig= trans_cts_mean[trans_cts_mean$gene %in%formatted_summary$sequence_name, ]
trans_cts_cluster_matching=as.data.frame(trans_cts_cluster_matching)
rownames(trans_cts_cluster_matching)=trans_cts_cluster_matching$sequence_name
merged_data= left_join(trans_cts_cluster_matching,trans_cts_cluster_matching_sig,by=c("sequence_name"="gene"))
rownames(merged_data)= merged_data$sequence_name
write.csv(merged_data,"new merge_data by differential expression.csv")


######Gene select####
trans_cts_cluster_matching_0.1= trans_cts_cluster_matching_sig[trans_cts_cluster_matching_sig$gene %in% candidate_genes,]
trans_cts_cluster_matching_0.1= trans_cts_cluster_matching_0.1[,-8]
cluster1_count= trans_cts_cluster_matching_0.1[trans_cts_cluster_matching_0.1$cluster==1,]
cluster2_count= trans_cts_cluster_matching_0.1[trans_cts_cluster_matching_0.1$cluster==2,]
cluster3_count= trans_cts_cluster_matching_0.1[trans_cts_cluster_matching_0.1$cluster==3,]
cluster4_count= trans_cts_cluster_matching_0.1[trans_cts_cluster_matching_0.1$cluster==4,]
cluster5_count= trans_cts_cluster_matching_0.1[trans_cts_cluster_matching_0.1$cluster==5,]
cluster6_count= trans_cts_cluster_matching_0.1[trans_cts_cluster_matching_0.1$cluster==6,]


cluster_counts= list(cluster1_count,cluster2_count,cluster3_count,cluster4_count,cluster5_count,cluster6_count)
# Inspect the results
cluster_counts


#####Gene Filter####
# Initialize a list to store gene statistics for each cluster
gene_stats_list <- list()

# Process each filtered cluster to compute statistics
for (i in 1:6) {
  # Get the current cluster count data
  cluster_count <- cluster_counts[[i]]  # Correct indexing without "_count"
  
  # Compute the mean and variance for each gene
  means <- apply(cluster_count[, 1:7], 1, mean)  # Select only numeric columns (D0h to R4h)
  vars <- apply(cluster_count[, 1:7], 1, var)
  
  # Compute coefficient of variation (CV)
  cv <- sqrt(vars) / means
  
  # Create a data frame to store the statistics
  gene_stats <- data.frame(
    mean = means,
    variance = vars,
    cv = cv,
    row.names = rownames(cluster_count)
  )
  
  # Sort genes by variance in decreasing order and select top 80% most variable genes
  vars_sorted <- sort(vars, decreasing = TRUE)
  top_var <- names(vars_sorted)[1:100]
  top_gene_stats <- gene_stats[top_var, ]
  
  # Store the top gene stats in the list
  gene_stats_list[[paste0("cluster", i, "_gene_stats")]] <- top_gene_stats
}

# Inspect the results
gene_stats_list


#####PCA####
# Initialize a list to store results for each cluster
pca_results_list <- list()
leading_genes_list <- list()
top_correlations_list <- list()

# Process each cluster
for (i in 1:length(cluster_counts)) {
  # Get the current cluster count data
  cluster_count <- cluster_counts[[i]] # Access the cluster data directly
  
  # Perform PCA on the top variable genes
  pca_result <- prcomp(t(cluster_count[, 1:7]), center = TRUE, scale. = TRUE) # Transpose for PCA
  
  # Calculate percentage variance explained
  PC_sd <- setNames(pca_result$sdev, paste0("PC", 1:length(pca_result$sdev)))
  PC_var_expl <- (PC_sd^2) / sum(PC_sd^2) * 100
  
  # Store the PCA result
  pca_results_list[[paste0("cluster", i)]] <- list(
    pca = pca_result,
    variance_explained = PC_var_expl
  )
  
  # Identify leading genes based on PC1 contributions
  leading_genes <- pca_result$rotation
  leading_PC1 <- sort(leading_genes[, 1], decreasing = TRUE)
  leadgene <- names(leading_PC1)
  leading_genes_list[[paste0("cluster", i)]] <- leadgene
  
  # Calculate cluster profile as the average across columns
  cluster_profile <- colMeans(cluster_count[, 1:7])
  
  # Compute correlations of each gene with the cluster profile
  gene_correlations <- apply(cluster_count[, 1:7], 1, function(gene) cor(gene, cluster_profile))
  
  # Get the top genes based on correlation values
  top_genes <- names(sort(gene_correlations, decreasing = TRUE))
  top_genes <- top_genes[top_genes %in% leadgene]
  top_correlations <- sort(gene_correlations, decreasing = TRUE)
  
  # Store the top correlations
  top_correlations_list[[paste0("cluster", i)]] <- list(
    top_genes = top_genes,
    correlations = top_correlations
  )
}


# Set up a 2x3 layout for plotting all clusters
par(mfrow = c(2, 3))

# Loop through each cluster to compute and plot
for (i in 1:length(cluster_counts)) {
  # Get the current cluster name and data
  cluster_name <- names(cluster_counts)[i]
  cluster_count <- cluster_counts[[i]]
  
  # Calculate the cluster profile as the average across columns
  cluster_profile <- colMeans(cluster_count)
  
  # Compute correlations for each gene with the cluster profile
  gene_correlations <- apply(cluster_count, 1, function(gene) cor(gene, cluster_profile))
  
  # Get the top 5 genes based on correlation values
  top_genes <- names(sort(gene_correlations, decreasing = TRUE))[1:5]
  
  # Set x-axis range with extra space for labels
  xlim_range <- c(0, length(cluster_profile) + 10)
  
  # Plot the cluster profile
  plot(cluster_profile, type = "l", col = "blue", 
       main = paste("Cluster Profile:", cluster_name), 
       ylab = "Expression", xlab = "Conditions", 
       ylim = c(0, max(c(cluster_profile, sapply(top_genes, function(g) max(cluster_count[g, ]))))), 
       xlim = xlim_range)
  
  # Overlay lines for each of the top 5 genes
  for (j in seq_along(top_genes)) {
    gene <- top_genes[j]
    gene_profile <- as.numeric(cluster_count[gene, ])
    
    # Plot the line for each gene
    lines(gene_profile, col = "red", lty = 2)
    
    # Calculate label position
    label_x_position <- length(gene_profile)
    label_y_position <- gene_profile[length(gene_profile)] + (j %% 2) * 5  # Offset for clarity
    
    # Add the text label at the end of each gene line
    text(x = label_x_position, y = label_y_position, labels = gene, pos = 4, col = "red", cex = 0.8)
  }
}

# Set up a 2x3 layout for plotting all clusters
par(mfrow = c(2, 3))

# Loop through each cluster
for (i in 1:length(cluster_counts)) {
  # Get the current cluster name and data
  cluster_name <- paste0("Cluster", i)
  cluster_count <- cluster_counts[[i]]
  
  # Split data into Drought (D) and Recovery (R) conditions
  drought_cols <- grep("^D", colnames(cluster_count))
  recovery_cols <- grep("^R", colnames(cluster_count))
  
  # Calculate the combined cluster profile (Drought + Recovery)
  cluster_profile <- colMeans(cluster_count)
  
  # Compute correlations for each gene with the combined cluster profile
  gene_correlations <- apply(cluster_count, 1, function(gene) cor(gene, cluster_profile))
  
  # Get the top 5 genes based on correlation values
  top_genes <- names(sort(gene_correlations, decreasing = TRUE))[1:5]
  
  # Calculate profiles for Drought and Recovery conditions
  cluster_profile_drought <- colMeans(cluster_count[, drought_cols])
  cluster_profile_recovery <- colMeans(cluster_count[, recovery_cols])
  
  # Set x-axis range with extra space for labels
  xlim_range <- c(0, length(cluster_profile_drought) + length(cluster_profile_recovery) + 10)
  
  # Plot the cluster profiles
  plot(cluster_profile_drought, type = "l", col = "blue", 
       main = paste("Top 5 Genes in Cluster", i), 
       ylab = "Expression", xlab = "Conditions", 
       ylim = range(
         c(
           cluster_profile_drought, 
           cluster_profile_recovery, 
           sapply(top_genes, function(g) cluster_count[g, ])
         ),
         na.rm = TRUE), 
       xlim = xlim_range)
  lines(seq_along(cluster_profile_recovery) + length(cluster_profile_drought), 
        cluster_profile_recovery, col = "green", lty = 1)
  
  # Add legend
  legend("topright", legend = c("Drought", "Recovery"), col = c("blue", "green"), lty = 1, bty = "n")
  
  # Overlay lines for the top 5 genes
  for (j in seq_along(top_genes)) {
    gene <- top_genes[j]
    
    # Get gene profiles for Drought and Recovery conditions
    gene_profile_drought <- as.numeric(cluster_count[gene, drought_cols])
    gene_profile_recovery <- as.numeric(cluster_count[gene, recovery_cols])
    
    # Plot the line for each gene under Drought
    lines(seq_along(gene_profile_drought), gene_profile_drought, col = "red", lty = 2)
    
    # Plot the line for each gene under Recovery
    lines(seq_along(gene_profile_recovery) + length(gene_profile_drought), 
          gene_profile_recovery, col = "orange", lty = 2)
    
    # Add labels at the end of the Recovery line
    label_x_position <- length(gene_profile_drought) + length(gene_profile_recovery)
    label_y_position <- gene_profile_recovery[length(gene_profile_recovery)] + (j %% 2) * 5  # Offset for clarity
    
    # Add the text label at the end of the Recovery line
    text(x = label_x_position, y = label_y_position, labels = gene, pos = 4, col = "red", cex = 0.8)
  }
}


#####Seperate by Rep####
# Set up a 2x3 layout for plotting all clusters
par(mfrow = c(2, 3))

# Loop through each cluster
for (i in 1:length(cluster_counts)) {
  # Get the current cluster name and data
  cluster_name <- names(cluster_counts)[i]
  cluster_count <- cluster_counts[[cluster_name]]
  
  # Split data into replicates
  rep1_cols <- grep("rep1$", colnames(cluster_count))
  rep2_cols <- grep("rep2$", colnames(cluster_count))
  
  # Calculate cluster profiles for each replicate
  cluster_profile_rep1 <- colMeans(cluster_count[, rep1_cols])
  cluster_profile_rep2 <- colMeans(cluster_count[, rep2_cols])
  
  # Combine profiles for correlation calculations
  combined_profile <- c(cluster_profile_rep1, cluster_profile_rep2)
  
  # Compute correlations for each gene with the combined profile
  gene_correlations <- apply(cluster_count, 1, function(gene) cor(gene, combined_profile))
  
  # Get the top 5 genes based on correlation values
  top_genes <- names(sort(gene_correlations, decreasing = TRUE))[1:5]
  
  # Set x-axis range with extra space for labels
  xlim_range <- c(0, length(cluster_profile_rep1) + length(cluster_profile_rep2) + 10)
  
  # Plot the cluster profiles
  plot(cluster_profile_rep1, type = "l", col = "blue", 
       main = paste("Top 5 Genes in Cluster", i), 
       ylab = "Expression", xlab = "Conditions (Replicates)", 
       ylim = c(0, max(c(cluster_profile_rep1, cluster_profile_rep2, 
                         sapply(top_genes, function(g) max(cluster_count[g, ]))))), 
       xlim = xlim_range)
  lines(seq_along(cluster_profile_rep2) + length(cluster_profile_rep1), 
        cluster_profile_rep2, col = "green", lty = 1)
  
  # Add legend
  legend("topright", legend = c("Rep1", "Rep2"), col = c("blue", "green"), lty = 1, bty = "n")
  
  # Overlay lines for the top 5 genes
  for (j in seq_along(top_genes)) {
    gene <- top_genes[j]
    
    # Get gene profiles for each replicate
    gene_profile_rep1 <- as.numeric(cluster_count[gene, rep1_cols])
    gene_profile_rep2 <- as.numeric(cluster_count[gene, rep2_cols])
    
    # Plot the line for each gene under Rep1
    lines(seq_along(gene_profile_rep1), gene_profile_rep1, col = "red", lty = 2)
    
    # Plot the line for each gene under Rep2
    lines(seq_along(gene_profile_rep2) + length(gene_profile_rep1), 
          gene_profile_rep2, col = "orange", lty = 2)
    
    # Add labels at the end of the Rep2 line
    label_x_position <- length(gene_profile_rep1) + length(gene_profile_rep2)
    label_y_position <- gene_profile_rep2[length(gene_profile_rep2)] + (j %% 2) * 5  # Offset for clarity
    
    # Add the text label at the end of the Rep2 line
    text(x = label_x_position, y = label_y_position, labels = gene, pos = 4, col = "red", cex = 0.8)
  }
}

#####Pulling out cluster data####
# Initialize a list to store filtered data for the top 10% most variable genes in each cluster
top_10_percent_data <- list()
#select top 10 gene
top_10_gene= list()
# Loop through each cluster
for (i in 1:6) {
  # Get the current cluster name and top 10% variable genes
  cluster_name <- names(cluster_counts)[i]
  cluster_count <- cluster_counts[[i]]
  top_var <- rownames(gene_stats_list[[paste0("cluster", i, "_gene_stats")]])
  
  # Filter the cluster data for the top 10% most variable genes
  filtered_data <- cluster_count[rownames(cluster_count) %in% top_var, ]
  
  # Store the filtered data
  top_10_percent_data[[i]] <- filtered_data
}

cluster1_tab= top_10_percent_data[[1]]
cluster2_tab= top_10_percent_data[[2]]
cluster3_tab= top_10_percent_data[[3]]
cluster4_tab= top_10_percent_data[[4]]
cluster5_tab= top_10_percent_data[[5]]
cluster6_tab= top_10_percent_data[[6]]
topgene_cluster <- do.call(rbind, list(cluster1_tab,cluster2_tab, cluster3_tab,cluster4_tab,cluster5_tab,cluster6_tab))
topgene_cluster$gene= rownames(topgene_cluster)
topgene_motifocuur= merged_data[merged_data$sequence_name %in% topgene_cluster$gene,]
write.csv(topgene_motifocuur,"top 100 significant gene in 6 cluster by differential expression.csv")

#top 10 gene
top_10_gene= list()
for (i in 1:6) {
  # Get the current cluster name and top 10% variable genes
  cluster_name <- names(cluster_counts)[i]
  cluster_count <- cluster_counts[[i]]
  top_var <- rownames(gene_stats_list[[paste0("cluster", i, "_gene_stats")]])
  
  # Filter the cluster data for the top 10% most variable genes
  filtered_data <- cluster_count[rownames(cluster_count) %in% top_var, ]
  
  # Store the filtered data
  top_10_gene[[i]] <- filtered_data
}
top10_cluster1= top_10_gene[[1]]
top10_cluster2= top_10_gene[[2]]
top10_cluster3= top_10_gene[[3]]
top10_cluster4= top_10_gene[[4]]
top10_cluster5= top_10_gene[[5]]
top10_cluster6= top_10_gene[[6]]
top10gene_cluster <- do.call(rbind, list(top10_cluster1,top10_cluster2, top10_cluster3,top10_cluster4,top10_cluster5,top10_cluster6))
top10gene_cluster$gene= row.names(top10gene_cluster)

#merged file
#make gene column for merging
top10gene= merged_data[merged_data$sequence_name %in% top10gene_cluster$gene,]
top_10_gene=left_join(top10gene,top10gene_cluster, by=c("sequence_name"="gene"))
write.csv(top10gene,"top 10 gene in 6 cluster significant.csv")

# check for duplicated#
topgene_cluster$gene= rownames(topgene_cluster)
table(duplicated(topgene_cluster$gene))

