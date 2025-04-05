####Package used####
library(DESeq2)
library(tidyverse)
library(airway)
library(readr)


#####open data
TPM10 <- read.delim("TPM10.txt", header=TRUE)
counts_data= TPM10
#Col data
colData=read.csv("samples.csv")


#####Data cleaning#####
#naming the columns name & remove the first row in counts data
TPM10[,-1]= lapply(TPM10[,-1],as.integer)
rownames(TPM10)=TPM10[,1]
TPM10= TPM10[,-1]
#convert to numeric
summary(counts_data)


#col_data
rownames(colData)= colData[,3]
colData=colData[,-1]
#Since the counts data doesn't have R0, needs to remove it
colData= colData[c(-9,-10),]



#Matching check
#row in col data should found in column of counts_ data
all(colnames(counts_data)%in% rownames(colData)) #should return TRUE
#check order
all(colnames(counts_data)==rownames(colData)) #return TRUE

conditions = factor(c("Drought_0h", "Drought_1h",
                   "Drought_2h", "Drought_4h",
                      "Recover_1h", "Recover_2h",
                      "Recover_4h"))
####construct DESEQ object####
dds=DESeqDataSetFromMatrix(countData =counts_data,
                           colData = colData,
                           design= ~conditions)

dds
####optional####
keep=rowSums(counts(dds)) >=10
dds= dds[keep,]
# set factor level
levels(dds$conditions)
dds$conditions= relevel(dds$conditions, ref = "Drought_0h")
#NOTE: collapse the replicate
####RUN DESEQ####
dds= DESeq(dds)
####result####
res= results(dds)
res
summary(res)
res[order(res$log2FoldChange),]
res_0.05= results(dds,alpha= 0.05)
summary(res_0.05)
#####Extract all pairwise ####
#automated all the extraction
reference_level <- "Drought_0h"
results_list <- list()

for (level in setdiff(levels(dds$conditions), reference_level)) {
  results_list[[paste(level, "vs", reference_level)]] <- results(dds, contrast = c("conditions", level, reference_level))
}

#Pairwise result
pairwise_results <- list()

for (i in seq_along(levels(dds$conditions))) {
  for (j in seq_along(levels(dds$conditions))) {
    if (i < j) { # Ensure unique pairs
      condition1 <- levels(dds$conditions)[i]
      condition2 <- levels(dds$conditions)[j]
      pairwise_results[[paste(condition1, "vs", condition2)]] <- results(dds, contrast = c("conditions", condition2, condition1))
    }
  }
}
#save the result 
for (name in names(results_list)) {
  write.csv(as.data.frame(results_list[[name]]), file = paste0(name, ".csv"))
}

#Combine the file
combined_results <- do.call(rbind, lapply(names(pairwise_results), function(name) {
  df <- as.data.frame(pairwise_results[[name]])
  df$comparison <- name
  return(df)
}))
#####Choose significant gene
significant_gene= rownames(combined_results [!is.na(combined_results$padj) & combined_results$padj < 0.1,])
significant_gene= unique(significant_gene)
#####Significant level#####
significant_results <- combined_results[rownames(combined_results) %in% significant_gene, ]
log2foldchange_gene= data.frame(
    gene= rownames(significant_results),
    log2FoldChange= significant_results$log2FoldChange,
    condition= significant_results$comparison
)
# Reshape the data into a wide format
gene_table <- reshape(
  data = log2foldchange_gene,
  idvar = "gene",
  timevar = "condition",
  direction = "wide"
)

# Add the Drought_0h reference column and set it to 0
gene_table$Drought_0h <- 0

# Reorder columns to place the reference first
gene_table <- gene_table[, c("gene", "Drought_0h", setdiff(colnames(gene_table), c("gene", "Drought_0h")))]

# View the resulting table
print(gene_table)


####Plot####
# Generate MA plot
plotMA(res, main="MA Plot of 7976 genes", ylim=c(-5, 5))
# Add labels to significantly differentially expressed genes (adjust as needed)
sigGenes <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 2)  # Example threshold
text(res$baseMean[sigGenes], res$log2FoldChange[sigGenes], 
     labels=rownames(res)[sigGenes], pos=3, cex=0.6, col="red")


####Modification####
#remove genes with zero counts
res[which(res$baseMean == 0),] %>% 
  data.frame() %>% 
  View()
ntd= normTransform(dds)
significant_thres= 0.1
significant_result= res[which(res$padj<=0.1),]
upregulated= rownames(significant_result[significant_result$log2FoldChange>0,])
downregulated= rownames(significant_result[significant_result$log2FoldChange<0,])

######differential expression#######
# Initialize a list to store log2FoldChange vectors
log2FC_list <- list()

# Loop through all conditions (including D0h as a reference with log2FoldChange = 0)
conditions <- levels(dds$conditions)

for (cond in conditions) {
  if (cond == "Drought_0h") {
    # Add Drought_0h with log2FoldChange = 0 for all genes
    log2FC_list[[cond]] <- rep(0, nrow(dds))
    names(log2FC_list[[cond]]) <- rownames(dds)
  } else {
    # Extract differential expression results for other conditions
    res <- results(dds, contrast = c("conditions", cond, "Drought_0h"))
    
    # Filter significant genes (padj <= 0.1)
    #significant_res <- res[which(res$padj <= 0.1), ]
    significant_res <- res
    
    # Store log2FoldChange with gene names
    log2FC <- significant_res$log2FoldChange
    names(log2FC) <- rownames(significant_res)
    log2FC_list[[cond]] <- log2FC
  }
}

# Combine log2 fold changes into a single data frame
log2FC_df <- do.call(cbind, lapply(log2FC_list, function(x) {
  # Ensure all gene names are included
  full_vector <- rep( 0,nrow(dds))
  names(full_vector) <- rownames(dds)
  full_vector[names(x)] <- x
  return(full_vector)
}))

# Convert to a data frame and set column names
log2FC_df <- as.data.frame(log2FC_df)
log2FC_df= log2FC_df[rownames(log2FC_df) %in% significant_gene,]
colnames(log2FC_df) <- conditions
log2FC_df$gene= rownames(log2FC_df)
significant_res$gene= row.names(significant_res)
# Filter out genes with any NA values
# Create a logical vector indicating whether each gene in log2FC_df is significant
log2FC_df_complete <- log2FC_df[log2FC_df$gene %in%rownames(significant_res),]
log2FC_df_complete=log2FC_df_complete[,-8]

# Convert to a data frame and set column names
log2FC_df <- as.data.frame(log2FC_df)
colnames(log2FC_df) <- conditions
# Filter out genes with any NA values
log2FC_df_complete <- log2FC_df[complete.cases(log2FC_df), ]
nrow(log2FC_df_complete)

# View significant genes
head(significant_genes)

# View results
head(lrt_res)

#Shrinkage method
resultsNames(dds)
resLFC= lfcShrink()

####Bug Fix########Bug Fix####data.frame()
str(counts_data)  # Check structure of the data
sapply(counts_data, class)  # Verify the class of each column
colData$conditions= gsub(" ","_",colData$conditions)
colData$conditions <- as.factor(colData$conditions)

####Normalized the dataset####
ntd= normTransform(dds)
library("pheatmap")
select= order(rowMeans(counts(dds,normalized=T)),decreasing = FALSE)[1:20]
df_anotate= as.data.frame(colData(dds)[,c("conditions","groups")])

gene_names= rownames(ntd)
rownames(assay(ntd))=gene_names
# Specify the order of the conditions
desired_order <- c("Drought_0h", "Drought_1h", "Drought_2h", "Drought_4h", "Recover_1h", "Recover_2h", "Recover_4h")
pheatmap(
  assay(ntd)[select,order(ntd$groups)],  # Ordered expression levels of selected genes
  cluster_rows = TRUE,   # Don't cluster genes
  show_rownames = TRUE,   # Show gene names
  cluster_cols = F,   # Don't cluster samples
  annotation_col = df_anotate,  # Add condition annotations
  main = "Gene Expression Heatmap by Ordered Conditions"
)

####transformation& visualiztion####
vsd= vst(dds,blind=FALSE)
rld= rlog(dds,blind= FALSE)
head(assay(vsd),3)
ntd= normTransform(dds)
#heatmap

sample_dis= dist(t(assay(ntd)))
sampleDis_matrix= as.matrix(sample_dis)
# Set row and column names based on your sample identifiers
rownames(sampleDis_matrix) <- colnames(normalized_counts)
colnames(sampleDis_matrix) <- colnames(normalized_counts)

