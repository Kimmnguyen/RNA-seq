# Load libraries
library(GenomicRanges)
library(vcfR)
# Load gene coordinates
gene_coordinates <- read.table("gene_coordinates.bed", header = FALSE, col.names = c("chromosome", "start", "end", "gene"))
# Convert gene coordinates to GRanges
genes_gr <- GRanges(
  seqnames = gene_coordinates$chromosome,
  ranges = IRanges(start = gene_coordinates$start, end = gene_coordinates$end),
  gene = gene_coordinates$gene
)

# Load the filtered VCF
vcf <- read.vcfR("filtered_genes.vcf.gz")

# Convert fixed fields to a data frame
fixed_fields <- as.data.frame(vcf@fix)

# View the column names
colnames(fixed_fields)

#Fixed the format data
library(reshape2)

# Reshape the matrix directly
genotype_long <- melt(vcf@gt, varnames = c("Row", "accession"), value.name = "value")

# Rename the first column to "FORMAT" (assumes it's the first column in your matrix)
genotype_long$FORMAT <- vcf@gt[genotype_long$Row, "FORMAT"]

# Drop the "Row" column (if unnecessary)
genotype_long <- genotype_long[, c("FORMAT", "accession", "value")]

# View the result
head(genotype_long)


# Convert genotype data to a data frame
genotype_data <- as.data.frame(vcf@gt)

#make it as one column
library(tidyr)

# Convert the genotype data into the desired format
genotype_long <- genotype_data %>%
  pivot_longer(
    cols = -FORMAT,        # Exclude the FORMAT column
    names_to = "accession", # New column for accession names
    values_to = "value"     # New column for genotype values
  )

# Combine the two data frames
combined_data <- cbind(fixed_fields, genotype_data)

# View the result
head(combined_data)



# Add accession IDs to the variants
variant_to_gene_with_info <- data.frame(
  variant_position = paste(seqnames(vcf_gr), start(vcf_gr), sep = ":")[variant_to_gene$variant_index],
  gene = variant_to_gene$gene
)

# Print the result
head(variant_to_gene_with_info)

variants <- extract.gt(vcf)
gene_vcf <- merge(
  gene_coordinates,
  vcf_variants,
  by.y = c("chromosome", "position"),  # Ensure `position` exists in `vcf_variants`
  by.x = c("chromosome", "start"),
  all.x = TRUE
)

# Save the merged data
write.csv(gene_vcf, "merged_variants_with_genes.csv", row.names = FALSE)
print(head(variants))

# Merge by chromosome and position
gene_vcf <- merge(
  gene_coordinates,
  vcf_variants,
  by.y = c("chromosome", "position"),
  by.x = c("chromosome", "start"),  # Adjust 'start' if exact match is required
  all.x = TRUE
)

# Save the result
write.csv(gene_vcf, "merged_variants_with_genes.csv", row.names = FALSE)

matched_data <- read.table("matched_variants.bed", header = FALSE)
# Inspect the matched data
head(matched_data)


#load the geographical data
geo <- read_delim("accesions.txt")
geo= geo[,c("V1","V2","V3","V4","V6","V7","V8","V10","V11")]
colnames(geo)= c("accession","Sequenced by","Name","Country","latitude","longitude"
                 ,"Collector","Admixture Group","Sequenced By")

vcf_samples <- read.table("vcf_samples_with_positions.tsv", header = FALSE)
merged_metadata <- merge(vcf_samples, metadata, by = "accession", all.x = TRUE)


#####Rainfall data####
library(geodata)
library(raster)

# Specify a directory to save the data
# Ensure this directory exists or create it
data_path <- "worldclim_data"
if (!dir.exists(data_path)) dir.create(data_path)
# Load WorldClim data
clim <- worldclim_global(var = "prec", res = 10, path = "worldclim_data")


# Extract rainfall for each coordinate
geo$rainfall <- raster::extract(clim, geo[, c("latitude","longitude")])


