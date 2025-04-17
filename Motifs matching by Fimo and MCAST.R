library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Load FIMO Data
fimo_df <- read_tsv("~/Desktop/fimo_out/fimo.tsv", col_types = cols())

# Select relevant columns
fimo_df <- fimo_df %>%
  select(sequence_name, motif_id, start, stop, strand, matched_sequence)

# Load MCAST Data
mcast_df <- read_tsv("~/Downloads/mcast (10).tsv", col_types = cols())

# Select relevant columns
mcast_df <- mcast_df %>%
  select(sequence_name, cluster_start = start, cluster_end = stop, cluster_score = score)

# Filter out sequences in FIMO that do not appear in MCAST
filtered_fimo <- fimo_df %>%
  semi_join(mcast_df, by = "sequence_name")

# Many-to-many join
merged_data <- filtered_fimo %>%
  left_join(mcast_df, by = "sequence_name", relationship = "many-to-many") %>%
  #filter(start >= cluster_start & stop <= cluster_end) %>%
  arrange(sequence_name, start)

# Ensure motifs do not start before the previous one stops
motif_locations <- merged_data %>%
  group_by(sequence_name) %>%
  mutate(previous_stop = lag(stop, default = first(start))) %>%
  filter(start >= previous_stop) %>%
  ungroup()

# 🚀 Check if there are too many sequences
num_sequences <- n_distinct(motif_locations$sequence_name)
print(paste("Total unique sequences in plot:", num_sequences))

# Optionally filter the top 50 most common sequences to avoid overcrowding
top_sequences <- motif_locations %>%
  count(sequence_name, sort = TRUE) %>%
  slice_head(n = 50) %>%
  pull(sequence_name)

filtered_plot_data <- motif_locations %>%
  filter(sequence_name %in% top_sequences)

# 🎨 Create the Scatter Plot
ggplot(filtered_plot_data, aes(x = start, y = sequence_name, color = motif_id)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Motif Positions Along Sequences",
       x = "Motif Start Position",
       y = "Sequence Name",
       color = "Motif ID") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))  # Adjust text size for readability



#######Reshape#####
library(dplyr)
library(tidyr)
library(readr)

# Count occurrences and store locations
motif_summary <- motif_locations %>%
  group_by(sequence_name, motif_id) %>%
  summarise(
    motif_count = n(),
    motif_location = paste(start,"-",stop, collapse = ","),
    .groups = "drop"
  )

# Pivot wider so each motif gets its own count and location columns
pivot_motif_df <- motif_summary %>%
  pivot_wider(
    names_from = motif_id,
    values_from = c(motif_count, motif_location),
    names_glue = "{motif_id}_{.value}",  # Correct naming convention
    values_fill = list(motif_count = 0, motif_location = "NA")
  )

# View results
print(pivot_motif_df)

##Rename the columns
colnames(pivot_motif_df)= gsub("_motif","",colnames(pivot_motif_df))

# Save to CSV
write_csv(pivot_motif_df, "17_accession_motifs.csv")

