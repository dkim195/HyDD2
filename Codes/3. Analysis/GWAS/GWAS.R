library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

# Set working directory and load the data
set.seed(1234)
setwd("/media/thomaskim/Data/")

# Load the data
#https://www.ebi.ac.uk/gwas/docs/file-downloads -> All associations v1.0 
data <- fread("GWAS/gwas_catalog_v1.0-associations_e112_r2024-07-08.tsv", header = TRUE, sep = "\t", quote = "")

# Extract the required columns
extracted_data <- data[, .(`DISEASE/TRAIT`, `MAPPED_GENE`)]

# Rename columns to remove special characters for easier handling
setnames(extracted_data, old = c("DISEASE/TRAIT", "MAPPED_GENE"), new = c("DISEASE_TRAIT", "GENES"))

# Determine the maximum number of genes
max_genes <- max(sapply(strsplit(extracted_data$GENES, " - "), length))

# Split the 'MAPPED_GENES' column into multiple columns
split_genes <- separate(extracted_data, GENES, into = paste0("GENE_", 1:max_genes), sep = " - ", fill = "right")

# Download the list of transcription factors
tf_list_url <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
tf_list <- fread(tf_list_url, header = FALSE, sep = "\t")
transcription_factors <- tf_list$V1

# Filter rows that contain transcription factors in any of the GENE_ columns
is_tf <- apply(split_genes[, paste0("GENE_", 1:max_genes)], 1, function(row) {
  any(row %in% transcription_factors, na.rm = TRUE)
})

filtered_genes <- split_genes[is_tf, ]

# Save the extracted data to a new TSV file
extracted_file_path <- "Figures/GWAS/GWAS_extracted_data.tsv"  # Replace with the desired save path
fwrite(filtered_genes, extracted_file_path, sep="\t")

####eRegulon####

SCENIC <- read.table("SCENIC/eRegulon/PATTERN_eRegulons_extended.tsv", header = TRUE, sep = "\t")

# Filter and prepare the data
filtered_data <- SCENIC %>%
  filter(grepl("extended_\\+/\\+", eRegulon_name) | grepl("extended_\\-/\\+", eRegulon_name)) %>%
  select(TF, eRegulon_name, Gene, importance_TF2G, triplet_rank) %>%
  mutate(
    TF = toupper(TF),
    type = ifelse(grepl("\\+/\\+", eRegulon_name), "Activator", ifelse(grepl("\\-/\\+", eRegulon_name), "Repressor", NA))
  )

####Matching####
tf_names <- unique(filtered_data$TF)
matched_rows <- apply(filtered_genes[, paste0("GENE_", 1:max_genes)], 1, function(row) {
  any(row %in% tf_names, na.rm = TRUE)
})
matched_genes <- filtered_genes[matched_rows, ]
matched_file_path <- "Figures/GWAS/matched_disease_trait_split_genes.tsv"
fwrite(matched_genes, matched_file_path, sep = "\t")

####Summarizing####
long_genes <- matched_genes %>%
  gather(key = "gene_num", value = "Gene", starts_with("GENE_")) %>%
  filter(!is.na(Gene))

# Group by DISEASE_TRAIT and summarize the unique genes
combined_genes <- long_genes %>%
  group_by(DISEASE_TRAIT) %>%
  summarize(Associated_Genes = paste(unique(Gene), collapse = ", "))

combined_file_path <- "Figures/GWAS/combined_disease_trait_genes.tsv"
fwrite(combined_genes, combined_file_path, sep = "\t")