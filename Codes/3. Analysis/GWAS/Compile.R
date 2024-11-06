####Step1####

#Identify values for plotting heatmap in SCENIC pythoin
#Save as CSV
#Python

#def export_heatmap_data(scplus_mudata, size_modality, color_modality, group_variable, eRegulon_metadata_key, size_feature_key, color_feature_key, feature_name_key, subset_feature_names=None):
#  size_matrix = scplus_mudata[size_modality].to_df()
#  color_matrix = scplus_mudata[color_modality].to_df()
#  group_by = scplus_mudata.obs[group_variable].tolist()
  
#  if subset_feature_names is None:
#    size_features, color_features, feature_names = scplus_mudata.uns[eRegulon_metadata_key][[size_feature_key, color_feature_key, feature_name_key]].drop_duplicates().values.T
#  else:
#    size_features, color_features, feature_names = scplus_mudata.uns[eRegulon_metadata_key][[size_feature_key, color_feature_key, feature_name_key]].drop_duplicates().query(f"{feature_name_key} in @subset_feature_names").values.T
  
#  plotting_df = generate_dotplot_df(
#    size_matrix=size_matrix,
#    color_matrix=color_matrix,
#    group_by=group_by,
#    size_features=size_features,
#    color_features=color_features,
#    feature_names=feature_names,
#    scale_size_matrix=True,    
#    scale_color_matrix=True,
#    group_name=group_variable,
#    size_name=size_modality,
#    color_name=color_modality,
#    feature_name=feature_name_key)
  
#  return plotting_df
#  dataframe_for_plot = export_heatmap_data(
#    scplus_mudata=scplus_mdata,
#    size_modality="extended_region_based_AUC",
#    color_modality="extended_gene_based_AUC",
#    group_variable="scRNA_counts:Clusters",
#    eRegulon_metadata_key="extended_e_regulon_metadata",
#    size_feature_key="Region_signature_name",
#    color_feature_key="Gene_signature_name",
#    feature_name_key="eRegulon_name"
#  )
  
#  dataframe_for_plot.to_csv("heatmap_plot_data.csv", index=False)

####Step2####
#From GWAS - extract traits and key TFs (human -> Mouse annotation)
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
data <- fread("Figures/GWAS/gwas_catalog_v1.0-associations_e112_r2024-07-27.tsv", header = TRUE, sep = "\t", quote = "")

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

####Step3####
SCENIC <-  read.csv("Figures/GWAS/PATTERN_heatmap_plot_data.csv")
# Filter and prepare the data
filtered_data <- SCENIC %>%
  # Filter based on eRegulon_name patterns for activators only
  filter(grepl("extended_\\+\\/\\+", eRegulon_name)) %>%
  select(scRNA_counts.Cluster_Pass2, eRegulon_name, extended_gene_based_AUC, extended_region_based_AUC) %>%
  mutate(
    # Correctly extract the TF name from eRegulon_name by removing the pattern and additional characters
    TF = toupper(gsub("_extended_[+/]+", "", eRegulon_name)),  # Simplified the regex by removing unnecessary escapes
    # Since only activators are selected, we directly assign the type
    type = "Activator"
  )

# View the first few rows of the filtered data to check your results
head(filtered_data)

####Matching####
# Reshape filtered_genes from wide to long format
filtered_genes_long <- filtered_genes %>%
  pivot_longer(
    cols = starts_with("GENE_"),
    names_to = "GENE_index",
    values_to = "Gene"
  ) %>%
  filter(!is.na(Gene))  # Remove NA values to reduce unnecessary joins

# Join with filtered_data based on the TF and Gene columns
augmented_data <- filtered_data %>%
  left_join(filtered_genes_long, by = c("TF" = "Gene"), relationship = "many-to-many") %>%
  select(-GENE_index)  # Removing the gene index column

####PLotting####
library(dplyr)
library(tidyr)
library(pheatmap)  # or library(ComplexHeatmap) for more complex heatmaps

# Calculate a weighted average AUC
weighted_avg_data <- augmented_data %>%
  group_by(scRNA_counts.Cluster_Pass2, DISEASE_TRAIT) %>%
  summarise(Weighted_Avg_AUC = sum(extended_gene_based_AUC * extended_gene_based_AUC) / sum(extended_gene_based_AUC), .groups = 'drop')

# Check the resulting summary
print(head(weighted_avg_data))

# First, aggregate the data by taking the maximum AUC
max_data <- weighted_avg_data %>%
  group_by(scRNA_counts.Cluster_Pass2, DISEASE_TRAIT) %>%
  summarise(Max_AUC = max(Weighted_Avg_AUC), .groups = 'drop')

# First, aggregate the data by taking the maximum AUC
#max_data <- augmented_data %>%
#  group_by(scRNA_counts.Cluster_Pass2, DISEASE_TRAIT) %>%
#  summarise(Max_AUC = max(extended_gene_based_AUC), .groups = 'drop')

# Create a wider format for the heatmap
heatmap_data <- max_data %>%
  pivot_wider(names_from = DISEASE_TRAIT, values_from = Max_AUC)

# Replace NA with 0 or appropriate minimal value if needed
heatmap_data[is.na(heatmap_data)] <- 0
#rownames(heatmap_data) <- heatmap_data$scRNA_counts.Cluster_Pass2
#heatmap_matrix <- as.matrix(heatmap_data[,-1])  # Remove the cluster names column for the heatmap

# Define color palette
#color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

##select
columns_of_interest <- c("Adult body size", "Adult onset asthma and/or BMI", "Age at first sexual intercourse", 
                         "Age of smoking initiation",     
                         "Age-related diseases and mortality", 
                         "Alcohol consumption (drinks per week)", 
                         "Alcohol consumption x hours spent using computers interaction", 
                         "Alcohol use disorder (consumption score)", 
                         "Alzheimer disease and age of onset", "Alzheimer's disease", 
                         "Alzheimer's disease (late onset)", 
                         "Anorexia nervosa (excluding migration to or from binge-eating disorder or bulimia nervosa)", 
                         "Anterior thigh muscle fat infiltration percentage", "Anti-saccade response", 
                         "Anti-thyroid drug induced agranulocytosis", "Antidepressant treatment resistance (> 2 drugs prescribed)", 
                         "Appendicular lean mass", "Aspartate aminotransferase levels", 
                         "Asthma", "Asthma (adult onset)", "Asthma (childhood onset)", 
                         "Biological sex", "Bipolar disorder",  "Body fat percentage", 
                         "Body mass index", 
                         "Body surface area", "Bone mineral density mean", "Bone mineral density variability", 
                         "Bone stiffness index", 
                         "Breast cancer", 
                         "Childhood aggressive behavior",
                         "Childhood body mass index", 
                         "Chronic kidney disease", "COVID-19 (covid vs negative)", 
                         "COVID-19 (hospitalized covid vs population)", 
                         "Depression", 
                         "Educational attainment",
                         "Emotional lability", "Emotional recognition","Erectile dysfunction", 
                         "Fasting blood glucose", "Fasting glucose", "Fasting insulin", 
                         "Fasting plasma glucose", 
                         "Glomerular filtration rate (creatinine)", "GLP-1 levels in response to oral glucose tolerance test (fasting)", 
                         "Glucose levels", 
                         "Headache", 
                         "Height",  "Hypertension", "Insomnia", "Intelligence", 
                         "Lewy body disease", "Lifetime smoking",  "Migraine", 
                         "Obesity-related traits",
                         "Offspring birth weight",
                         "Pericardial fat", "Preterm delivery",
                         "Prostate cancer",
                         "Schizophrenia", "Seasonality and depression", 
                         "Sleep disturbance (decreased sleep) in depressive disorders", 
                         "Smoking initiation",
                         "Systolic blood pressure",
                         "Testosterone levels", "Testosterone levels in postmenopausal women", 
                         "Thyroid stimulating hormone levels",
                         "Total cholesterol levels", 
                         "Triglyceride levels",
                         "Type 1 diabetes", "Type 2 diabetes", "Waist-hip index", 
                         "Waist-hip ratio")


# If you want to dynamically select columns containing 'diabetes'
#additional_columns <- grep("diabetes", colnames(heatmap_data), value = TRUE)
# Combine all specific columns you want to include
#all_columns <- c(columns_of_interest, additional_columns)

all_columns <- columns_of_interest

# Ensure all columns exist in the dataset
existing_columns <- intersect(all_columns, colnames(heatmap_data))

# Subset the data using the corrected column names
subset_heatmap_data <- heatmap_data %>%
  select(c("scRNA_counts.Cluster_Pass2", all_of(existing_columns)))

rownames(subset_heatmap_data) <- subset_heatmap_data$scRNA_counts.Cluster_Pass2

#if selection is needed
#selected_groups <-c("AntIDID",  "MMN", 
#                    "NPCGliogenic", "PMN", "PreThal", "PVNSON", 
#                    "SMN", "SMNProg", "SSTDMH", "TubARCVMH", "TubPMNProg")
#subset_heatmap_data <- heatmap_data %>%
#  filter(scRNA_counts.Cluster_Pass2 %in% selected_groups) %>%
#  select(c("scRNA_counts.Cluster_Pass2", all_of(existing_columns)))
#rownames(subset_heatmap_data) <- subset_heatmap_data$scRNA_counts.Cluster_Pass2
# Remove the cluster column before converting to matrix
# This assumes that 'scRNA_counts.Cluster_Pass2' is the first column
heatmap_matrix <- as.matrix(subset_heatmap_data[,-1])

# Verify that row names are still present
print(head(rownames(heatmap_matrix)))  # Should print the first few cluster names

# If the row names are not showing up correctly, let's diagnose:
if(is.null(rownames(heatmap_matrix))) {
  print("Row names are null. Setting them again.")
  rownames(heatmap_matrix) <- subset_heatmap_data$scRNA_counts.Cluster_Pass2
}

# Define color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# Draw the heatmap
pheatmap(heatmap_matrix,
         color = color_palette,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         border_color = NA  # No borders around cells
)

heatmap_matrixT <- t(heatmap_matrix)
# Draw the heatmap
pheatmap(heatmap_matrixT,
         color = color_palette,
         show_rownames = TRUE,
         show_colnames = TRUE,  # Ensure column names are shown, adjust if there are too many
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row",  # Scale rows to have zero mean and unit variance
         border_color = NA,  # No borders around cells
         fontsize = 8, fontsize_col = 8, cellwidth = 10,
         cellheight = 8, angle_col = 45,
         annotation_legend = TRUE
         )


