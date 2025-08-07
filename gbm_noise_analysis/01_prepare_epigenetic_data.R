# ==============================================================================
# SCRIPT: 01_prepare_epigenetic_data.R
#
# PURPOSE:
#   This script serves as the foundational first step of the analysis pipeline.
#   It gathers and prepares all necessary epigenetic data, creating a unified
#   master annotation file. This involves:
#     1. Reading methylation classifications from Cahn et al. (2024) and
#        Bewick et al. (2016).
#     2. Standardizing gene identifiers and classification terminology.
#     3. Integrating H2A.Z enrichment status for each gene.
#
#   The output is a single, clean CSV file that associates each gene with its
#   relevant epigenetic annotations, which will be used in all subsequent
#   downstream analyses.
# ==============================================================================

# --- 1. LOAD AND PROCESS CAHN ET AL. (2024) DATA ---

# Read Cahn data, which provides methylation status (gbM, TE-like, Unmethylated)
df_cahn <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/Cahn_et_al_2024.csv", stringsAsFactors = FALSE)
# Assign clear column names for merging
colnames(df_cahn) <- c("Gene_ID", "Cahn_Methylation_status")


# --- 2. LOAD AND PROCESS BEWICK ET AL. (2016) DATA ---

# Read Bewick data, which provides an independent methylation classification
df_bewick <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/Bewick_et_al_2016.csv", stringsAsFactors = FALSE)
# Assign clear column names for the first two columns
colnames(df_bewick)[1:2] <- c("Gene", "Bewick_Classification")

# Standardize gene IDs by removing isoform suffixes (e.g., ".1", ".2") to match other datasets
df_bewick$Gene_ID <- sub("\\..*", "", df_bewick$Gene)

# Standardize terminology: Ensure "Unmethylated" is consistent across datasets
df_bewick$Bewick_Classification[df_bewick$Bewick_Classification == "UM"] <- "Unmethylated"

# Keep only the essential columns for the final table
df_bewick <- df_bewick[, c("Gene_ID", "Bewick_Classification")]


# --- 3. LOAD H2A.Z ENRICHMENT DATA ---

# Read the list of genes identified as being depleted of H2A.Z in their gene body
dep <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/H2A.Z Body-Depleted Genes journal.csv", stringsAsFactors = FALSE)[,1]
# Read the list of genes identified as being enriched with H2A.Z in their gene body
enr <- read.csv("/group/sms029/mnieuwenh/Methylation_Data/H2A.Z Body-Enriched Genes journal.csv", stringsAsFactors = FALSE)[,1]


# --- 4. MERGE ALL DATASETS INTO A SINGLE TABLE ---

# Perform a full outer join to merge the two methylation datasets, keeping all genes from both
merged <- merge(df_cahn, df_bewick, by = "Gene_ID", all = TRUE)


# --- 5. ADD AGREEMENT AND H2A.Z STATUS COLUMNS ---

# Create a boolean column to easily identify genes classified as gbM by both studies.
# This is set to NA if data from either study is missing for a given gene.
merged$Both_gbM <- ifelse(is.na(merged$Cahn_Methylation_status) | is.na(merged$Bewick_Classification), NA,
                          merged$Cahn_Methylation_status == "gbM" & merged$Bewick_Classification == "gbM")

# Create a boolean column to identify genes with TE-like methylation patterns.
# Note that Bewick et al. splits this into mCHG and mCHH, so we check for either.
merged$Both_TE_like <- ifelse(is.na(merged$Cahn_Methylation_status) | is.na(merged$Bewick_Classification), NA,
                              merged$Cahn_Methylation_status == "TE-like methylation" &
                              (merged$Bewick_Classification == "mCHG" | merged$Bewick_Classification == "mCHH"))

# Create boolean columns for H2A.Z status. A gene is marked TRUE if in one list,
# FALSE if in the other, and NA if it is in neither list.
merged$H2AZ_Depleted <- ifelse(merged$Gene_ID %in% dep, TRUE,
                              ifelse(merged$Gene_ID %in% enr, FALSE, NA))
merged$H2AZ_Enriched <- ifelse(merged$Gene_ID %in% enr, TRUE,
                              ifelse(merged$Gene_ID %in% dep, FALSE, NA))


# --- 6. SAVE THE FINAL OUTPUT ---

# Define the output directory and create it if it doesn't already exist
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_prepare_epigenetic_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define the final output filename
output_file <- file.path(output_dir, "01_epigenetic_annotations.csv")

# Write the final, merged table to a CSV file
write.csv(merged, output_file, row.names = FALSE)

message(paste("Script finished. Master epigenetic annotation file saved to:", output_file))