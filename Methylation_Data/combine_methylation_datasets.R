# Script: combine_methylation_datasets.R
# Purpose: Combine Cahn_et_al_2024 and Bewick_et_al_2016 methylation data, add agreement columns, H2A.Z enrichment columns, and output CSV file (base R only).

# Read Cahn data
df_cahn <- read.csv("/group/sms029/mnieuwenh/gbM_data/Cahn_et_al_2024.csv", stringsAsFactors = FALSE)
colnames(df_cahn) <- c("Gene_ID", "Cahn_Methylation_status")

# Read Bewick data
df_bewick <- read.csv("/group/sms029/mnieuwenh/gbM_data/Bewick_et_al_2016.csv", stringsAsFactors = FALSE)
colnames(df_bewick)[1:2] <- c("Gene", "Bewick_Classification")

# Standardize gene IDs (remove .x isoform from Bewick)
df_bewick$Gene_ID <- sub("\\..*", "", df_bewick$Gene)

# Make Bewick naming consistent: UM -> Unmethylated
df_bewick$Bewick_Classification[df_bewick$Bewick_Classification == "UM"] <- "Unmethylated"

# Only keep relevant columns from Bewick
df_bewick <- df_bewick[, c("Gene_ID", "Bewick_Classification")]

# Read H2A.Z depleted and enriched gene lists
dep <- read.csv("/group/sms029/mnieuwenh/gbM_data/H2A.Z Body-Depleted Genes journal.csv", stringsAsFactors = FALSE)[,1]
enr <- read.csv("/group/sms029/mnieuwenh/gbM_data/H2A.Z Body-Enriched Genes journal.csv", stringsAsFactors = FALSE)[,1]

# Merge datasets
merged <- merge(df_cahn, df_bewick, by = "Gene_ID", all = TRUE)

# Agreement columns (NA if missing data)
merged$Both_gbM <- ifelse(is.na(merged$Cahn_Methylation_status) | is.na(merged$Bewick_Classification), NA,
                          merged$Cahn_Methylation_status == "gbM" & merged$Bewick_Classification == "gbM")
merged$Both_TE_like <- ifelse(is.na(merged$Cahn_Methylation_status) | is.na(merged$Bewick_Classification), NA,
                              merged$Cahn_Methylation_status == "TE-like methylation" & 
                              (merged$Bewick_Classification == "mCHG" | merged$Bewick_Classification == "mCHH"))

# H2A.Z columns (NA if not present in either list)
merged$H2AZ_Depleted <- ifelse(merged$Gene_ID %in% dep, TRUE,
                              ifelse(merged$Gene_ID %in% enr, FALSE, NA))
merged$H2AZ_Enriched <- ifelse(merged$Gene_ID %in% enr, TRUE,
                              ifelse(merged$Gene_ID %in% dep, FALSE, NA))

# Output CSV
write.csv(merged, "/group/sms029/mnieuwenh/gbM_data/combined_methylation_data.csv", row.names = FALSE)
