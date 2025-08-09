# ==============================================================================
# SCRIPT: 01_prepare_epigenetic_data.R
#
# PURPOSE:
#   This script gathers and prepares all necessary epigenetic and promoter
#   architecture data, creating a unified master annotation file. It reads
#   methylation classifications, H2A.Z enrichment status, and includes an
#   analysis of TATA box presence in promoters.
# ==============================================================================

# --- 1. SETUP and DATA LOADING ---
cahn_path <- "/group/sms029/mnieuwenh/Methylation_Data/Cahn_et_al_2024.csv"
bewick_path <- "/group/sms029/mnieuwenh/Methylation_Data/Bewick_et_al_2016.csv"
h2az_dep_path <- "/group/sms029/mnieuwenh/Methylation_Data/H2A.Z Body-Depleted Genes journal.csv"
h2az_enr_path <- "/group/sms029/mnieuwenh/Methylation_Data/H2A.Z Body-Enriched Genes journal.csv"
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/01_prepare_epigenetic_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "01_epigenetic_annotations.csv")
summary_file_path <- file.path(output_dir, "01_epigenetic_data_summary.txt")


# --- 2. LOAD AND PROCESS CAHN ET AL. (2024) DATA ---
df_cahn <- read.csv(cahn_path, stringsAsFactors = FALSE)
colnames(df_cahn) <- c("Gene_ID", "Cahn_Methylation_status")

# --- 3. LOAD AND PROCESS BEWICK ET AL. (2016) DATA ---
df_bewick <- read.csv(bewick_path, stringsAsFactors = FALSE)
colnames(df_bewick)[1:2] <- c("Gene", "Bewick_Classification")
df_bewick$Gene_ID <- sub("\\..*", "", df_bewick$Gene) # Remove isoform suffixes
df_bewick$Bewick_Classification[df_bewick$Bewick_Classification == "UM"] <- "Unmethylated"
df_bewick <- df_bewick[, c("Gene_ID", "Bewick_Classification")]

# --- 4. LOAD H2A.Z ENRICHMENT DATA ---
dep <- read.csv(h2az_dep_path, stringsAsFactors = FALSE)[,1]
enr <- read.csv(h2az_enr_path, stringsAsFactors = FALSE)[,1]


# --- 4.5. IDENTIFY TATA-CONTAINING GENES ---
message("Identifying TATA-box containing promoters...")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org/")
packages_to_install <- c("BSgenome.Athaliana.TAIR.TAIR9", "TxDb.Athaliana.BioMart.plantsmart28", "Biostrings", "dplyr", "GenomeInfoDb")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

genome <- BSgenome.Athaliana.TAIR.TAIR9
txdb <- TxDb.Athaliana.BioMart.plantsmart28

# **ROBUST CORRECTION 1: Manually rename seqlevels for a guaranteed match.**
# We now include the chloroplast/plastid mapping ('Pt' -> 'ChrC').
message("Original seqlevels from TxDb: ", paste(seqlevels(txdb), collapse=", "))
message("Target seqlevels from BSgenome: ", paste(seqlevels(genome), collapse=", "))
mapping <- c("1"="Chr1", "2"="Chr2", "3"="Chr3", "4"="Chr4", "5"="Chr5", "Mt"="ChrM", "Pt"="ChrC")
txdb <- renameSeqlevels(txdb, mapping)
message("Corrected seqlevels in TxDb: ", paste(seqlevels(txdb), collapse=", "))


promoter_regions <- promoters(txdb, upstream = 100, downstream = 50, columns = c("tx_name", "gene_id"))

# **ROBUST CORRECTION 2: Trim out-of-bounds ranges IMMEDIATELY after creation.**
# This is the key fix for the "out of bounds" error.
message("Trimming promoter ranges to fit within chromosome boundaries...")
promoter_regions <- trim(promoter_regions)

promoter_regions$gene_id <- unlist(lapply(promoter_regions$gene_id, function(x) if(is.null(x)) NA else x[1]))

# Filter for the chromosomes of interest after all corrections.
promoter_regions <- promoter_regions[seqnames(promoter_regions) %in% c("Chr1","Chr2","Chr3","Chr4","Chr5") & !is.na(promoter_regions$gene_id)]

promoter_regions_unique <- as.data.frame(promoter_regions) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# The getSeq call should now succeed without any errors or warnings.
promoter_seqs <- getSeq(genome, promoter_regions_unique)
tata_counts <- vcountPattern("TATAA", promoter_seqs, max.mismatch = 1)
df_tata <- data.frame(
  Gene_ID = promoter_regions_unique$gene_id,
  has_TATA_box = tata_counts > 0
)


# --- 5. MERGE ALL DATASETS AND CREATE ANNOTATION COLUMNS ---
message("Merging all epigenetic and promoter datasets...")
merged <- merge(df_cahn, df_bewick, by = "Gene_ID", all = TRUE)
merged <- merge(merged, df_tata, by = "Gene_ID", all.x = TRUE)
merged$has_TATA_box[is.na(merged$has_TATA_box)] <- FALSE
merged$H2AZ_Depleted <- merged$Gene_ID %in% dep
merged$H2AZ_Enriched <- merged$Gene_ID %in% enr

# --- 6. SAVE THE FINAL OUTPUT ---
write.csv(merged, output_file, row.names = FALSE)
message(paste("Script finished. Master epigenetic annotation file saved to:", output_file))

# --- 7. GENERATE STATISTICAL SUMMARY ---
message("Generating statistical summary file...")
summary_text <- c(
  "=================================================================",
  "         STATISTICAL SUMMARY: 01_prepare_epigenetic_data.R",
  "=================================================================",
  "", paste("Summary generated on:", Sys.Date()), "",
  "--- Input Data Dimensions ---",
  paste("Cahn et al. (2024) records:", nrow(df_cahn)),
  paste("Bewick et al. (2016) records:", nrow(df_bewick)),
  paste("H2A.Z Depleted gene list size:", length(dep)),
  paste("H2A.Z Enriched gene list size:", length(enr)),
  paste("Total unique promoters scanned for TATA box:", nrow(df_tata)), "",
  "--- Merged Data Summary ---",
  paste("Total unique genes in merged file:", nrow(merged)), "",
  "--- Gene Counts by Category (in final merged file) ---",
  "Cahn et al. Classifications:",
  paste("  - ", names(table(merged$Cahn_Methylation_status, useNA = "ifany")), ": ", table(merged$Cahn_Methylation_status, useNA = "ifany"), " genes", collapse = "\n"), "",
  "Bewick et al. Classifications:",
  paste("  - ", names(table(merged$Bewick_Classification, useNA = "ifany")), ": ", table(merged$Bewick_Classification, useNA = "ifany"), " genes", collapse = "\n"), "",
  "H2A.Z Status:",
  paste("  - H2A.Z Depleted:", sum(merged$H2AZ_Depleted, na.rm = TRUE)),
  paste("  - H2A.Z Enriched:", sum(merged$H2AZ_Enriched, na.rm = TRUE)),
  paste("  - Overlap (Enriched & Depleted):", sum(merged$H2AZ_Enriched & merged$H2AZ_Depleted, na.rm = TRUE)), "",
  "Promoter Architecture:",
  paste("  - Genes with a TATA box in promoter:", sum(merged$has_TATA_box, na.rm = TRUE)),
  paste("  - Genes without a TATA box in promoter:", sum(!merged$has_TATA_box, na.rm = TRUE)), ""
)
writeLines(summary_text, con = summary_file_path)
message(paste("Statistical summary saved to:", summary_file_path))
message("Script 01 finished successfully.")