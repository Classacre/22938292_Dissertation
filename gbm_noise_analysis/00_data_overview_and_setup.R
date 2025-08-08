# ==============================================================================
# SCRIPT: 00_data_overview_and_setup.R
#
# PURPOSE:
#   This script serves as a preliminary step to provide essential context for the
#   entire analysis pipeline. It addresses two key areas:
#
#   1. Epigenetic Data Concordance: It visualizes the agreement between the
#      two independent methylation datasets (Cahn et al., 2024 and
#      Bewick et al., 2016), specifically for the 'gbM' classification,
#      using a Venn diagram. This clarifies how consistently gbM is defined.
#
#   2. Single-Cell Data Overview: It provides a high-level summary and
#      visualization of the unpublished Seurat dataset. This includes:
#        - A UMAP plot to visualize cell type clusters.
#        - A bar chart showing the number of cells per cell type.
#        - A summary of data dimensions (cells, genes).
#
#   This script does not generate data for downstream processing but creates
#   key summary figures and statistics that contextualize the entire project.
#
# DISCLAIMER:
#   The single-cell RNA-seq (Seurat) data is from an unpublished study and is
#   used here with permission. It is not the work of the author of this
#   pipeline.
# ==============================================================================


# --- 1. SETUP: Load libraries and define paths ---

# Load packages required for data manipulation and plotting
# ggvenn is used for creating the Venn diagram
packages_to_load <- c("Seurat", "dplyr", "ggplot2", "patchwork", "ggvenn")
for (pkg in packages_to_load) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Install from CRAN; Seurat is also on CRAN now.
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# --- Define I/O Paths ---
# INPUTS
# Paths to the raw epigenetic data
cahn_path <- "/group/sms029/mnieuwenh/Methylation_Data/Cahn_et_al_2024.csv"
bewick_path <- "/group/sms029/mnieuwenh/Methylation_Data/Bewick_et_al_2016.csv"
# Path to the primary Seurat object
seurat_path <- "/group/sms029/Oliva_dataset/integrated_col_trajectories_colonly.rds"

# OUTPUTS
# A dedicated directory for this script's summary outputs
output_dir <- "/group/sms029/mnieuwenh/gbm_noise_analysis/results/00_data_overview"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define full paths for the output figures
venn_diagram_path <- file.path(output_dir, "00_gbm_classification_agreement.png")
seurat_overview_path <- file.path(output_dir, "00_seurat_data_overview.png")
summary_text_path <- file.path(output_dir, "00_data_summary.txt")


# --- 2. VISUALIZE AGREEMENT OF GBM CLASSIFICATIONS ---

message("Loading and processing epigenetic data for Venn diagram...")

# Load Cahn et al. (2024) data
df_cahn <- read.csv(cahn_path, stringsAsFactors = FALSE)
colnames(df_cahn) <- c("Gene_ID", "Cahn_Methylation_status")
cahn_gbm_genes <- df_cahn %>%
  filter(Cahn_Methylation_status == "gbM") %>%
  pull(Gene_ID)

# Load Bewick et al. (2016) data and standardize Gene IDs
df_bewick <- read.csv(bewick_path, stringsAsFactors = FALSE)
colnames(df_bewick)[1:2] <- c("Gene", "Bewick_Classification")
df_bewick$Gene_ID <- sub("\\..*", "", df_bewick$Gene)
bewick_gbm_genes <- df_bewick %>%
  filter(Bewick_Classification == "gbM") %>%
  pull(Gene_ID)

# Create a list of the two gene sets for the Venn diagram
gbm_gene_lists <- list(
  `Cahn et al. (2024)` = cahn_gbm_genes,
  `Bewick et al. (2016)` = bewick_gbm_genes
)

message("Generating Venn diagram for gbM gene overlap...")

# Generate the Venn diagram using ggvenn
p_venn <- ggvenn(
  gbm_gene_lists,
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 5
) +
  labs(title = "Agreement of gbM Gene Classification",
       subtitle = "Overlap of genes classified as 'gbM' in two independent studies") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# Save the plot
ggsave(venn_diagram_path, plot = p_venn, width = 8, height = 6, bg = "white")
message(paste("Venn diagram saved to:", venn_diagram_path))


# --- 3. SUMMARIZE AND VISUALIZE SEURAT DATASET ---

message("Loading Seurat object for data overview...")
# Check for file existence before loading
if (!file.exists(seurat_path)) {
  stop("Error: Seurat object file not found. Please check the path.")
}
seurat_obj <- readRDS(seurat_path)
message("Seurat object loaded successfully.")

# --- Panel A: UMAP Visualization ---
message("Generating UMAP plot...")

# ============================================================================ #
# === FIX APPLIED HERE ===
# The reduction name was changed from "umap" to "umap.rpca" based on the
# diagnostic script output.
# ============================================================================ #
p_umap <- DimPlot(seurat_obj, reduction = "umap.rpca", group.by = "identity", label = TRUE, repel = TRUE) +
  labs(title = "UMAP of Cell Types",
       subtitle = "Visualization of single-cell clusters (using 'umap.rpca')") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# --- Panel B: Cell Type Distribution ---
message("Generating cell count summary plot...")
cell_counts <- as.data.frame(table(seurat_obj$identity))
colnames(cell_counts) <- c("CellType", "Count")

p_counts <- ggplot(cell_counts, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Count), vjust = -0.3, size = 3.5) +
  labs(title = "Cell Counts per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

# --- Assemble the overview figure ---
message("Assembling final Seurat overview figure...")
p_seurat_overview <- p_umap + p_counts +
  plot_annotation(
    title = "Overview of the Single-Cell RNA-Seq Dataset",
    subtitle = "Unpublished data used for expression and noise analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

# Save the combined plot
ggsave(seurat_overview_path, plot = p_seurat_overview, width = 14, height = 7, bg = "white")
message(paste("Seurat overview figure saved to:", seurat_overview_path))


# --- 4. WRITE SUMMARY TEXT FILE ---

message("Writing summary statistics to text file...")

# Capture summary information
summary_text <- c(
  "=================================================================",
  "            DATASET OVERVIEW & CONTEXT",
  "=================================================================",
  "",
  paste("Summary generated on:", Sys.Date()),
  "",
  "--- Epigenetic Data Summary ---",
  paste("Number of gbM genes in Cahn et al.:", length(cahn_gbm_genes)),
  paste("Number of gbM genes in Bewick et al.:", length(bewick_gbm_genes)),
  paste("Number of gbM genes common to both:", length(intersect(cahn_gbm_genes, bewick_gbm_genes))),
  paste("A Venn diagram visualizing this overlap has been saved to:", venn_diagram_path),
  "",
  "--- Single-Cell (Seurat) Data Summary ---",
  "DISCLAIMER: This is unpublished data from an external source.",
  paste("Total number of cells:", ncol(seurat_obj)),
  paste("Total number of genes:", nrow(seurat_obj)),
  paste("Number of identified cell types:", length(unique(seurat_obj$identity))),
  "",
  "Cell counts per type:",
  paste("  - ", cell_counts$CellType, ": ", cell_counts$Count, " cells", sep = "", collapse = "\n"),
  "",
  paste("A figure summarizing the Seurat data has been saved to:", seurat_overview_path)
)

# Write to file
writeLines(summary_text, con = summary_text_path)
message(paste("Summary text file saved to:", summary_text_path))

message("Script 00 finished successfully. Contextual figures and summaries are generated.")