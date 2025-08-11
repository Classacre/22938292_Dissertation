#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(Seurat)
  library(eulerr)
})

opt <- list(
  make_option("--cahn", type="character"),
  make_option("--bewick", type="character"),
  make_option("--seurat", type="character"),
  make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cahn <- read.csv(opt$cahn, stringsAsFactors = FALSE)
colnames(cahn) <- c("Gene_ID", "Cahn_Methylation_status")
cahn_gbm <- cahn %>% filter(Cahn_Methylation_status == "gbM") %>% pull(Gene_ID)

bew <- read.csv(opt$bewick, stringsAsFactors = FALSE)
colnames(bew)[1:2] <- c("Gene", "Bewick_Classification")
bew$Gene_ID <- sub("\\..*", "", bew$Gene)
bew_gbm <- bew %>% filter(Bewick_Classification == "gbM") %>% pull(Gene_ID)

fit <- eulerr::euler(list(`Cahn et al.` = cahn_gbm, `Bewick et al.` = bew_gbm))
png(file.path(opt$outdir, "00_gbm_overlap.png"), width=1000, height=800, res=150)
plot(fit, fills = list(fill = c("#0073C2FF","#EFC000FF"), alpha=0.6),
     labels = list(font=2), quantities = TRUE, edges = list(lwd=1.5))
dev.off()

if (!file.exists(opt$seurat)) stop("Seurat RDS not found")
obj <- readRDS(opt$seurat)
idents <- if ("identity" %in% colnames(obj@meta.data)) as.character(obj$identity) else Idents(obj)

cell_counts <- sort(table(idents), decreasing = TRUE)
sum_lines <- c(
  "=================================================================",
  "DATASET OVERVIEW & CONTEXT",
  "=================================================================",
  paste("Summary generated on:", Sys.Date()),
  "",
  "--- Epigenetic Data Summary ---",
  paste("Cahn gbM:", length(unique(cahn_gbm))),
  paste("Bewick gbM:", length(unique(bew_gbm))),
  paste("Overlap (intersection):", length(intersect(cahn_gbm, bew_gbm))),
  paste("Venn saved:", file.path(opt$outdir, "00_gbm_overlap.png")),
  "",
  "--- Seurat Object ---",
  paste("Cells:", ncol(obj)),
  paste("Genes:", nrow(obj)),
  paste("Identities:", length(unique(idents))),
  "Cell counts per identity:"
)
sum_lines <- c(sum_lines, paste(sprintf("  - %s: %d", names(cell_counts), as.integer(cell_counts)), collapse = "\n"))

writeLines(sum_lines, file.path(opt$outdir, "00_data_summary.txt"))