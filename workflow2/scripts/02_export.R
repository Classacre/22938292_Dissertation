#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(Seurat); library(Matrix); library(readr); library(dplyr)
})

opt <- list(
  make_option("--seurat", type="character"),
  make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$seurat)
DefaultAssay(obj) <- "RNA"

# counts recommended for noise modeling
counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(GetAssayData(obj, "counts"), "dgCMatrix")

meta <- obj@meta.data
meta <- meta %>%
  mutate(
    identity = if ("identity" %in% colnames(meta)) as.character(identity) else as.character(Idents(obj)),
    ers_state = case_when(
      grepl("ERS\\+$", identity) ~ "ERS+",
      grepl("ERS-$", identity) ~ "ERS-",
      TRUE ~ "NA"
    ),
    base_cell_type = trimws(gsub("ERS\\+?/?-?$", "", identity))
  )

# genes.tsv
genes <- rownames(counts)
gz1 <- gzfile(file.path(opt$outdir, "genes.tsv.gz"), "wb")
writeLines(genes, gz1)
close(gz1)

# mtx
Matrix::writeMM(counts, file.path(opt$outdir, "expression_matrix.mtx"))
system(paste("gzip -f", shQuote(file.path(opt$outdir, "expression_matrix.mtx"))))

# cell metadata
write.csv(meta, file.path(opt$outdir, "cell_metadata.csv"), row.names = TRUE)