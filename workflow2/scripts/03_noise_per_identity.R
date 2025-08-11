#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(Matrix); library(data.table)
  library(SingleCellExperiment); library(scuttle); library(scran)
  library(sctransform); library(dplyr); library(readr)
})

opt <- list(
  make_option("--mtx", type="character"),
  make_option("--genes", type="character"),
  make_option("--cells", type="character"),
  make_option("--anno", type="character"),
  make_option("--outdir", type="character"),
  make_option("--min_cells", type="integer", default=50),
  make_option("--max_cells", type="integer", default=500),
  make_option("--min_mean", type="double", default=0.01),
  make_option("--min_pct", type="double", default=0.10)
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

message("Reading matrix and metadata...")
mtx <- readMM(opt$mtx)
genes <- fread(cmd=paste("zcat", opt$genes), header=FALSE)$V1
cells <- read.csv(opt$cells, row.names = 1, check.names = FALSE)
stopifnot(nrow(mtx) == length(genes))
colnames(mtx) <- rownames(cells); rownames(mtx) <- genes

identity <- as.character(cells$identity)
stopifnot(length(identity) == ncol(mtx))

# Cap cells per identity
set.seed(1)
split_cells <- split(colnames(mtx), identity)
split_cells <- lapply(split_cells, function(v) if (length(v) > opt$max_cells) sample(v, opt$max_cells) else v)
identities <- names(split_cells)
keep <- vapply(split_cells, length, integer(1)) >= opt$min_cells
identities <- identities[keep]
split_cells <- split_cells[keep]

sce_all <- SingleCellExperiment(assays = list(counts = mtx))
rm(mtx); gc()

noise_list <- vector("list", length(split_cells))
names(noise_list) <- identities

message("Computing per-identity noise metrics...")
for (id in identities) {
  cols <- split_cells[[id]]
  sce <- sce_all[, cols, drop=FALSE]

  # Normalize and model variance
  sce <- scuttle::computeLibraryFactors(sce)
  sce <- scuttle::logNormCounts(sce)
  mv <- scran::modelGeneVar(sce)

  counts  <- assay(sce, "counts")
  mu <- Matrix::rowMeans(counts)
  vr <- apply(counts, 1, var)
  pct <- Matrix::rowMeans(counts > 0)

  resid_var <- rep(NA_real_, nrow(sce))
  if (ncol(sce) >= 100) {
    vst_out <- tryCatch(sctransform::vst(as.matrix(counts), return_gene_attr = TRUE, verbosity = 0), error=function(e) NULL)
    if (!is.null(vst_out)) {
      resid_var <- vst_out$gene_attr$residual_variance
      names(resid_var) <- rownames(sce)
    }
  }

  noise_list[[id]] <- data.frame(
    gene = rownames(sce),
    identity = id,
    n_cells = ncol(sce),
    mean_expr = mu / ncol(sce),
    variance = vr,
    fano = vr / pmax(mu, 1e-8),
    cv2 = (vr / pmax(mu, 1e-8)^2),
    scran_resid_var = mv$bio,
    scran_total_var = mv$mean,
    pct_expressing = as.numeric(pct),
    sct_resid_var = resid_var[rownames(sce)],
    stringsAsFactors = FALSE
  )
}

noise_dt <- dplyr::bind_rows(noise_list)
rm(noise_list, sce_all); gc()

anno <- read.csv(opt$anno, stringsAsFactors = FALSE) %>%
  dplyr::distinct(Gene_ID, .keep_all = TRUE) %>%
  dplyr::rename(gene = Gene_ID)

noise_dt <- noise_dt %>%
  dplyr::left_join(
    anno %>%
      dplyr::transmute(
        gene,
        cahn_group = dplyr::case_when(
          Cahn_Methylation_status == "gbM" ~ "gbM",
          Cahn_Methylation_status == "TE-like methylation" ~ "TE-like",
          Cahn_Methylation_status == "Unmethylated" ~ "Unmethylated",
          TRUE ~ NA_character_
        ),
        bewick_group = Bewick_Classification,
        h2az_group = dplyr::case_when(
          H2AZ_Depleted ~ "H2A.Z-Depleted",
          H2AZ_Enriched ~ "H2A.Z-Enriched",
          TRUE ~ "Other"
        ),
        has_TATA_box = has_TATA_box,
        gbm_union = gbm_union,
        gbm_intersection = gbm_intersection
      ), by = "gene"
  )

filtered <- noise_dt %>%
  dplyr::filter(is.finite(mean_expr), is.finite(pct_expressing),
                mean_expr > opt$min_mean, pct_expressing > opt$min_pct)

gz <- gzfile(file.path(opt$outdir, "noise_per_identity.csv.gz"), "wb")
write.csv(noise_dt, gz, row.names = FALSE)
close(gz)

sink(file.path(opt$outdir, "noise_summary.txt"))
cat("Per-identity noise summary\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Identities processed:", length(unique(noise_dt$identity)), "\n")
cat("Rows total:", nrow(noise_dt), "\n")
cat("Rows after base filter:", nrow(filtered), "\n")
cat("scran_resid_var finite:", sum(is.finite(noise_dt$scran_resid_var)), "\n")
cat("sct_resid_var finite:", sum(is.finite(noise_dt$sct_resid_var)), "\n")
sink()