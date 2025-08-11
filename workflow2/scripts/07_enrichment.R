#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(clusterProfiler); library(org.At.tair.db); library(enrichplot); library(GOSemSim); library(ggplot2)
})

opt <- list(
  make_option("--csv", type="character"),
  make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(opt$csv) %>% filter(is.finite(scran_resid_var), !is.na(gene))
# Rank by centered average scran_resid_var
rank_global <- df %>%
  group_by(gene) %>%
  summarise(avg_noise = mean(scran_resid_var, na.rm=TRUE), .groups="drop") %>%
  filter(!is.na(avg_noise))
z <- as.numeric(scale(rank_global$avg_noise, center=TRUE, scale=TRUE))
names(z) <- rank_global$gene
z <- sort(z, decreasing = TRUE)

set.seed(123)
gsea <- tryCatch(gseGO(geneList = z, OrgDb = org.At.tair.db, ont="BP", keyType="TAIR",
                       minGSSize=10, maxGSSize=5000, pvalueCutoff=0.05, pAdjustMethod="BH", eps=0, verbose=FALSE),
                 error=function(e) NULL)

if (!is.null(gsea) && nrow(gsea@result) > 0) {
  p <- tryCatch(dotplot(gsea, showCategory=20, split=".sign") + facet_grid(~.sign) +
                  labs(title="GSEA on corrected noise (global)"), error=function(e) NULL)
  if (!is.null(p)) ggsave(file.path(opt$outdir, "gsea_global.png"), p, width=12, height=10, dpi=300)
}

sink(file.path(opt$outdir, "gsea_summary.txt"))
cat("GSEA summary (clusterProfiler::gseGO, BP)\n")
cat("Date:", as.character(Sys.Date()), "\n")
if (!is.null(gsea)) {
  cat("Significant terms:", sum(gsea@result$p.adjust < 0.05, na.rm=TRUE), "\n")
} else {
  cat("gseGO failed or returned no results.\n")
}
sink()