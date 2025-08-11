#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(readr)
})

opt <- list(
  make_option("--seurat", type="character"),
  make_option("--noise", type="character"),
  make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$seurat)
meta <- obj@meta.data %>%
  mutate(
    identity_char = if ("identity" %in% colnames(.)) as.character(identity) else as.character(Idents(obj)),
    ers_state = case_when(grepl("ERS\\+$", identity_char) ~ "ERS+",
                          grepl("ERS-$", identity_char) ~ "ERS-",
                          TRUE ~ "NA"),
    base_cell_type = trimws(gsub("ERS\\+?/?-?$", "", identity_char))
  )
obj@meta.data <- meta

noise <- read.csv(opt$noise)
noise <- noise %>% filter(is.finite(scran_resid_var))

# Per-cell-type DE (ERS+ vs ERS-), then global responsive set = union of significant in any type
Idents(obj) <- obj$base_cell_type
cts <- sort(unique(obj$base_cell_type))
LOGFC <- 0.25; PADJ <- 0.05

responsive_sets <- list()
for (ct in cts) {
  cells_ct <- rownames(obj@meta.data)[obj$base_cell_type == ct & obj$ers_state %in% c("ERS+","ERS-")]
  if (length(cells_ct) < 100) next
  sub <- subset(obj, cells = cells_ct)
  if (length(unique(sub$ers_state)) < 2) next
  Idents(sub) <- sub$ers_state
  de <- tryCatch(FindMarkers(sub, ident.1 = "ERS+", ident.2 = "ERS-", logfc.threshold = 0,
                             min.pct = 0.1, test.use = "wilcox"), error=function(e) NULL)
  if (is.null(de)) next
  de$gene <- rownames(de)
  resp <- de %>% filter(p_val_adj < PADJ, avg_log2FC > LOGFC) %>% pull(gene)
  responsive_sets[[ct]] <- resp
}
responsive_genes <- unique(unlist(responsive_sets))

noise2 <- noise %>%
  mutate(is_responsive = gene %in% responsive_genes) %>%
  filter(mean_expr > 0.01, pct_expressing > 0.10)

# Global comparison
w <- tryCatch(wilcox.test(scran_resid_var ~ is_responsive, data = noise2), error=function(e) NULL)

pA <- ggplot(noise2, aes(x = is_responsive, y = scran_resid_var, fill = is_responsive)) +
  geom_violin(trim=FALSE, alpha=0.6) +
  geom_boxplot(width=0.12, outlier.shape = NA, fill="white") +
  scale_fill_manual(values = c("FALSE"="cornflowerblue","TRUE"="coral1")) +
  scale_x_discrete(labels=c("FALSE"="Non-responsive","TRUE"="Responsive")) +
  labs(title="Corrected noise vs responsiveness", subtitle = ifelse(!is.null(w), paste0("Wilcoxon p=", signif(w$p.value,3)), "NA"),
       x="", y="scran residual variance") + theme_bw() + theme(legend.position = "none")

ggsave(file.path(opt$outdir, "global_plots.png"), pA, width = 8, height = 6, dpi = 300)

gz <- gzfile(file.path(opt$outdir, "noise_with_responsiveness.csv.gz"), "wb")
write.csv(noise2, gz, row.names = FALSE)
close(gz)

sink(file.path(opt$outdir, "global_summary.txt"))
cat("Responsive gene analysis (per cell type; union as global responsive set)\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Cell types tested:", length(responsive_sets), "\n")
cat("Total responsive genes (union):", length(responsive_genes), "\n")
if (!is.null(w)) { print(w) }
sink()