#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2); library(patchwork); library(viridisLite)
})

opt <- list(
  make_option("--noise", type="character"),
  make_option("--outdir", type="character"),
  make_option("--min_mean", type="double", default=0.01),
  make_option("--min_pct", type="double", default=0.10)
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(opt$noise, check.names = FALSE)
df <- df %>% mutate(across(where(is.numeric), ~replace(., !is.finite(.), NA_real_)))

base <- df %>% filter(is.finite(mean_expr), is.finite(pct_expressing),
                      mean_expr > opt$min_mean, pct_expressing > opt$min_pct)

p1d <- base %>% filter(is.finite(cv2), cv2 > 0, mean_expr > 0)
p1 <- if (nrow(p1d) > 0) {
  ggplot(p1d, aes(x = log10(mean_expr), y = log10(cv2))) +
    stat_bin2d(bins=60) + scale_fill_viridis_c() +
    geom_smooth(method = "gam", se = FALSE, color = "red") +
    labs(title="Raw noise (CV2) vs mean", x="log10(mean counts)", y="log10(CV2)") +
    theme_bw()
} else ggplot() + theme_void() + labs(title="CV2 panel: insufficient data")

p2d <- base %>% filter(is.finite(scran_resid_var), mean_expr > 0)
p2 <- if (nrow(p2d) > 0) {
  ggplot(p2d, aes(x = log10(mean_expr), y = scran_resid_var)) +
    stat_bin2d(bins=60) + scale_fill_viridis_c() +
    geom_smooth(method = "gam", se = FALSE, color = "red") +
    labs(title="Corrected noise (scran residual variance) vs mean",
         x="log10(mean counts)", y="scran residual variance") +
    theme_bw()
} else ggplot() + theme_void() + labs(title="scran residual panel: insufficient data")

p3d <- base %>% filter(is.finite(scran_resid_var))
p3 <- if (nrow(p3d) > 0) {
  ggplot(p3d, aes(x = scran_resid_var)) +
    geom_histogram(bins=100, fill="skyblue", color="black") +
    geom_density(color="navy") +
    labs(title="Distribution of corrected noise", x="scran residual variance", y="Density") +
    theme_bw()
} else ggplot() + theme_void() + labs(title="No corrected noise available")

p <- (p1 | p2) / p3
ggsave(file.path(opt$outdir, "filter_qc.png"), plot = p, width = 14, height = 10, dpi = 300)

lines <- c(
  "=================================================================",
  "FILTERING & QC SUMMARY",
  "=================================================================",
  paste("Date:", Sys.Date()),
  paste("Base filter:", sprintf("mean_expr > %g & pct_expressing > %g", opt$min_mean, opt$min_pct)),
  paste("Total rows:", nrow(df)),
  paste("Rows after base filter:", nrow(base)),
  paste("Finite scran_resid_var:", sum(is.finite(base$scran_resid_var)))
)
writeLines(lines, file.path(opt$outdir, "filtering_summary.txt"))