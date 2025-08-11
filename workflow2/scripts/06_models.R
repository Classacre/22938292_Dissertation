#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2); library(lme4); library(MatchIt); library(metafor); library(patchwork)
})

opt <- list(
  make_option("--csv", type="character"),
  make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(opt$csv)
df <- df %>% filter(is.finite(scran_resid_var), mean_expr > 0.01, pct_expressing > 0.10, !is.na(cahn_group))

# Mixed-effects model across identities and genes
df$cahn_group <- factor(df$cahn_group, levels = c("Unmethylated","gbM","TE-like"))
df$h2az_group <- factor(df$h2az_group, levels = c("Other","H2A.Z-Enriched","H2A.Z-Depleted"))

m1 <- lmer(scran_resid_var ~ log10(mean_expr) + pct_expressing + has_TATA_box + h2az_group + cahn_group + (1|identity) + (1|gene), data = df, REML = TRUE)

# Per-identity contrasts: gbM vs Unmethylated, adjusted
effect_list <- list()
idents <- unique(df$identity)
for (id in idents) {
  di <- df %>% filter(identity == id, cahn_group %in% c("Unmethylated","gbM"))
  if (nrow(di) < 200 || length(unique(di$cahn_group)) < 2) next
  fit <- tryCatch(lm(scran_resid_var ~ log10(mean_expr) + pct_expressing + has_TATA_box + h2az_group + cahn_group, data = di), error=function(e) NULL)
  if (is.null(fit)) next
  co <- coef(summary(fit))
  if ("cahn_groupgbM" %in% rownames(co)) {
    est <- co["cahn_groupgbM","Estimate"]; se <- co["cahn_groupgbM","Std. Error"]
    effect_list[[id]] <- data.frame(identity = id, yi = est, sei = se, stringsAsFactors = FALSE)
  }
}
eff <- dplyr::bind_rows(effect_list)
meta <- if (nrow(eff) >= 2) tryCatch(rma(yi = eff$yi, sei = eff$sei, method = "REML"), error=function(e) NULL) else NULL

# Matching analysis (propensity on covariates) within identity, pooled paired diffs
match_effects <- list()
for (id in idents) {
  di <- df %>% filter(identity == id, cahn_group %in% c("Unmethylated","gbM")) %>%
    mutate(treat = as.integer(cahn_group == "gbM"))
  if (nrow(di) < 200 || length(unique(di$treat)) < 2) next
  # propensity: mean_expr, pct_expressing, has_TATA_box, h2az_group
  m <- tryCatch(MatchIt::matchit(treat ~ log10(mean_expr) + pct_expressing + has_TATA_box + h2az_group,
                                 data = di, method = "nearest", ratio = 1), error=function(e) NULL)
  if (is.null(m)) next
  md <- MatchIt::match.data(m)
  # paired difference by nearest neighbor groups; approximate via regression on matched set
  fit <- lm(scran_resid_var ~ treat + log10(mean_expr) + pct_expressing + has_TATA_box + h2az_group, data = md)
  co <- coef(summary(fit))
  if ("treat" %in% rownames(co)) {
    est <- co["treat","Estimate"]; se <- co["treat","Std. Error"]
    match_effects[[id]] <- data.frame(identity = id, yi = est, sei = se, stringsAsFactors = FALSE)
  }
}
meff <- dplyr::bind_rows(match_effects)
meta_match <- if (nrow(meff) >= 2) tryCatch(rma(yi = meff$yi, sei = meff$sei, method = "REML"), error=function(e) NULL) else NULL

# Figure: global distributions
p1 <- ggplot(df, aes(x = cahn_group, y = scran_resid_var, fill = cahn_group)) +
  geom_violin(trim=FALSE) + geom_boxplot(width=0.12, outlier.shape = NA, fill="white") +
  scale_fill_manual(values=c("Unmethylated"="grey70","gbM"="#0072B2","TE-like"="#D55E00")) +
  labs(title="Noise by methylation group", x="", y="scran residual variance") + theme_bw() + theme(legend.position="none")

ggsave(file.path(opt$outdir, "global_figure.png"), p1, width = 8, height = 6, dpi = 300)

# Write summaries
sink(file.path(opt$outdir, "global_model_summary.txt"))
cat("Global mixed-effects model on identity-specific corrected noise\n")
cat("Date:", as.character(Sys.Date()), "\n\n")
print(summary(m1))
sink()

sink(file.path(opt$outdir, "meta_analysis.txt"))
cat("Per-identity gbM vs Unmethylated contrasts (adjusted), random-effects meta-analysis\n")
if (!is.null(meta)) {
  print(meta)
} else {
  cat("Meta-analysis could not be computed (insufficient identities).\n")
}
cat("\nMatching-based contrasts meta-analysis\n")
if (!is.null(meta_match)) {
  print(meta_match)
} else {
  cat("Matching meta-analysis could not be computed.\n")
}
sink()