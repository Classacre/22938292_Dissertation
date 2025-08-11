#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2); library(purrr); library(broom); library(readr)
})

opt <- list(
  make_option("--csv", type="character"),
  make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(opt$csv) %>%
  filter(is.finite(scran_resid_var), mean_expr > 0.01, pct_expressing > 0.10, !is.na(cahn_group))

agg <- df %>%
  group_by(gene, cahn_group, has_TATA_box) %>%
  summarise(
    scran_resid_var = mean(scran_resid_var, na.rm=TRUE),
    log10_mean_expr = log10(mean(mean_expr, na.rm=TRUE)),
    .groups = "drop"
  ) %>% filter(is.finite(log10_mean_expr))
agg$cahn_group <- factor(agg$cahn_group, levels = c("Unmethylated","gbM","TE-like"))

N <- 1000
safe_lm <- safely(lm)
boot <- map_dfr(seq_len(N), function(i) {
  bs <- dplyr::sample_n(agg, size = nrow(agg), replace = TRUE)
  fit <- safe_lm(scran_resid_var ~ log10_mean_expr + has_TATA_box + cahn_group, data = bs)
  if (is.null(fit$error)) {
    broom::tidy(fit$result) %>% mutate(bootstrap_id = i, .before = 1)
  } else NULL
})

sumtab <- boot %>% group_by(term) %>%
  summarise(
    n_success = n(),
    mean_estimate = mean(estimate, na.rm=TRUE),
    std_error = sd(estimate, na.rm=TRUE),
    ci_2.5 = quantile(estimate, 0.025, na.rm=TRUE),
    ci_97.5 = quantile(estimate, 0.975, na.rm=TRUE),
    .groups = "drop"
  )

p <- ggplot(boot, aes(x = estimate)) +
  geom_density(fill="skyblue", alpha=0.7) +
  geom_vline(data = sumtab, aes(xintercept = mean_estimate), color="red", linetype="dashed") +
  geom_vline(data = sumtab, aes(xintercept = ci_2.5), color="black", linetype="dotted") +
  geom_vline(data = sumtab, aes(xintercept = ci_97.5), color="black", linetype="dotted") +
  facet_wrap(~ term, scales = "free", ncol = 2) +
  labs(title="Bootstrap distributions of coefficients", subtitle=paste(N, "iterations"),
       x="Estimate", y="Density") + theme_bw()
ggsave(file.path(opt$outdir, "bootstrap_coeffs.png"), p, width=10, height=12, dpi=300)

write.csv(sumtab, file.path(opt$outdir, "bootstrap_table.csv"), row.names = FALSE)

sink(file.path(opt$outdir, "bootstrap_summary.txt"))
cat("Bootstrap summary\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Iterations:", N, "\n\n")
print(sumtab, n = nrow(sumtab))
sink()