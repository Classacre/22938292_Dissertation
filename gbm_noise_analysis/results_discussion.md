# Results and Discussion

## Overview
This section presents a detailed analysis of gene expression noise, methylation, and chromatin features in Arabidopsis root single-cell RNA-seq data. The results address all project objectives and test the stated hypotheses, using extensive statistical and visual outputs. All figures referenced are available in the `gbm_noise_analysis/results/high_low_noise/` and `gbm_noise_analysis/results/responsive_genes/` directories.

---

## 1. Identification of High- and Low-Noise Genes

Gene expression noise was quantified as the coefficient of variation (CV) for each gene across all single cells. Genes were ranked by CV, with the top 10% classified as high-noise and the bottom 10% as low-noise. This classification is visualized in:
- **Boxplots:**
  - `gbm_noise_boxplot_cahn.png` (Cahn: gbM, TE-like, Unmethylated)
  - `gbm_noise_boxplot_bewick.png` (Bewick: gbM, mCHG, mCHH, Unmethylated)
  - `gbm_noise_boxplot_combined.png` (all groups)
  - `gbm_noise_boxplot_te_like.png`, `gbm_noise_boxplot_by_celltype.png`, `gbm_noise_boxplot_by_lineage.png`
- **Scatter plots:**
  - `gbm_noise_scatter_cahn.png`, `gbm_noise_scatter_bewick.png`, `gbm_noise_scatter_combined.png`
  - `mean_vs_sd_by_cahn_group.png`, `mean_vs_cv_by_cahn_group.png`
  - `mean_vs_sd_by_bewick_group.png`, `mean_vs_cv_by_bewick_group.png`
  - `mean_vs_sd_by_methylation_group.png`, `mean_vs_cv_by_methylation_group.png`

These plots show that high-noise genes exhibit substantial variability, while low-noise genes are more stably expressed.

---

## 2. Identification and Characterization of Responsive Genes

Responsive genes were defined as those highly expressed in one cell type but low in others. The top 10% by responsiveness score were classified as highly responsive. Visualizations include:
- `responsive_vs_nonresponsive_noise.png` (CV for responsive vs non-responsive genes)
- `cv_by_cahn_and_responsive.png`, `cv_by_bewick_and_responsive.png` (CV by methylation group and responsiveness)
- `cv_by_h2az_and_responsive.png` (CV by H2A.Z group and responsiveness)
- `cv_by_meth_h2az_and_responsive.png` (CV by methylation + H2A.Z group and responsiveness)

These figures demonstrate that highly responsive genes tend to have higher expression noise.

---

## 3. Genomic Feature Analysis: Methylation and H2A.Z Status

### a. Methylation Status (gbM and TE-like)

Genes were annotated with methylation status from both Cahn and Bewick datasets. Key findings:
- **gbM genes:** Consistently exhibited lower expression noise (see `gbm_noise_boxplot_cahn.png`, `gbm_noise_boxplot_bewick.png`).
- **TE-like methylated genes:** Displayed higher variability (see `gbm_noise_boxplot_te_like.png`).
- **Unmethylated genes:** Showed the highest noise.

#### Statistical Results
- **Cahn gbM vs non-gbM:** p = 0
- **Bewick gbM vs non-gbM:** p = 1.08 × 10⁻²¹⁴
- **Bewick TE-like vs other:** p = 2.03 × 10⁻⁷³

All differences are highly significant, confirming the hypothesis that gbM is associated with reduced noise.

### b. H2A.Z Status

Genes were annotated for H2A.Z enrichment or depletion. Results:
- **H2A.Z-depleted genes:** Lower noise (see `gbm_noise_boxplot_h2az.png`).
- **H2A.Z-enriched genes:** Higher noise.
- **Combined methylation + H2A.Z:** Genes with both gbM and H2A.Z depletion had the lowest noise (see `gbm_noise_boxplot_meth_h2az.png`).

#### Statistical Results
- **H2A.Z-Depleted vs H2A.Z-Enriched:** p = 0

---

## 4. Overlap and Shared Features Between Noisy and Responsive Genes

Overlap analysis (see `genomic_features_overlap_summary.csv`):
- Genes both high-noise and highly responsive are enriched for unmethylated and H2A.Z-enriched status.
- Low-noise, low-responsive genes are enriched for gbM and H2A.Z depletion.

Fisher's exact tests (see `genomic_features_fisher_summary.csv`) confirm these enrichments are statistically significant.

---

## 5. Responsive Gene Enrichment and Statistical Results

- **Responsive vs Non-Responsive (CV):** p = 1.08 × 10⁻¹²
- **Fisher's exact (Cahn):** p = 1.56 × 10⁻¹⁷
- **Fisher's exact (Bewick):** p = 3.67 × 10⁻¹³
- **Fisher's exact (H2A.Z):** p = 3.09 × 10⁻⁷

These results show that responsive genes are significantly less likely to be gbM or H2A.Z-depleted, and more likely to be noisy.

---

## 6. Discussion and Interpretation

### a. Hypothesis Evaluation
The results strongly support the hypothesis that gbM and H2A.Z depletion are associated with reduced gene expression noise. All statistical comparisons yielded extremely low p-values, indicating robust differences between groups. The visualizations further reinforce these findings, with clear separation of noise distributions by methylation and H2A.Z status.

Responsive genes, which are cell-type-specific, are more likely to be noisy and less likely to be gbM or H2A.Z-depleted. This suggests that epigenetic stability is a hallmark of housekeeping genes, while regulatory genes are more dynamic.

### b. Broader Implications
The integration of methylation and chromatin features provides a nuanced view of gene regulation in plant single cells. The results highlight the importance of epigenetic mechanisms in stabilizing gene expression and shaping cell identity.

### c. Limitations and Future Directions
- The analysis is limited to a single, unpublished Seurat object. Broader validation is needed.
- Responsive gene classification could be refined with more sophisticated models.
- Future work should explore the mechanistic basis of these associations.

---

## 7. Figures and Tables

All referenced figures are available in the results folders:
- `gbm_noise_analysis/results/high_low_noise/` (noise and methylation plots)
- `gbm_noise_analysis/results/responsive_genes/` (responsive gene plots)
- `gbm_noise_analysis/results/genomic_features/` (overlap and Fisher's test tables)

Key figures include:
- `gbm_noise_boxplot_cahn.png`, `gbm_noise_boxplot_bewick.png`, `gbm_noise_boxplot_h2az.png`, `gbm_noise_boxplot_meth_h2az.png`
- `mean_vs_sd_by_cahn_group.png`, `mean_vs_cv_by_cahn_group.png`, `mean_vs_sd_by_bewick_group.png`, `mean_vs_cv_by_bewick_group.png`, `mean_vs_sd_by_h2az_group.png`, `mean_vs_cv_by_h2az_group.png`
- `responsive_vs_nonresponsive_noise.png`, `cv_by_cahn_and_responsive.png`, `cv_by_bewick_and_responsive.png`, `cv_by_h2az_and_responsive.png`, `cv_by_meth_h2az_and_responsive.png`

Key tables include:
- `gbm_noise_stats_summary.csv`, `responsive_stats_summary.csv`, `genomic_features_overlap_summary.csv`, `genomic_features_fisher_summary.csv`

---

## Conclusion

This study provides a comprehensive, integrative analysis of gene expression noise, methylation, and chromatin features in Arabidopsis root single cells. The results confirm that gbM and H2A.Z depletion are associated with stable, low-noise gene expression, while cell-type-specific and highly variable genes are depleted for these features. These findings advance our understanding of the epigenetic regulation of gene expression variability in plants.
