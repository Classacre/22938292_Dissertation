# Combined Methylation Data README

This directory contains a combined methylation annotation file generated from two studies:

- **Cahn_et_al_2024.csv**: Contains gene methylation status as classified by Cahn et al. (2024). Columns: `Gene_ID`, `Methylation_status` (values: gbM, TE-like methylation, Unmethylated).
- **Bewick_et_al_2016.csv**: Contains gene methylation status as classified by Bewick et al. (2016). Columns: `Gene`, `Classification` (values: gbM, mCHG, mCHH, Unmethylated), plus methylation percentages (not included in the combined file).

## Combined File

- **combined_methylation_data.csv**: Merges the above datasets by gene, with the following columns:
  - `Gene_ID`: Standardized gene identifier (no isoform suffix)
  - `Cahn_Methylation_status`: Methylation status from Cahn et al. (gbM, TE-like methylation, Unmethylated)
  - `Bewick_Classification`: Methylation status from Bewick et al. (gbM, mCHG, mCHH, Unmethylated)
  - `Both_gbM`: TRUE if both studies classify as gbM, FALSE if not, NA if missing in either
  - `Both_TE_like`: TRUE if Cahn is TE-like methylation and Bewick is mCHG or mCHH, FALSE if not, NA if missing in either
  - `H2AZ_Depleted`: TRUE if gene is in the H2A.Z body-depleted list, FALSE if in enriched, NA if in neither
  - `H2AZ_Enriched`: TRUE if gene is in the H2A.Z body-enriched list, FALSE if in depleted, NA if in neither

## How to Reproduce

Run the script `combine_methylation_datasets.R` in this directory. It will generate the combined file using only the relevant columns and logic described above.

---

This combined file can be used for downstream analyses comparing methylation status assignments between studies, for identifying consensus gbM or TE-like genes, and for integrating with other genomic features such as H2A.Z enrichment.
