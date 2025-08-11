#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(GenomeInfoDb)
  library(GenomicFeatures); library(GenomicRanges); library(BSgenome)
  library(Biostrings); library(AnnotationHub); library(TFBSTools); library(JASPAR2024)
})

opt <- list(
  make_option("--cahn", type="character"),
  make_option("--bewick", type="character"),
  make_option("--h2az_dep", type="character"),
  make_option("--h2az_enr", type="character"),
  make_option("--outdir", type="character"),
  make_option("--genome", type="character"),      # e.g., BSgenome.Athaliana.TAIR.TAIR9
  make_option("--txdb", type="character"),        # e.g., TxDb.Athaliana.BioMart.plantsmart28
  make_option("--use_araport11", type="character", default="true"),
  make_option("--jaspar_id", type="character", default="MA0108.2"),
  make_option("--pwm_min_score", type="character", default="85%")
)
opt <- parse_args(OptionParser(option_list = opt))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# 1) Load methylation and H2A.Z lists
cahn <- read.csv(opt$cahn, stringsAsFactors = FALSE)
colnames(cahn) <- c("Gene_ID", "Cahn_Methylation_status")
bew <- read.csv(opt$bewick, stringsAsFactors = FALSE)
colnames(bew)[1:2] <- c("Gene", "Bewick_Classification")
bew$Gene_ID <- sub("\\..*", "", bew$Gene)
dep <- read.csv(opt$h2az_dep, stringsAsFactors = FALSE)[,1]
enr <- read.csv(opt$h2az_enr, stringsAsFactors = FALSE)[,1]

# 2) Build TxDb and load genome
suppressWarnings({ library(opt$genome, character.only = TRUE) })
genome <- getExportedValue(opt$genome, sub(".*\\.", "", opt$genome)) # e.g., BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9

txdb <- NULL
use_araport11 <- tolower(opt$use_araport11) %in% c("true","t","1","yes","y")
if (use_araport11) {
  message("Attempting to use Araport11 GFF via AnnotationHub...")
  ah <- AnnotationHub::AnnotationHub()
  # Query Araport11 GFF3
  hits <- AnnotationHub::query(ah, c("Athaliana", "Araport11", "GFF3"))
  if (length(hits) > 0) {
    gff <- hits[[1]]
    txdb <- GenomicFeatures::makeTxDbFromGFF(gff)
  } else {
    warning("Araport11 GFF not found in AnnotationHub. Falling back to provided TxDb.")
  }
}
if (is.null(txdb)) {
  suppressWarnings({ library(opt$txdb, character.only = TRUE) })
  txdb <- get(opt$txdb)
}

# 3) Harmonize seqlevels to BSgenome style
mapping <- c("1"="Chr1", "2"="Chr2", "3"="Chr3", "4"="Chr4", "5"="Chr5", "Mt"="ChrM", "Pt"="ChrC")
common <- intersect(seqlevels(txdb), names(mapping))
txdb2 <- txdb
suppressWarnings({
  seqlevelsStyle(txdb2) <- seqlevelsStyle(genome)
  if (length(common) > 0) txdb2 <- renameSeqlevels(txdb2, mapping[names(mapping) %in% seqlevels(txdb2)])
})
seq_keep <- intersect(seqlevels(txdb2), seqlevels(genome))
txdb2 <- keepSeqlevels(txdb2, seq_keep, pruning.mode = "coarse")

# 4) Promoter extraction around TSS and TATA PWM scanning
# Promoters: [-200, +50] from TSS, strand-aware
proms <- promoters(txdb2, upstream = 200, downstream = 50, use.names = TRUE)
proms <- trim(proms)
proms <- keepStandardChromosomes(proms, pruning.mode = "coarse")
proms <- proms[seqnames(proms) %in% paste0("Chr", c(1:5, "M","C"))]

# gene_id retrieval: prefer gene_id if present, else map transcript to gene via transcriptsBy(txdb)
tx2gene <- tryCatch(GenomicFeatures::select(txdb2, keys = names(proms), keytype = "TXNAME", columns = c("GENEID")), error=function(e) NULL)
if (!is.null(tx2gene)) {
  m <- match(names(proms), tx2gene$TXNAME)
  gene_ids <- tx2gene$GENEID[m]
} else {
  gene_ids <- names(proms)
}
proms$gene_id <- gene_ids
proms <- proms[!is.na(proms$gene_id)]
# de-duplicate by gene: keep the most upstream promoter per gene
ord <- order(proms$gene_id)
proms <- proms[ord]
proms <- proms[!duplicated(proms$gene_id)]

seqs <- getSeq(genome, proms)

# JASPAR PWM for TBP
opts <- list()
pfm <- tryCatch({
  getMatrixSet(JASPAR2024, opts = list(ID = opt$jaspar_id))[[1]]
}, error=function(e) NULL)
if (is.null(pfm)) {
  warning("JASPAR PWM not found; using simple TATAWAWR consensus PWM.")
  # Construct simple IUPAC consensus PWM
  pwm <- matrix(0.25, nrow=4, ncol=8, dimnames=list(c("A","C","G","T"), NULL))
  # T A T A W A W R = T/A bias etc; keep it simple
  consensus <- c("T","A","T","A","W","A","W","R")
  iupac_map <- list(
    A=c(1,0,0,0), C=c(0,1,0,0), G=c(0,0,1,0), T=c(0,0,0,1),
    W=c(0.5,0,0,0.5), R=c(0.5,0,0.5,0)
  )
  for (i in seq_along(consensus)) pwm[,i] <- iupac_map[[consensus[i]]]
  pfm <- TFBSTools::PFMatrix(ID="CONS_TATA", name="CONS_TATA", strand="*", matrix=pwm)
}
pwm <- toPWM(pfm, pseudocounts=0.1)

min_score <- opt$pwm_min_score
scores <- vapply(seqs, function(s) {
  m1 <- searchSeq(pwm, s, min.score = min_score, strand = "+")
  m2 <- searchSeq(pwm, s, min.score = min_score, strand = "-")
  as.integer(length(m1) + length(m2))
}, integer(1))
df_tata <- data.frame(Gene_ID = proms$gene_id, has_TATA_box = scores > 0, stringsAsFactors = FALSE)

# 5) Merge all annotations
merged <- full_join(cahn, bew %>% select(Gene_ID, Bewick_Classification), by = "Gene_ID") %>%
  left_join(df_tata, by = "Gene_ID") %>%
  mutate(
    Cahn_group = case_when(
      Cahn_Methylation_status == "gbM" ~ "gbM",
      Cahn_Methylation_status == "TE-like methylation" ~ "TE-like",
      Cahn_Methylation_status == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    ),
    Bewick_group = case_when(
      Bewick_Classification == "gbM" ~ "gbM",
      Bewick_Classification %in% c("mCHG","mCHH") ~ "TE-like",
      Bewick_Classification == "Unmethylated" ~ "Unmethylated",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    H2AZ_Depleted = Gene_ID %in% dep,
    H2AZ_Enriched = Gene_ID %in% enr,
    gbm_intersection = (Cahn_group == "gbM") & (Bewick_group == "gbM"),
    gbm_union = (Cahn_group == "gbM") | (Bewick_group == "gbM")
  ) %>%
  distinct(Gene_ID, .keep_all = TRUE)

write.csv(merged, file.path(opt$outdir, "epigenetic_annotations.csv"), row.names = FALSE)

# 6) Summary
sink(file.path(opt$outdir, "summary.txt"))
cat("Epigenetic annotation summary\n")
cat("Date:", as.character(Sys.Date()), "\n\n")
cat("Input sizes:\n")
cat("  Cahn:", nrow(cahn), "\n")
cat("  Bewick:", nrow(bew), "\n")
cat("  H2A.Z depleted:", length(unique(dep)), "\n")
cat("  H2A.Z enriched:", length(unique(enr)), "\n\n")
cat("Promoters scanned:", nrow(df_tata), " (window -200:+50)\n")
cat("TATA PWM:", opt$jaspar_id, "min.score:", opt$pwm_min_score, "\n\n")
cat("Merged genes:", nrow(merged), "\n")
cat("\nCounts:\n")
print(table(merged$Cahn_group, useNA = "ifany"))
print(table(merged$Bewick_group, useNA = "ifany"))
cat("\nH2A.Z:\n")
cat("  Depleted:", sum(merged$H2AZ_Depleted, na.rm=TRUE), "\n")
cat("  Enriched:", sum(merged$H2AZ_Enriched, na.rm=TRUE), "\n")
cat("  Overlap:", sum(merged$H2AZ_Depleted & merged$H2AZ_Enriched, na.rm=TRUE), "\n")
cat("\nPromoter TATA:\n")
cat("  TRUE:", sum(merged$has_TATA_box, na.rm=TRUE), "\n")
cat("  FALSE:", sum(!merged$has_TATA_box, na.rm=TRUE), "\n")
cat("\nGBM union:", sum(merged$gbm_union, na.rm=TRUE), " intersection:", sum(merged$gbm_intersection, na.rm=TRUE), "\n")
sink()