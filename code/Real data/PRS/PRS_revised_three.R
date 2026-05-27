# ============================================================
# Lung cancer PRS — three SNP sets (ILCCO GWAS / Meta / Replication)
# 1. Build SNP lists for each method
# 2. Subset UKB GDS → PLINK bed/bim/fam per chromosome
# 3. (shell) Run PRSice per chromosome × method
# 4. Combine per-chr scores
# 5. Evaluate on repeated 80/20 splits: AUC (CI), OR/SD, decile ORs
#
# NOTE: all saved files include "three_way" to keep them separate
#       from the earlier two-way version's outputs.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(magrittr)
  library(data.table); library(here); library(stringr)
  library(ggplot2); library(cowplot); library(ggpubr); library(broom)
  library(pROC); library(SeqArray)
})

setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/PGS.R")

codePath <- here::here("SupportingCode")
invisible(lapply(list.files(codePath, "\\.R$", full.names = TRUE), source))

# Paths
outputDir     <- here::here("Lung", "output")
outputDir_Fig <- "/home/ychang11/csmGmm_reproduce/Lung/PGS/Fig"
dataDir       <- here::here("Data")
geno_dir      <- "/home/ychang11/csmGmm_reproduce/Lung/PGS/data"
score_dir     <- geno_dir
gds_template  <- "/rsrch3/home/biostatistics/UK_Biobank/raw/ukb_imp_chr%s_v3_qc.gds"
pheno_path    <- file.path(geno_dir, "pheno.txt")
cov_path      <- file.path(geno_dir, "cov.txt")

# Resampling settings
N_REPS     <- 100
TRAIN_FRAC <- 0.80
SEED0      <- 20250518

# ============================================================
# 1. SNP sets for the three methods
# ============================================================
three_Z <- fread(here::here(dataDir, "lc_overall_three.txt"))
tempNew <- fread(file.path(outputDir, "three_rep_data_aID1_newlfdr.txt"),
                 header = TRUE, data.table = FALSE)

tempDat <- three_Z %>%
  mutate(origIdx = seq_len(n()),
         newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr),
         rejNew = as.integer(cumNew < 0.1)) %>%
  arrange(origIdx)

rep_SNPs <- tempDat %>%
  filter(rejNew == 1) %>%
  transmute(SNP = RS, CHR = chrom, BP = pos,
            A1 = A1_UKB, A2 = A2_UKB, BETA = Beta_ukb, P = pUKB)

gwas_SNPs <- tempDat %>%
  filter(pILCCO < 5e-8) %>%
  transmute(SNP = RS, CHR = chrom, BP = pos,
            A1 = A1_UKB, A2 = A2_UKB, BETA = Beta_ukb, P = pUKB)

calc_neff <- function(n_case, n_ctrl) 4 / (1 / n_case + 1 / n_ctrl)
N_ILCCO <- calc_neff(29266, 56450)
N_MVP   <- calc_neff(10398, 62708)
N_UKB   <- calc_neff(2404, 330018)
w1 <- sqrt(N_ILCCO); w2 <- sqrt(N_MVP); w3 <- sqrt(N_UKB)

meta_SNPs <- three_Z %>%
  mutate(Z_meta = (zILCCO * w1 + zMVP * w2 + zUKB * w3) / sqrt(w1^2 + w2^2 + w3^2),
         P_meta = 2 * pnorm(abs(Z_meta), lower.tail = FALSE)) %>%
  filter(P_meta < 5e-8) %>%
  transmute(SNP = RS, CHR = chrom, BP = pos,
            A1 = A1_UKB, A2 = A2_UKB, BETA = Beta_ukb, P = pUKB)

snp_sets <- list(
  gwas = list(tag = "gwas", label = "GWAS-significant PRS",  snps = gwas_SNPs),
  meta = list(tag = "meta", label = "Meta-analysis PRS",     snps = meta_SNPs),
  rep  = list(tag = "rep",  label = "Replication-based PRS", snps = rep_SNPs)
)
# for (s in snp_sets) {
#   message(sprintf("%-25s : %d SNPs across %d chromosomes",
#                   s$label, nrow(s$snps), dplyr::n_distinct(s$snps$CHR)))
#   fwrite(s$snps,
#          file.path(geno_dir, paste0("three_way_base_", s$tag, ".txt")),
#          sep = "\t", quote = FALSE)
# }

# ============================================================
# 2. GDS → PLINK BED
# ============================================================
white_id <- fread("/rsrch3/home/biostatistics/UK_Biobank/ychang/lung_cancer/ukb_white_british_ids.txt")
white_id$eid <- paste0(white_id$V1, "_", white_id$V1)

export_chr_bed <- function(chr, snp_df, out_prefix, sample_ids) {
  gdsfile <- seqOpen(sprintf(gds_template, chr), allow.duplicate = TRUE)
  on.exit(seqClose(gdsfile), add = TRUE)
  
  gds_df <- data.frame(
    CHR   = as.numeric(seqGetData(gdsfile, "chromosome")),
    BP    = as.numeric(seqGetData(gdsfile, "position")),
    REF   = as.character(seqGetData(gdsfile, "$ref")),
    ALT   = as.character(seqGetData(gdsfile, "$alt")),
    SNP   = seqGetData(gdsfile, "annotation/id"),
    Index = seq_along(seqGetData(gdsfile, "annotation/id")),
    stringsAsFactors = FALSE
  )
  matched <- snp_df %>% inner_join(gds_df, by = "SNP") %>% distinct()
  if (!nrow(matched)) { message("CHR ", chr, ": no matches"); return(invisible(NULL)) }
  
  variant_keep <- seqGetData(gdsfile, "variant.id")[matched$Index]
  seqSetFilter(gdsfile, sample.id = sample_ids, variant.id = variant_keep)
  if (!length(seqGetData(gdsfile, "variant.id"))) return(invisible(NULL))
  
  tryCatch(seqGDS2BED(gdsfile, out.fn = out_prefix, verbose = FALSE),
           error = function(e) message("CHR ", chr, " export error: ", e$message))
}

# for (s in snp_sets) {
#   for (chr in sort(unique(s$snps$CHR))) {
#     export_chr_bed(chr, s$snps,
#                    file.path(geno_dir,
#                              sprintf("%s_three_way_ukb_chr%s_filtered", s$tag, chr)),
#                    white_id$eid)
#   }
# }

# ============================================================
# 3. PRSice (shell) — see PGS_run_prsice.sh
#    Outputs expected: ${tag}_three_way_prs_chr${chr}_r2.all_score
# ============================================================

# ============================================================
# 4. Combine per-chromosome PRS into a single wide table
# ============================================================
combine_chr_scores <- function(tag,
                               pattern_template = "^%s_three_way_prs_chr[0-9]+_r2\\.all_score$") {
  files <- list.files(score_dir,
                      pattern = sprintf(pattern_template, tag),
                      full.names = TRUE)
  stopifnot("No PRS files found" = length(files) > 0)
  dt <- fread(files[1]); setnames(dt, 1:2, c("FID","IID"))
  for (j in 3:ncol(dt)) set(dt, j = j, value = as.numeric(dt[[j]]))
  for (f in files[-1]) {
    d <- fread(f); setnames(d, 1:2, c("FID","IID"))
    for (j in 3:ncol(d)) set(d, j = j, value = as.numeric(d[[j]]))
    dt <- merge(dt, d, by = c("FID","IID"), all = TRUE, suffixes = c("", ".y"))
    for (y in grep("\\.y$", names(dt), value = TRUE)) {
      x <- sub("\\.y$", "", y)
      if (x %in% names(dt)) {
        dt[[x]][is.na(dt[[x]])] <- 0
        dt[[y]][is.na(dt[[y]])] <- 0
        dt[[x]] <- dt[[x]] + dt[[y]]; dt[[y]] <- NULL
      } else setnames(dt, y, x)
    }
  }
  for (nm in setdiff(names(dt), c("FID","IID"))) dt[[nm]][is.na(dt[[nm]])] <- 0
  dt
}

pull_prs <- function(dt, tag, thresh = "Pt_0.05") {
  stopifnot(thresh %in% names(dt))
  out <- dt[, .(IID, score = get(thresh))]
  setnames(out, "score", paste0(tag, "_PRS"))
  out
}

ph  <- fread(pheno_path)[, .(FID, IID, Phenotype)]
cov <- fread(cov_path)

prs_list <- lapply(snp_sets, function(s) pull_prs(combine_chr_scores(s$tag), s$tag))
prs_wide <- Reduce(function(a, b) merge(a, b, by = "IID", all = FALSE), prs_list)

final <- prs_wide %>%
  inner_join(ph,  by = "IID") %>%
  inner_join(cov, by = "IID") %>%
  as.data.frame()

# Keep complete cases on all evaluation columns
prs_cols <- paste0(vapply(snp_sets, `[[`, "", "tag"), "_PRS")
need <- c("Phenotype", "age", "sex", paste0("PC", 1:10), prs_cols)
final <- final[complete.cases(final[, need]), ]
final$Phenotype <- as.integer(final$Phenotype)
if (!is.numeric(final$sex)) final$sex <- factor(final$sex)

cat(sprintf("Evaluation set: N = %d (cases = %d, controls = %d)\n",
            nrow(final), sum(final$Phenotype == 1), sum(final$Phenotype == 0)))

# ============================================================
# 5. Repeated 80/20 split evaluation
# ============================================================
covar_terms <- c("age", "sex", paste0("PC", 1:10))
f_base <- as.formula(paste("Phenotype ~", paste(covar_terms, collapse = " + ")))

# z-score on TRAIN only, apply same mean/sd to TEST
zscore_train_test <- function(train_vec, test_vec) {
  mu <- mean(train_vec, na.rm = TRUE); sg <- sd(train_vec, na.rm = TRUE)
  if (!is.finite(sg) || sg == 0) sg <- 1
  list(train = (train_vec - mu) / sg, test = (test_vec - mu) / sg)
}

# Stratified split index by Phenotype to keep case prevalence stable
stratified_split <- function(y, train_frac, seed) {
  set.seed(seed)
  idx <- seq_along(y)
  i1  <- idx[y == 1]; i0 <- idx[y == 0]
  s1  <- sample(i1, floor(length(i1) * train_frac))
  s0  <- sample(i0, floor(length(i0) * train_frac))
  list(train = sort(c(s1, s0)), test = sort(setdiff(idx, c(s1, s0))))
}

tags   <- vapply(snp_sets, `[[`, "", "tag")
labels <- vapply(snp_sets, `[[`, "", "label")

# Storage
auc_mat <- matrix(NA_real_, nrow = N_REPS, ncol = length(tags) + 1,
                  dimnames = list(NULL, c("base", tags)))
or_mat  <- matrix(NA_real_, nrow = N_REPS, ncol = length(tags),
                  dimnames = list(NULL, tags))
delong_p <- matrix(NA_real_, nrow = N_REPS, ncol = length(tags),
                   dimnames = list(NULL, tags))
decile_records <- vector("list", N_REPS)

for (r in seq_len(N_REPS)) {
  spl <- stratified_split(final$Phenotype, TRAIN_FRAC, SEED0 + r)
  tr  <- final[spl$train, , drop = FALSE]
  te  <- final[spl$test,  , drop = FALSE]
  
  # z-score every PRS using TRAIN moments
  for (tg in tags) {
    col <- paste0(tg, "_PRS")
    z   <- zscore_train_test(tr[[col]], te[[col]])
    tr[[paste0(col, "_z")]] <- z$train
    te[[paste0(col, "_z")]] <- z$test
  }
  
  # base model (covariates only)
  m_base   <- glm(f_base, data = tr, family = binomial())
  p_te     <- predict(m_base, newdata = te, type = "response")
  roc_base <- roc(te$Phenotype, p_te, quiet = TRUE)
  auc_mat[r, "base"] <- as.numeric(auc(roc_base))
  
  rocs_full <- list()
  for (tg in tags) {
    z_col  <- paste0(tg, "_PRS_z")
    m_full <- glm(update(f_base, as.formula(paste(". ~ . +", z_col))),
                  data = tr, family = binomial())
    p_te_f <- predict(m_full, newdata = te, type = "response")
    roc_f  <- roc(te$Phenotype, p_te_f, quiet = TRUE)
    auc_mat[r, tg] <- as.numeric(auc(roc_f))
    or_mat[r, tg]  <- exp(coef(m_full)[z_col])
    rocs_full[[tg]] <- roc_f
  }
  
  # DeLong test: each PRS-augmented model vs covariates-only base, on TEST
  for (tg in tags) {
    delong_p[r, tg] <- roc.test(rocs_full[[tg]], roc_base,
                                method = "delong")$p.value
  }
  
  # decile ORs on TEST
  per_rep <- lapply(tags, function(tg) {
    z_col <- paste0(tg, "_PRS_z")
    te$decile <- ntile(te[[z_col]], 10)
    fml <- as.formula(paste("Phenotype ~ factor(decile) +",
                            paste(covar_terms, collapse = " + ")))
    fit <- glm(fml, data = te, family = binomial())
    td  <- tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    or_tbl <- td %>%
      filter(grepl("^factor\\(decile\\)", term)) %>%
      transmute(decile = as.integer(gsub("factor\\(decile\\)", "", term)),
                or = estimate, lo = conf.low, hi = conf.high) %>%
      bind_rows(tibble(decile = 1L, or = 1, lo = 1, hi = 1)) %>%
      mutate(tag = tg, rep = r)
    prev <- te %>%
      mutate(decile = ntile(.data[[z_col]], 10)) %>%
      group_by(decile) %>%
      summarise(n = n(),
                cases = sum(Phenotype == 1),
                prev = cases / n, .groups = "drop") %>%
      mutate(tag = tg, rep = r)
    list(or = or_tbl, prev = prev)
  })
  decile_records[[r]] <- list(
    or   = bind_rows(lapply(per_rep, `[[`, "or")),
    prev = bind_rows(lapply(per_rep, `[[`, "prev"))
  )
  
  if (r %% 10 == 0) message("rep ", r, "/", N_REPS, " done")
}

# ============================================================
# 6. Summarize across reps
# ============================================================
ci <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)

# DeLong: median p across reps + fraction of reps with p < 0.05
delong_med <- apply(delong_p, 2, median, na.rm = TRUE)
delong_sig <- apply(delong_p, 2, function(x) mean(x < 0.05, na.rm = TRUE))

summary_tab <- data.frame(
  Method   = c("Covariates only", labels),
  N_SNPs   = c(NA, vapply(snp_sets, function(s) nrow(s$snps), 0L)),
  AUC_mean = c(mean(auc_mat[, "base"]),
               sapply(tags, function(tg) mean(auc_mat[, tg]))),
  AUC_lo   = c(ci(auc_mat[, "base"])[1],
               sapply(tags, function(tg) ci(auc_mat[, tg])[1])),
  AUC_hi   = c(ci(auc_mat[, "base"])[3],
               sapply(tags, function(tg) ci(auc_mat[, tg])[3])),
  Delta_AUC_mean = c(0,
                     sapply(tags, function(tg) mean(auc_mat[, tg] - auc_mat[, "base"]))),
  OR_per_SD_mean = c(NA,
                     sapply(tags, function(tg) mean(or_mat[, tg]))),
  OR_lo = c(NA, sapply(tags, function(tg) ci(or_mat[, tg])[1])),
  OR_hi = c(NA, sapply(tags, function(tg) ci(or_mat[, tg])[3])),
  DeLong_median_p = c(NA, delong_med[tags]),
  DeLong_frac_p05 = c(NA, delong_sig[tags]),
  row.names = NULL
)
print(summary_tab, row.names = FALSE, digits = 4)
fwrite(summary_tab, file.path(outputDir_Fig, "three_way_AUC_summary_cv.csv"))

cat("\nDeLong median p (PRS vs covariates-only):\n"); print(round(delong_med, 4))
cat("\nFraction of reps with DeLong p<0.05:\n");      print(round(delong_sig, 3))

# Save raw rep-level matrices too (handy for later plots / sensitivity checks)
fwrite(as.data.frame(auc_mat),
       file.path(outputDir_Fig, "three_way_AUC_per_rep.csv"))
fwrite(as.data.frame(or_mat),
       file.path(outputDir_Fig, "three_way_OR_per_rep.csv"))
fwrite(as.data.frame(delong_p),
       file.path(outputDir_Fig, "three_way_DeLong_per_rep.csv"))

# Pool decile records (mean OR and prevalence across reps)
or_all <- bind_rows(lapply(decile_records, `[[`, "or"))

or_df <- or_all %>%
  group_by(tag, decile) %>%
  summarise(or_mean = mean(or, na.rm = TRUE),
            lo      = quantile(or, 0.025, na.rm = TRUE),
            hi      = quantile(or, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  rename(or = or_mean) %>%
  mutate(PRS = factor(setNames(labels, tags)[tag], levels = labels))

fwrite(or_df, file.path(outputDir_Fig, "three_way_decile_OR_cv.csv"))

message("Done. Outputs in: ", outputDir_Fig)