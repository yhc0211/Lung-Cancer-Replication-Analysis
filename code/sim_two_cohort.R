# Simulation for Figure 2 (revised; optional MAMBA)
# Two-cohort setting
# Methods included:
#   1. Threshold methods
#   2. Proposed method
#   3. Meta-analysis (fixed threshold)
#   4. Meta-analysis + BH
#   5. repfdr
#   6. MAMBA (optional; controlled by runMamba)

setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/sim_two.R")

## =========================
## load libraries
## =========================
library(MASS)
library(mvtnorm)
library(data.table)
library(parallel)
library(bindata)
library(dplyr)
library(ks)
library(csmGmm)
library(here)
library(purrr)
library(repfdr)
## MAMBA is optional and can be turned off below.
## Do not load it here unconditionally, because the package may be unavailable
## or slow on some runs.

## =========================
## command line args
## =========================
args <- commandArgs(trailingOnly = TRUE)
aID  <- as.numeric(args[1])
Snum <- as.numeric(args[2])

## =========================
## source supporting code
## =========================
codePath <- here::here("SupportingCode")
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(file.path(codePath, toBeSourced), source)

## =========================
## output
## =========================
outputDir <- here::here("Lung", "output_revised")


outName <- file.path(
  outputDir,
  paste0("rep_two_sim_Snum", Snum, "_aID", aID, ".txt")
)

fpOutName <- file.path(
  outputDir,
  paste0("rep_two_fpbreak_Snum", Snum, "_aID", aID, ".txt")
)

rejOutName <- file.path(
  outputDir,
  paste0("rep_two_rejtruth_Snum", Snum, "_aID", aID, ".txt")
)



## =========================
## options
## =========================
loadData <- FALSE
saveData <- FALSE
testStatsName <- here::here(outputDir, "rep_sim_allZ")

## Set this to TRUE only when you want to run MAMBA.
## When FALSE, all MAMBA-specific outputs remain NA and MAMBA rows are omitted
runMamba <- FALSE

## =========================
## simulation parameters
## =========================
qvalue <- 0.1
case   <- Snum * 9
nSNPs  <- 6744757
K      <- 2

## =========================
## scenario proportions
## =========================
sProp1 <- c(0.9,    0.00675, 0.00675, 0.04305, 0.04305, 0.0001, 0.0001, 0.0001, 0.0001)
sProp2 <- c(0.8,    0.03085, 0.03085, 0.06715, 0.06715, 0.0010, 0.0010, 0.0010, 0.0010)
sProp3 <- c(0.6,    0.07685, 0.07685, 0.11315, 0.11315, 0.0050, 0.0050, 0.0050, 0.0050)

sProp <- c(sProp1, sProp2, sProp3)

## =========================
## build Hmat
## =========================
Hmat <- expand.grid(rep(list(-1:1), K))

blVec <- rep(0, nrow(Hmat))
slVec <- rep(0, nrow(Hmat))
for (k_it in K:1) {
  blVec <- blVec + 2^(K - k_it) * abs(Hmat[, k_it])
  slVec <- slVec + abs(Hmat[, k_it])
}

symAltVec <- ifelse(
  apply(Hmat[, 1:K], 1, sum) == K | apply(Hmat[, 1:K], 1, sum) == -K,
  1, 0
)

Hmat <- Hmat %>%
  mutate(
    bl = blVec,
    sl = slVec,
    symAlt = symAltVec
  ) %>%
  arrange(bl, Var1, Var2) %>%
  mutate(l = 0:(nrow(.) - 1)) %>%
  relocate(l, .before = symAlt)

Hmat <- Hmat %>%
  bind_rows(Hmat, Hmat) %>%
  cbind(sProp) %>%
  dplyr::select(Var1, Var2, sProp)

mu1 <- c(Hmat$Var1 * 2,
         Hmat$Var1 * 3,
         Hmat$Var1 * 4,
         Hmat$Var1 * 5,
         Hmat$Var1 * 6)

mu2 <- c(Hmat$Var2 * 2,
         Hmat$Var2 * 3,
         Hmat$Var2 * 4,
         Hmat$Var2 * 5,
         Hmat$Var2 * 6)

Hmat <- Hmat %>%
  bind_rows(Hmat, Hmat, Hmat, Hmat) %>%
  bind_cols(mu1, mu2) %>%
  rename(mu1 = "...4", mu2 = "...5") %>%
  mutate(number = as.integer(nSNPs * sProp))

Hmat <- Hmat[(case - 8):case, ]

## =========================
## define truth and causal
## replicable = same nonzero direction in both studies
## =========================
Hmat <- Hmat %>%
  mutate(
    truth = case_when(
      Var1 == 0 & Var2 == 0 ~ "null",
      Var1 != 0 & Var2 != 0 & Var1 == Var2 ~ "replicable",
      Var1 != 0 & Var2 != 0 & Var1 != Var2 ~ "discordant",
      xor(Var1 == 0, Var2 == 0) ~ "partial"
    ),
    causal = as.integer(truth == "replicable")
  )

## =========================
## record results
## =========================
powerRes <- data.frame(
  nCausal = 0,
  propCausal = 0,
  mu = 0,
  seed = NA,
  
  nRejNew = NA,
  nRejMeta = NA,
  nRejMetaBH = NA,
  nRejRepfdr = NA,
  nRejMamba = NA,
  nRejThreshold1 = NA,
  nRejThreshold2 = NA,
  nRejThreshold3 = NA,
  
  powerNew = NA,
  powerMeta = NA,
  powerMetaBH = NA,
  powerRepfdr = NA,
  powerMamba = NA,
  powerThreshold1 = NA,
  powerThreshold2 = NA,
  powerThreshold3 = NA,
  
  fdpNew = NA,
  fdpMeta = NA,
  fdpMetaBH = NA,
  fdpRepfdr = NA,
  fdpMamba = NA,
  fdpThreshold1 = NA,
  fdpThreshold2 = NA,
  fdpThreshold3 = NA,
  
  inconNew = NA,
  inconRepfdr = NA,
  inconMamba = NA
)

## =========================
## set seed and scenario-level summaries
## =========================
set.seed(aID * 10^3)
powerRes$seed <- aID * 10^3
powerRes$nCausal <- sum(Hmat$number[Hmat$causal == 1])
powerRes$propCausal <- sum(Hmat$sProp[Hmat$causal == 1])
powerRes$mu <- Hmat$mu1[nrow(Hmat)]

## =========================
## simulate Z
## =========================
allZ <- NULL
for (i in 1:nrow(Hmat)) {
  tempZ <- MASS::mvrnorm(
    n = Hmat$number[i],
    mu = c(Hmat$mu1[i], Hmat$mu2[i]),
    Sigma = diag(K)
  )
  allZ <- rbind(allZ, tempZ)
}

## truncate extreme Z to avoid p-values = 0
for (col_it in 1:ncol(allZ)) {
  tooBig <- which(allZ[, col_it] > 8.1)
  tooSmall <- which(allZ[, col_it] < -8.1)
  if (length(tooBig) > 0) allZ[tooBig, col_it] <- 8.1
  if (length(tooSmall) > 0) allZ[tooSmall, col_it] <- -8.1
}

colnames(allZ) <- c("Z1", "Z2")

## SNP-level truth
truthVec  <- rep(Hmat$truth, Hmat$number)
causalVec <- rep(Hmat$causal, Hmat$number)

## p-values
allP <- 1 - pchisq(as.matrix(allZ)^2, df = 1)

## =========================
## create beta-hat / SE for MAMBA
## =========================
if (runMamba) {
  summaryStatDir <- here::here("Data")
  real_dat  <- fread(here::here(summaryStatDir, "lc_overall_three.txt"))
  se_dat <- real_dat %>%
    dplyr::select(SE1_real = SE,
                  SE2_real = SE_ukb) %>%
    dplyr::filter(is.finite(SE1_real), is.finite(SE2_real),
                  SE1_real > 0, SE2_real > 0)
  rm(real_dat)
  set.seed(aID * 10^3 + 99)
  idx_se <- sample(seq_len(nrow(se_dat)), size = nrow(allZ), replace = TRUE)
  
  SE1 <- se_dat$SE1_real[idx_se]
  SE2 <- se_dat$SE2_real[idx_se]
  
  B1 <- allZ[, 1] * SE1
  B2 <- allZ[, 2] * SE2
}

totOut <- data.frame(
  Z1 = allZ[, 1],
  Z2 = allZ[, 2],
  P1 = allP[, 1],
  P2 = allP[, 2],
  origIdx = seq_len(nrow(allP)),
  truth = truthVec,
  causal = causalVec
)

if (runMamba) {
  totOut$B1  <- B1
  totOut$B2  <- B2
  totOut$SE1 <- SE1
  totOut$SE2 <- SE2
}

## =========================
## helper functions
## =========================
get_metrics <- function(sig, causalVec) {
  nRej <- sum(sig == 1)
  nFP  <- sum(sig == 1 & causalVec == 0)
  nTP  <- sum(sig == 1 & causalVec == 1)
  
  data.frame(
    nRej = nRej,
    fdp = ifelse(nRej == 0, 0, nFP / nRej),
    power = nTP / sum(causalVec)
  )
}

fp_breakdown <- function(sig, truth, method_name, Snum, aID) {
  idx <- which(sig == 1 & truth != "replicable")
  levs <- c("null", "partial", "discordant")
  
  if (length(idx) == 0) {
    return(data.frame(
      method = method_name,
      truth = levs,
      n = c(0, 0, 0),
      prop_among_fp = c(0, 0, 0),
      Snum = Snum,
      aID = aID
    ))
  }
  
  tb <- table(factor(truth[idx], levels = levs))
  data.frame(
    method = method_name,
    truth = levs,
    n = as.numeric(tb),
    prop_among_fp = as.numeric(tb) / sum(tb),
    Snum = Snum,
    aID = aID
  )
}

rej_by_truth <- function(sig, truth, method_name, Snum, aID) {
  tmp <- aggregate(as.numeric(sig) ~ truth, FUN = mean)
  names(tmp)[2] <- "rej_rate"
  tmp$method <- method_name
  tmp$Snum <- Snum
  tmp$aID <- aID
  tmp[, c("method", "truth", "rej_rate", "Snum", "aID")]
}


## =========================
## 1. Threshold methods
## =========================
same_dir_idx <- (totOut$Z1 > 0 & totOut$Z2 > 0) | (totOut$Z1 < 0 & totOut$Z2 < 0)

totOut <- totOut %>%
  mutate(
    threshold1 = as.integer(P1 < 1e-8 & P2 < 1e-8 & same_dir_idx),
    threshold2 = as.integer(P1 < 1e-6 & P2 < 1e-6 & same_dir_idx),
    threshold3 = as.integer(P1 < 1e-5 & P2 < 1e-5 & same_dir_idx)
  )

tmp <- get_metrics(totOut$threshold1, totOut$causal)
powerRes$fdpThreshold1   <- tmp$fdp
powerRes$powerThreshold1 <- tmp$power
powerRes$nRejThreshold1  <- tmp$nRej

tmp <- get_metrics(totOut$threshold2, totOut$causal)
powerRes$fdpThreshold2   <- tmp$fdp
powerRes$powerThreshold2 <- tmp$power
powerRes$nRejThreshold2  <- tmp$nRej

tmp <- get_metrics(totOut$threshold3, totOut$causal)
powerRes$fdpThreshold3   <- tmp$fdp
powerRes$powerThreshold3 <- tmp$power
powerRes$nRejThreshold3  <- tmp$nRej

## =========================
## 2. Proposed method
## =========================
zmat <- as.matrix(totOut[, c("Z1", "Z2")])
initPiList <- list(
  c(0.82),
  c(0.02, 0.02),
  c(0.02, 0.02),
  c(0.1)
)

initMuList <- list(
  matrix(data = 0, nrow = 2, ncol = 1),
  matrix(data = c(0, 3, 0, 6), nrow = 2),
  matrix(data = c(3, 0, 6, 0), nrow = 2),
  matrix(data = c(8, 8), nrow = 2)
)

newRes <- symm_fit_ind_EM(
  testStats = allZ,
  initMuList = initMuList,
  initPiList = initPiList,
  sameDirAlt = TRUE,
  eps = 1e-5
)

totOut <- totOut %>%
  mutate(newLfdr = newRes$lfdrResults) %>%
  arrange(newLfdr) %>%
  mutate(newAvg = cummean(newLfdr)) %>%
  arrange(origIdx) %>%
  mutate(newSig = as.integer(newAvg < qvalue))

tmp <- get_metrics(totOut$newSig, totOut$causal)
powerRes$fdpNew   <- tmp$fdp
powerRes$powerNew <- tmp$power
powerRes$nRejNew  <- tmp$nRej
powerRes$inconNew <- length(check_incongruous(zMatrix = zmat, lfdrVec = totOut$newLfdr))
## =========================
## 3. Meta-analysis
## =========================
w1 <- sqrt(85716)
w2 <- sqrt(332422)

totOut <- totOut %>%
  mutate(
    Z_meta    = (Z1 * w1 + Z2 * w2) / sqrt(w1^2 + w2^2),
    P_meta    = 2 * pnorm(abs(Z_meta), lower.tail = FALSE),
    meta      = as.integer(P_meta < 5e-8),
    P_meta_BH = p.adjust(P_meta, method = "BH"),
    meta_bh   = as.integer(P_meta_BH <= qvalue)
  )

tmp <- get_metrics(totOut$meta, totOut$causal)
powerRes$fdpMeta   <- tmp$fdp
powerRes$powerMeta <- tmp$power
powerRes$nRejMeta  <- tmp$nRej

tmp <- get_metrics(totOut$meta_bh, totOut$causal)
powerRes$fdpMetaBH   <- tmp$fdp
powerRes$powerMetaBH <- tmp$power
powerRes$nRejMetaBH  <- tmp$nRej

## =========================
## 4. repfdr
## =========================


bin_out <- ztobins(
  zmat,
  df = 15,
  plot.diagnostics = FALSE
)

H <- hconfigs(
  n.studies = 2,
  n.association.status = 3
)

non.null.rows <- which(
  apply(H, 1, function(x) all(x == -1) || all(x == 1))
)

rep_out <- repfdr(
  bin_out$pdf.binned.z,
  bin_out$binned.z.mat,
  non.null = "user.defined",
  non.null.rows = non.null.rows,
  control = em.control(
    tol = 1e-4,
    verbose = FALSE
  )
)

totOut <- totOut %>%
  mutate(
    repfdr_fdr = rep_out$mat[, "fdr"],
    repfdr_Fdr = rep_out$mat[, "Fdr"],
    repfdr_sig = as.integer(repfdr_Fdr <= qvalue)
  )

tmp <- get_metrics(totOut$repfdr_sig, totOut$causal)

powerRes$fdpRepfdr   <- tmp$fdp
powerRes$powerRepfdr <- tmp$power
powerRes$nRejRepfdr  <- tmp$nRej

powerRes$inconRepfdr <- length(
  check_incongruous(
    zMatrix = zmat,
    lfdrVec = totOut$repfdr_fdr
  )
)

## =========================
## 5. MAMBA, optional
## =========================
if (runMamba) {
  requireNamespace("mamba", quietly = TRUE)
  
  betajk <- as.matrix(totOut[, c("B1", "B2")])
  sjk2   <- as.matrix(totOut[, c("SE1", "SE2")])^2
  fit_mamba        <- mamba::mamba(betajk = betajk, sjk2 = sjk2)
  totOut$mamba_ppr  <- fit_mamba$ppr
  totOut$mamba_lfdr <- 1 - totOut$mamba_ppr
  
  ppr_rank <- rank(-totOut$mamba_ppr, ties.method = "max")
  ord      <- order(-totOut$mamba_ppr)
  Fdr_hat  <- cumsum(1 - totOut$mamba_ppr[ord]) / ppr_rank[ord]
  
  # Official monotonicity fix: max within tied groups only
  totOut$mamba_Fdr <- NA_real_
  totOut$mamba_Fdr[ord] <- Fdr_hat
  totOut$mamba_Fdr <- ave(totOut$mamba_Fdr, ppr_rank, FUN = max)
  
  totOut$mamba_sig <- as.integer(totOut$mamba_Fdr <= qvalue)
  
  tmp <- get_metrics(totOut$mamba_sig, totOut$causal)
  powerRes$fdpMamba   <- tmp$fdp
  powerRes$powerMamba <- tmp$power
  powerRes$nRejMamba  <- tmp$nRej
  incon_mamba <- get_incon_table(
    zmat = as.matrix(totOut[, c("Z1", "Z2")]),
    lfdrVec = totOut$mamba_lfdr,
    method_name = "MAMBA",
    Snum = Snum,
    aID = aID
  )
  
  powerRes$inconMamba <- nrow(incon_mamba)
}

## =========================
## false positive breakdown
## =========================
fpRes <- bind_rows(
  fp_breakdown(totOut$threshold1,  totOut$truth, "Threshold 1e-8", Snum, aID),
  fp_breakdown(totOut$threshold2,  totOut$truth, "Threshold 1e-6", Snum, aID),
  fp_breakdown(totOut$threshold3,  totOut$truth, "Threshold 1e-5", Snum, aID),
  fp_breakdown(totOut$newSig,      totOut$truth, "Proposed",       Snum, aID),
  fp_breakdown(totOut$meta,        totOut$truth, "Meta",           Snum, aID),
  fp_breakdown(totOut$meta_bh,     totOut$truth, "Meta-BH",        Snum, aID),
  fp_breakdown(totOut$repfdr_sig,  totOut$truth, "repfdr",         Snum, aID)
)

if (runMamba) {
  fpRes <- bind_rows(
    fpRes,
    fp_breakdown(totOut$mamba_sig, totOut$truth, "MAMBA", Snum, aID)
  )
}

## =========================
## rejection rate by truth class
## =========================
rejTruthRes <- bind_rows(
  rej_by_truth(totOut$threshold1,  totOut$truth, "Threshold 1e-8", Snum, aID),
  rej_by_truth(totOut$threshold2,  totOut$truth, "Threshold 1e-6", Snum, aID),
  rej_by_truth(totOut$threshold3,  totOut$truth, "Threshold 1e-5", Snum, aID),
  rej_by_truth(totOut$newSig,      totOut$truth, "Proposed",       Snum, aID),
  rej_by_truth(totOut$meta,        totOut$truth, "Meta",           Snum, aID),
  rej_by_truth(totOut$meta_bh,     totOut$truth, "Meta-BH",        Snum, aID),
  rej_by_truth(totOut$repfdr_sig,  totOut$truth, "repfdr",         Snum, aID)
)

if (runMamba) {
  rejTruthRes <- bind_rows(
    rejTruthRes,
    rej_by_truth(totOut$mamba_sig, totOut$truth, "MAMBA", Snum, aID)
  )
}


## =========================
## save result
## =========================
write.table(
  powerRes,
  file = outName,
  append = FALSE,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

write.table(
  fpRes,
  file = fpOutName,
  append = FALSE,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

write.table(
  rejTruthRes,
  file = rejOutName,
  append = FALSE,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)