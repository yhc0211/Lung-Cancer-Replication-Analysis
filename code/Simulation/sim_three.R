# Simulation for Figure 2: Two-cohort simulation
# Three-cohort setting
# Methods included:
#   1. Threshold methods
#   2. Proposed method
#   3. Meta-analysis 
#   4. Meta-analysis + BH
#   5. repfdr


setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/sim_three.R")

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
library(locfdr)

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
  paste0("rep_three_sim_Snum", Snum, "_aID", aID, ".txt")
)

fpOutName <- file.path(
  outputDir,
  paste0("rep_three_fpbreak_Snum", Snum, "_aID", aID, ".txt")
)





## =========================
## options
## =========================
loadData <- FALSE
saveData <- FALSE

## Set to TRUE only when you want to run MAMBA.
## When FALSE, all MAMBA-specific input preparation, fitting, and output rows are skipped.
runMamba <- FALSE

testStatsName <- here::here(outputDir, "rep_sim_allZ")

## =========================
## simulation parameters
## =========================
qvalue <- 0.1
nSNPs  <- 6744757
K      <- 3

nConfig <- 3^K
case    <- Snum * nConfig

## =========================
## scenario proportions
## =========================
sProp1 <- c(
  0.9002,
  rep(0.0055, 18),
  rep(0.0001, 8)
)

sProp2 <- c(
  0.794,
  rep(0.011, 18),
  rep(0.001, 8)
)

sProp3 <- c(
  0.6,
  rep(0.02, 18),
  rep(0.005, 8)
)

## =========================
## build Hmat
## =========================
Hbase <- expand.grid(rep(list(-1:1), K))

blVec <- rep(0, nrow(Hbase))
slVec <- rep(0, nrow(Hbase))

for (k_it in K:1) {
  blVec <- blVec + 2^(K - k_it) * abs(Hbase[, k_it])
  slVec <- slVec + abs(Hbase[, k_it])
}

symAltVec <- ifelse(
  apply(Hbase[, 1:K], 1, sum) == K |
    apply(Hbase[, 1:K], 1, sum) == -K,
  1, 0
)

Hbase <- Hbase %>%
  mutate(
    bl = blVec,
    sl = slVec,
    symAlt = symAltVec
  ) %>%
  arrange(bl, Var1, Var2, Var3) %>%
  mutate(l = 0:(nrow(.) - 1)) %>%
  relocate(l, .before = symAlt)

Hmat <- bind_rows(
  Hbase %>% mutate(sProp = sProp1, propScenario = 1),
  Hbase %>% mutate(sProp = sProp2, propScenario = 2),
  Hbase %>% mutate(sProp = sProp3, propScenario = 3)
) %>%
  dplyr::select(Var1, Var2, Var3, sProp, propScenario)

mu1 <- c(
  Hmat$Var1 * 2,
  Hmat$Var1 * 3,
  Hmat$Var1 * 4,
  Hmat$Var1 * 5,
  Hmat$Var1 * 6
)

mu2 <- c(
  Hmat$Var2 * 2,
  Hmat$Var2 * 3,
  Hmat$Var2 * 4,
  Hmat$Var2 * 5,
  Hmat$Var2 * 6
)

mu3 <- c(
  Hmat$Var3 * 2,
  Hmat$Var3 * 3,
  Hmat$Var3 * 4,
  Hmat$Var3 * 5,
  Hmat$Var3 * 6
)

Hmat <- Hmat %>%
  bind_rows(Hmat, Hmat, Hmat, Hmat) %>%
  bind_cols(mu1, mu2, mu3) %>%
  rename(
    mu1 = "...6",
    mu2 = "...7",
    mu3 = "...8"
  ) %>%
  mutate(number = as.integer(nSNPs * sProp))

Hmat <- Hmat[(case - nConfig + 1):case, ]

## =========================
## define truth and causal
## replicable = same nonzero direction in both studies
## =========================
Hmat <- Hmat %>%
  mutate(
    n_nonzero = (Var1 != 0) + (Var2 != 0) + (Var3 != 0),
    n_zero    = (Var1 == 0) + (Var2 == 0) + (Var3 == 0),
    same_dir  = Var1 == Var2 & Var2 == Var3,
    
    truth = case_when(
      n_nonzero == 0 ~ "null",
      
      n_nonzero == 3 & same_dir ~ "replicable",
      
      n_nonzero == 3 & !same_dir ~ "discordant",
      
      n_nonzero %in% c(1, 2) ~ "partial",
      
      TRUE ~ NA_character_
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
    mu = c(Hmat$mu1[i], Hmat$mu2[i], Hmat$mu3[i]),
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

colnames(allZ) <- c("Z1", "Z2", "Z3")

## SNP-level truth
truthVec  <- rep(Hmat$truth, Hmat$number)
causalVec <- rep(Hmat$causal, Hmat$number)

## p-values
allP <- 1 - pchisq(as.matrix(allZ)^2, df = 1)

totOut <- data.frame(
  Z1 = allZ[, 1],
  Z2 = allZ[, 2],
  Z3 = allZ[, 3],
  P1 = allP[, 1],
  P2 = allP[, 2],
  P3 = allP[, 3],
  origIdx = seq_len(nrow(allP)),
  truth = truthVec,
  causal = causalVec
)

## =========================
## create beta-hat / SE for MAMBA only if requested
## approximate summary-level parameterization
## =========================
if (runMamba) {
  if (!requireNamespace("mamba", quietly = TRUE)) {
    stop("runMamba is TRUE, but the mamba package is not installed.")
  }
  
  summaryStatDir <- here::here("Data")
  real_dat  <- fread(here::here(summaryStatDir, "lc_overall_three.txt"))
  se_dat <- real_dat %>%
    dplyr::select(SE1_real = SE,
                  SE2_real = SE_ukb,
                  SE3_real = SE_mvp) %>%
    dplyr::filter(is.finite(SE1_real), is.finite(SE2_real), is.finite(SE3_real),
                  SE1_real > 0, SE2_real > 0, SE3_real > 0)
  rm(real_dat)
  
  set.seed(aID * 10^3 + 99)
  idx_se <- sample(seq_len(nrow(se_dat)), size = nrow(allZ), replace = TRUE)
  
  totOut <- totOut %>%
    mutate(
      SE1 = se_dat$SE1_real[idx_se],
      SE2 = se_dat$SE2_real[idx_se],
      SE3 = se_dat$SE3_real[idx_se],
      B1 = Z1 * SE1,
      B2 = Z2 * SE2,
      B3 = Z3 * SE3
    )
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



## =========================
## 1. Threshold methods
## =========================
same_dir_idx <- (totOut$Z1 > 0 & totOut$Z2 > 0 & totOut$Z3 > 0) | (totOut$Z1 < 0 & totOut$Z2 < 0 & totOut$Z3 < 0 )

totOut <- totOut %>%
  mutate(
    threshold1 = as.integer(P1 < 1e-8 & P2 < 1e-8 & P3 < 1e-8 & same_dir_idx),
    threshold2 = as.integer(P1 < 1e-6 & P2 < 1e-6 & P3 < 1e-6 & same_dir_idx),
    threshold3 = as.integer(P1 < 1e-5 & P2 < 1e-5 & P3 < 1e-5 & same_dir_idx)
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
zmat <- as.matrix(totOut[, c("Z1", "Z2","Z3")])

initPiList <- list(c(0.82))
for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
initPiList[[8]] <- c(0.1)
tempH <- expand.grid(c(0, 1), c(0, 1), c(0, 1)) %>%
  mutate(s = Var1 + Var2 + Var3) %>%
  arrange(s) %>%
  dplyr::select(-s) %>%
  as.matrix(.)
initMuList <- list(matrix(data=0, nrow=3, ncol=1))
for (i in 2:7) {
  initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
}
initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, sameDirAlt=TRUE, eps=10^(-5))


totOut <- totOut %>% mutate(newSig = newRes$lfdrResults) %>%
  arrange(newSig) %>%
  mutate(newAvg = cummean(newSig)) %>%
  arrange(origIdx)

powerRes$fdpNew <- length(which(totOut$newAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$newAvg < 0.1))
powerRes$powerNew <- length(which(totOut$newAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
powerRes$nRejNew <- length(which(totOut$newAvg < 0.1))

powerRes$inconNew <- powerRes$inconNew <- length(check_incongruous(zMatrix = zmat, lfdrVec = totOut$newLfdr))
## =========================
## 3. Meta-analysis
## =========================
calc_neff <- function(n_case, n_ctrl) {
  4 / (1 / n_case + 1 / n_ctrl)
}

# Overall lung cancer case/control numbers
N_ILCCO <- calc_neff(29266, 56450)
N_UKB   <- calc_neff(2404, 330018)
N_MVP   <- calc_neff(10398, 62708)

w1 <- sqrt(N_ILCCO)  # ILCCO
w2 <- sqrt(N_UKB)    # UKB
w3 <- sqrt(N_MVP)    # MVP

totOut <- totOut %>%
  mutate(
    Z_meta    = (Z1 * w1 + Z2 * w2 + Z3 * w3) / sqrt(w1^2 + w2^2 + w3^2),
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
zmat <- as.matrix(totOut[, c("Z1", "Z2", "Z3")])

bin_out <- ztobins(
  zmat,
  df = 15,
  plot.diagnostics = FALSE
)

H <- hconfigs(
  n.studies = 3,
  n.association.status = 3
)

## only (-1,-1,-1) and (1,1,1) are considered replicated
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
  betajk <- as.matrix(totOut[, c("B1", "B2", "B3")])
  sjk2   <- as.matrix(totOut[, c("SE1", "SE2", "SE3")])^2
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
    zmat = as.matrix(totOut[, c("Z1", "Z2", "Z3")]),
    lfdrVec = totOut$mamba_lfdr,
    method_name = "MAMBA",
    Snum = Snum,
    aID = aID
  )
  
  powerRes$inconMamba <- nrow(incon_mamba)
}

## =========================
## 6.kernel
## =========================
oldResKernel <- emp_bayes_framework(summary_tab = zmat, sameDirAlt=TRUE, kernel = TRUE, joint=FALSE, ind = TRUE,
                                    dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)

totOut <- totOut %>%
  mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
  arrange(kernelLfdr) %>%
  mutate(kernelAvg = cummean(kernelLfdr)) %>%
  arrange(origIdx) %>%
  mutate(kernelSig = as.integer(kernelAvg <= qvalue))

tmp <- get_metrics(totOut$kernelSig, totOut$causal)

powerRes$fdpKernel   <- tmp$fdp
powerRes$powerKernel <- tmp$power
powerRes$nRejKernel  <- tmp$nRej

powerRes$inconKernel <- length(
  check_incongruous(
    zMatrix = zmat,
    lfdrVec = totOut$kernelLfdr
  )
)
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
  fp_breakdown(totOut$repfdr_sig,  totOut$truth, "repfdr",         Snum, aID),
  fp_breakdown(totOut$kernelSig,   totOut$truth, "kernel",         Snum, aID)
)

if (runMamba) {
  fpRes <- bind_rows(
    fpRes,
    fp_breakdown(totOut$mamba_sig, totOut$truth, "MAMBA", Snum, aID)
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

