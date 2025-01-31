# Analysis of real data - pleiotropy and replication studies

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/replication_lung.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)
library(locfdr)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
fnameRoot <- paste0(outputDir, "/rep_data_aID", aID)

# where is the data
summaryStatDir <- here::here("Data")

# controls convergence of EM algorithms
oldEps <- 0.01
newEps <- 10^(-5)

# controls which analysis to perform
replication <- FALSE
threeway <- FALSE
cor <- FALSE
if (aID == 1) {
  # CAD, BMI
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  testDat <- cleanUKB %>% select(Zcad, Zbmi, p_CAD, pBMI)
} else if (aID == 2) {
  # Overall LC, CAD
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  testDat <- cleanUKB %>% select(Zoverall, Zcad, pOverall, p_CAD)
} else if (aID == 3) {
  # Overall LC, BMI
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  testDat <- cleanUKB %>% select(Zoverall, Zbmi, pOverall, pBMI)
} else if (aID == 4) {
  # Replication - overall LC, UKB LC
  replication <- TRUE 
  cleanUKB <- fread(here::here(summaryStatDir, "replication_with_lcoverall.txt"))
  testDat <- cleanUKB %>% select(Zoverall, Zlcukb, pOverall, pLCukb)
} else if (aID == 5) {
  # Replication - CAD, UKB CAD
  replication <- TRUE
  cleanUKB <- fread(here::here(summaryStatDir, "cad_for_replication.txt"))
  testDat <- cleanUKB %>% select(Zcad_cardio, Zcadukb, p_CAD_cardio, pCADukb)
} else if (aID == 6) {
  # three way - overall, cad, bmi
  threeway <- TRUE
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  testDat <- cleanUKB %>% select(Zoverall, Zcad, Zbmi, pOverall, p_CAD, pBMI)
} else if (aID == 7) {
  # Correlation - LC UKB, CAD UKB 
  cor <- TRUE 
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  # pvalue doesn't matter 
  testDat <- cleanUKB %>% select(Zlcukb, Zcadukb, pLC, p_CAD)
} else if (aID == 8) {
  # Correlation - LC UKB, bmi
  cor <- TRUE 
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  # pvalue doesn't matter
  testDat <- cleanUKB %>% select(Zlcukb, Zbmi, pLC, pBMI)
} else if (aID == 9) {
  # Correlation - CAD UKB, bmi
  cor <- TRUE 
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  # pvalue doesn't matter
  testDat <- cleanUKB %>% select(Zcadukb, Zbmi, pCAD, pBMI)
}
testDat <- testDat %>% as.matrix(.)

# adjust so test statistics are not too large for R
testDat <- testDat %>% as.matrix(.)
for (col_it in 1:(ncol(testDat)/2)) {
  tooBig <- which(testDat[, col_it] > 8.1)
  tooSmall <- which(testDat[, col_it] < -8.1)
  if (length(tooBig) > 0) {
    testDat[tooBig, col_it] <- 8.1
  }
  if (length(tooSmall) > 0) {
    testDat[tooSmall, col_it] <- -8.1
  }
}
# adjust so p-values are not 0 or 1
for (col_it in ((ncol(testDat)/2)+1):ncol(testDat)) {
  tooBig <- which(testDat[, col_it] == 1)
  tooSmall <- which(testDat[, col_it] == 0)
  minVal <- min(testDat[which(testDat[, col_it] > 0), col_it]) 
  if (length(tooBig) > 0) {                
    testDat[tooBig, col_it] <- 0.999        
  }        
  if (length(tooSmall) > 0) {                
    testDat[tooSmall, col_it] <- minVal          
  }
}

# run HDMT and DACT first for non-three-way analyses
if (!replication & !threeway & !cor) {
  # sped up version of HDMT - produces exactly same output as original file but faster for use on large datasets
  nullprop <- tryCatch(null_estimation(testDat[, 3:4]), error=function(e) e, warning=function(w) w)
  if (class(nullprop)[1] == "list") {
    hdmtOut <- tryCatch(fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,testDat[,3:4]),
                        error = function(e) e, warning = function(w) w)
  } else {hdmtOut <- rep(NA, nrow(testDat))}
  if (class(nullprop)[1] != "list" | class(hdmtOut)[1] != "data.frame") {
    hdmtOut <- rep(NA, nrow(testDat))
  } 
  # save
  write.table(hdmtOut, paste0(fnameRoot, "_hdmt.txt"), append=F, quote=F, row.names=F, col.names=T)
  
  # sped up version of DACT - negligible differences but much faster
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  DACTout <- tryCatch(DACT(p_a = testDat[, 3], p_b = testDat[, 4], nullEst=nullprop, correction="JC"), 
                      error = function(e) e, warning=function(w) w)
  if (class(DACTout)[1] != "list") {
    DACTfreqp <- rep(NA, nrow(testDat))
  } else {
    DACTfreqp <- DACTout$p_dact
  }
  
  # save
  write.table(DACTfreqp, paste0(fnameRoot, "_DACTp.txt"), append=F, quote=F, row.names=F, col.names=T)
  
} 

# run kernel
if (!cor) {
  oldResKernel <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = TRUE, joint=FALSE, ind = TRUE,
                                      dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldResKernel)[1] != "list") {
    kernelLfdr <- rep(NA, nrow(testDat))
  } else {
    kernelLfdr <- oldResKernel$lfdrVec
  }
  # save
  write.table(kernelLfdr, paste0(fnameRoot, "_kernel.txt"), append=F, quote=F, row.names=F, col.names=T)
}

# run 7 df
if (!cor) {
  oldRes7df <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = FALSE, joint=FALSE, ind = TRUE,
                                   dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes7df)[1] != "list") {
    df7Lfdr <- rep(NA, nrow(testDat))
  } else {
    df7Lfdr <- oldRes7df$lfdrVec
  }
  # save
  write.table(df7Lfdr, paste0(fnameRoot, "_df7.txt"), append=F, quote=F, row.names=F, col.names=T)
}

# run 50 df
if (!cor) {
  oldRes50df <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = FALSE, joint=FALSE, ind = TRUE,
                                    dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes50df)[1] != "list") {
    df50Lfdr <- rep(NA, nrow(testDat))
  } else {
    df50Lfdr <- oldRes50df$lfdrVec
  }
  # save
  write.table(df50Lfdr, paste0(fnameRoot, "_df50.txt"), append=F, quote=F, row.names=F, col.names=T)
}

# 2D cases pleiotropy
if (!threeway) {
  initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
  initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                     matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
  if (!cor) { 
    newRes <- symm_fit_ind_EM(testStats = testDat[, 1:2], sameDirAlt = replication, initMuList = initMuList, initPiList = initPiList, eps=newEps)
  } else {
    tempCorVal <- cor(testDat[, 1], testDat[, 2])
    tempCorMat <- matrix(data=c(1, tempCorVal, tempCorVal, 1), nrow=2)
    newRes <- symm_fit_cor_EM(testStats = testDat[, 1:2], corMat = tempCorMat, initMuList = initMuList, initPiList = initPiList, eps=newEps)
  }
} else {
  # 3D independent 
  initPiList <- list(c(0.82))
  for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
  initPiList[[8]] <- c(0.1)
  tempH <- expand.grid(c(0, 1), c(0, 1), c(0, 1)) %>%
    mutate(s = Var1 + Var2 + Var3) %>%
    arrange(s) %>%
    select(-s) %>%
    as.matrix(.)
  initMuList <- list(matrix(data=0, nrow=3, ncol=1))
  for (i in 2:7) {
    initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
  }
  initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
  
  newRes <- symm_fit_ind_EM(testStats = testDat[, 1:3], initMuList = initMuList, initPiList = initPiList, eps = newEps)
  
}

# save
write.table(newRes$lfdrResults, paste0(fnameRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$muInfo), paste0(fnameRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$piInfo), paste0(fnameRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)





