# Analysis of real data - three-way adeno replication studies

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/three_rep_adeno.R")

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
fnameRoot <- paste0(outputDir, "/rep_three_adeno_aID", aID)

# where is the data
summaryStatDir <- here::here("Data")

# controls convergence of EM algorithms
oldEps <- 0.01
newEps <- 10^(-5)

# controls which analysis to perform

cleanUKB <- fread(here::here(summaryStatDir, "/lc_adeno_three.txt"))
testDat <- cleanUKB %>% select(zILCCO, zUKB, zMVP)

testDat <- testDat %>% as.matrix(.)

# adjust so test statistics are not too large for R
# testDat <- testDat %>% as.matrix(.)
# for (col_it in 1:ncol(testDat)) {
#   tooBig <- which(testDat[, col_it] > 8.1)
#   tooSmall <- which(testDat[, col_it] < -8.1)
#   if (length(tooBig) > 0) {
#     testDat[tooBig, col_it] <- 8.1
#   }
#   if (length(tooSmall) > 0) {
#     testDat[tooSmall, col_it] <- -8.1
#   }
#}



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

newRes <- symm_fit_ind_EM(testStats = testDat[, 1:3], sameDirAlt = TRUE,initMuList = initMuList, initPiList = initPiList, eps = newEps)



# save
write.table(newRes$lfdrResults, paste0(fnameRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$muInfo), paste0(fnameRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$piInfo), paste0(fnameRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)

