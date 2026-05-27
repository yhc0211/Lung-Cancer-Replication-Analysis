# Analysis of real data - LUAD three-way analysis (kernel & repfdr)

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/three_repfdr_adeno.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)
library(locfdr)
library(repfdr)

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
fnameRoot <- paste0(outputDir, "/three_rep_adeno_data_aID", aID)

# where is the data
summaryStatDir <- here::here("Data")

# controls convergence of EM algorithms
oldEps <- 0.01
newEps <- 10^(-5)

#read and clean three gwas summary statistics
cleanUKB <- fread(here::here(summaryStatDir, "lc_adeno_three.txt"))
testDat <- cleanUKB %>% select(zILCCO, zUKB, zMVP)
testDat <- testDat %>% as.matrix(.)



# oldResKernel <- emp_bayes_framework(summary_tab = testDat[,1:3], sameDirAlt=TRUE, kernel = TRUE, joint=FALSE, ind = TRUE,
#                                     dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
# 
# 
# # save
# write.table(oldResKernel$lfdrVec, paste0(fnameRoot, "_newResKernel.txt"), append=F, quote=F, row.names=F, col.names=T)
# 



## =========================
## 4. repfdr
## =========================


bin_out <- ztobins(
  testDat[,1:3],
  df = 7,
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

# save
write.table(rep_out$mat[, "Fdr"], paste0(fnameRoot, "_newResRepfdr.txt"), append=F, quote=F, row.names=F, col.names=T)

