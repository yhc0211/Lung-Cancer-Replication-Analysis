# Analysis of two-cohort overall lung cancer real data - incongruous ranking between z-score ranking and lfdr ranking (kernel & repfdr)

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/two_incon_overall.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)
library(repfdr)
library(locfdr)


# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

## =========================
## command line args
## =========================
args <- commandArgs(trailingOnly = TRUE)
aID  <- as.numeric(args[1])
Snum <- as.numeric(args[2])


# set output directory 
outputDir <- here::here("Lung", "output")
fnameRoot <- paste0(outputDir, "/two_rep_overall_aID", aID)

# where is the data
summaryStatDir <- here::here("Data")



#read and clean three gwas summary statistics
cleanUKB <- fread(here::here(summaryStatDir, "lc_overall_three.txt"))
testDat <- cleanUKB %>% select(zILCCO,zUKB,zMVP)



library(dplyr)
library(data.table)

run_incon_for_aID <- function(aID, testDat, fnameRoot) {
  
  # Select GWAS summary statistic pair
  if (aID == 1) {
    pair_name <- "ILCCO_UKB"
    zMat <- testDat %>% select(zILCCO, zUKB) %>% as.matrix()
  } else if (aID == 2) {
    pair_name <- "ILCCO_MVP"
    zMat <- testDat %>% select(zILCCO, zMVP) %>% as.matrix()
  } else if (aID == 3) {
    pair_name <- "UKB_MVP"
    zMat <- testDat %>% select(zUKB, zMVP) %>% as.matrix()
  } else {
    stop("aID must be 1, 2, or 3")
  }
  
  # Read lfdr results
  lfdrKernel <- fread(paste0(fnameRoot, "_newResKernel.txt"))
  lfdrRepfdr <- fread(paste0(fnameRoot, "_newResRepfdr.txt"))
  lfdrnew    <- fread(paste0(fnameRoot, "_newlfdr.txt"))
  
  # Calculate incongruous numbers
  num_incon_Kernel <- length(
    check_incongruous(
      zMatrix = zMat[, 1:2],
      lfdrVec = lfdrKernel$x
    )
  )
  
  num_incon_Repfdr <- length(
    check_incongruous(
      zMatrix = zMat[, 1:2],
      lfdrVec = lfdrRepfdr$x
    )
  )
  
  num_incon_new <- length(
    check_incongruous(
      zMatrix = zMat[, 1:2],
      lfdrVec = lfdrnew$x
    )
  )
  
  data.frame(
    aID = aID,
    pair = pair_name,
    num_incon_Kernel = num_incon_Kernel,
    num_incon_Repfdr = num_incon_Repfdr,
    num_incon_new = num_incon_new
  )
}

res_all <- bind_rows(
  lapply(1:3, function(id) {
    run_incon_for_aID(
      aID = id,
      testDat = testDat,
      fnameRoot = fnameRoot
    )
  })
)

write.table(
  res_all,
  paste0(fnameRoot, "_two_incon_overall.txt"),
  append = FALSE,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)