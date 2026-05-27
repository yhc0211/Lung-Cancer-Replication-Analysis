# Analysis of three-cohort overall lung cancer real data - incongruous ranking between z-score ranking and lfdr ranking (kernel & repfdr)

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/three_incon_overall.R")

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
fnameRoot <- paste0(outputDir, "/three_rep_data_aID", aID)
# where is the data
summaryStatDir <- here::here("Data")



#read and clean three gwas summary statistics
cleanUKB <- fread(here::here(summaryStatDir, "lc_overall_three.txt"))
testDat <- cleanUKB %>% select(zILCCO,zUKB,zMVP) %>% as.matrix()

  
# Read lfdr results
lfdrKernel <- fread(paste0(fnameRoot, "_newResKernel.txt"))
lfdrRepfdr <- fread(paste0(fnameRoot, "_newResRepfdr.txt"))
lfdrnew    <- fread(paste0(fnameRoot, "_newlfdr.txt"))

# Calculate incongruous numbers
num_incon_Kernel <- length(
  check_incongruous(
    zMatrix = testDat[, 1:3],
    lfdrVec = lfdrKernel$x
  )
)

num_incon_Repfdr <- length(
  check_incongruous(
    zMatrix = testDat[, 1:3],
    lfdrVec = lfdrRepfdr$x
  )
)

num_incon_new <- length(
  check_incongruous(
    zMatrix = testDat[, 1:3],
    lfdrVec = lfdrnew$x
  )
)

res_all <- data.frame(
  num_incon_Kernel = num_incon_Kernel,
  num_incon_Repfdr = num_incon_Repfdr,
  num_incon_new = num_incon_new
)




write.table(
  res_all,
  paste0(fnameRoot, "_three_incon_overall.txt"),
  append = FALSE,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)