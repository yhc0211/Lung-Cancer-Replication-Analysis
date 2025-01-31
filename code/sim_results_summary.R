here::i_am("Lung/summarize_replication_results.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- 4
Snum <- 1
file_it <- 4

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
dataDir <- here::here("Data")
fnameOut <- paste0(outputDir, "/processed_ukb_data_S", Snum, ".txt")

# how the raw output files are named
fnameRoot <- paste0(outputDir, "/rep_sim_aID", 1:5)
fnameNew <- c(paste0(fnameRoot, "_newlfdr.txt"))

# hold temporary results
tempRes <- data.frame(Method=c("Sim1", "Sim2", "Sim3", "Sim4", "Sim5"), numReject=NA)


# new 
for (file_it in 1:5) {
  tempNew <- fread(fnameNew[file_it], header=T, data.table=F)
  tempDat <- tempNew %>%
    mutate(origIdx = 1:nrow(.)) %>%   
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr)) %>%
    mutate(rejNew = ifelse(cumNew < 0.1, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[file_it] <- sum(tempDat$rejNew)
  
}
