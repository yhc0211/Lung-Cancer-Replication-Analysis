#make a three way Z statistics for overall Lung Cancer
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/clean_three_overall.R")

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(data.table)
library(xtable)
library(devtools)
library(wesanderson)
library(here)
library(manhattan)


# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
dataDir <- here::here("Data")

# read MVP data and (ILCCO, UKB) data
rawMVP <- fread(paste0(dataDir,"/lung_cancer_summary_statistics/MVP/LC-MVP-EA.tsv.gz"))
allZ <- fread(here::here(dataDir, "replication_with_lcoverall.txt")) 

#clean MVP data
cleanMVP <- rawMVP %>%
  select(SNP_ID, Chromosome, Position,`Allele 1`,`Allele 2`, `P value`,`Estimate effect`) %>%
  mutate(chrpos = paste0(Chromosome,":",Position)) %>%
  rename("A1mvp" =`Allele 1`,"A2mvp" =`Allele 2`, "Beta_mvp" = `Estimate effect`, "Pmvp" = `P value`) %>%
  select(chrpos, Pmvp, A1mvp,  A2mvp, SNP_ID, Beta_mvp)


#merge MVP with (ILCCO, UKB) data
final_z <- allZ %>%
  inner_join(cleanMVP,by="chrpos")

#Find the duplicated chrpos and retain the lists that have the same alt and ref alleles as those in the UKB and ILCCO datasets.
dup_z <- dup_z %>%
  mutate(merge_overall_1 = paste0(A1overall,A2overall)) %>%
  mutate(merge_overall_2 = paste0(A2overall,A1overall)) %>%
  mutate(merge_mvp = paste0(A1mvp,A2mvp)) 

dup_z <- final_z[duplicated(final_z$chrpos) | duplicated(final_z$chrpos, fromLast = TRUE), ]
dup_z = dup_z[!(dup_z$merge_mvp == dup_z$merge_overall_1 | dup_z$merge_mvp == dup_z$merge_overall_2), ]
uniq_z = dup_z[dup_z$merge_mvp == dup_z$merge_overall_1 | dup_z$merge_mvp == dup_z$merge_overall_2, ]

#Create a list of entries that should be removed.
rm_z = dup_z[,-c("merge_overall_1","merge_overall_2","merge_mvp")]

#remove them
final <- anti_join(final_z, rm_z)

#check the direction of effect (ref allele) and calculate the Z for MVP
final <- final %>%
  mutate(same_mvp = ifelse(A1mvp == A1overall & A2mvp == A2overall ,1,0)) %>%
  mutate(flipped_mvp = ifelse(!(A1mvp == A1overall) ,1,0)) %>%
  mutate(Zmvp = qnorm(1 - Pmvp / 2)) %>%
  mutate(Zmvp = ifelse(flipped_mvp == 1  , -Zmvp*sign(Beta_mvp) , Zmvp))  

write.table(final, paste0(dataDir, "/lc_overall_three.txt"), append=F, quote=F, row.names=F, col.names=T)
