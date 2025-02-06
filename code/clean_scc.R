#transfer MVP hg38 to hg19.
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
library(tidyr)


# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
dataDir <- here::here("Data")

#Read UKB
rawSCCUKB <- fread(paste0(dataDir,"/ukb_lc/gwasSCC_maf00001_chr1.txt"))
rawSCCUKB <- rawSCCUKB %>%
  select(chromosome,position,allele,RS,beta ,pvalue) %>%
  mutate(chrpos = paste0(chromosome,":",position))


for (i in 2:22) {
  temp <- fread(paste0(dataDir,"/ukb_lc/gwasSCC_maf00001_chr",i,".txt"))
  temp <- temp %>%
    select(chromosome,position,allele,RS,beta ,pvalue) %>%
    mutate(chrpos = paste0(chromosome,":",position))
    rawSCCUKB <- rbind(rawSCCUKB,temp)
}


#Read ILCCO
rawSCCILCCO <- fread(paste0(dataDir,"/replication_with_lcoverall.txt"))
rawSCCILCCO <- rawSCCILCCO %>%
  select(chrpos,A1,A2,p_LC_ilcco,Zlc_ilcco)

#merge UKB and ILCCO
allZ <- rawSCCILCCO %>%
  inner_join(rawSCCUKB, by = "chrpos") %>%
  rename(beta_ukb = beta,pUKB = pvalue)

#Find the duplicated chrpos and retain the lists that have the same alt and ref alleles.

dup_z <- allZ[duplicated(allZ$chrpos) | duplicated(allZ$chrpos, fromLast = TRUE), ]
dup_z <- dup_z %>%
  mutate(merge_ukb_allele_1 = paste0(A1,',',A2)) %>%
  mutate(merge_ukb_allele_2 = paste0(A2,',',A1)) 

uniq_z = dup_z[dup_z$allele == dup_z$merge_ukb_allele_1 | dup_z$allele == dup_z$merge_ukb_allele_2, ]
dup_z = dup_z[!(dup_z$allele == dup_z$merge_ukb_allele_1 | dup_z$allele == dup_z$merge_ukb_allele_2), ]
#Create a list of entries that should be removed.
rm_z = dup_z[,-c("merge_ukb_allele_1","merge_ukb_allele_2")]

#remove them
allZ <- anti_join(allZ, rm_z)

#check the direction of effect (ref allele) and calculate the Z for MVP
allZ <- allZ %>%
  separate(allele, into = c("A1_ukb", "A2_ukb"), sep = ",") %>%
  mutate(same_ukb = ifelse(A1 == A1_ukb & A2 == A2_ukb ,1,0)) %>%
  mutate(flipped_ukb = ifelse(!(A1 == A1_ukb) ,1,0)) %>%
  mutate(Zukb = qnorm(1 - pUKB / 2)) %>%
  mutate(Zukb = ifelse(flipped_ukb == 1  , -Zukb*sign(beta_ukb) , Zukb))  


#read and clean MVP data
rawMVP <- fread(paste0(dataDir,"/lung_cancer_summary_statistics/MVP/LUSC-MVP-EA.tsv.gz"))

cleanMVP <- rawMVP %>%
  select(SNP_ID, Chromosome, Position,`Allele 1`,`Allele 2`, `P value`,`Estimate effect`) %>%
  mutate(chrpos = paste0(Chromosome,":",Position)) %>%
  rename("A1mvp" =`Allele 1`,"A2mvp" =`Allele 2`, "Beta_mvp" = `Estimate effect`, "Pmvp" = `P value`) %>%
  select(chrpos, Pmvp, A1mvp,  A2mvp, SNP_ID, Beta_mvp)



#merge all of them
final_z <- allZ %>%
  inner_join(cleanMVP,by="chrpos")

#Find the duplicated chrpos and retain the lists that have the same alt and ref alleles.
dup_z <- final_z[duplicated(final_z$chrpos) | duplicated(final_z$chrpos, fromLast = TRUE), ]
dup_z <- dup_z %>%
  mutate(merge_ukb_allele_1 = paste0(A1,A2)) %>%
  mutate(merge_ukb_allele_2 = paste0(A2,A1)) %>%
  mutate(merge_mvp = paste0(A1mvp,A2mvp)) 

uniq_z = dup_z[dup_z$merge_mvp == dup_z$merge_ukb_allele_1 | dup_z$merge_mvp == dup_z$merge_ukb_allele_2, ]
dup_z = dup_z[!(dup_z$merge_mvp == dup_z$merge_ukb_allele_1 | dup_z$merge_mvp == dup_z$merge_ukb_allele_2), ]
rm_z = dup_z[, -which(names(dup_z) %in% c("merge_ukb_allele_1", "merge_ukb_allele_2", "merge_mvp"))]


final <- anti_join(final_z, rm_z)

#check the direction of effect (ref allele) and calculate the Z for MVP
final <- final %>%
  mutate(same_mvp = ifelse(A1mvp == A1 & A2mvp == A2 ,1,0)) %>%
  mutate(flipped_mvp = ifelse(!(A1mvp == A1) ,1,0)) %>%
  mutate(Zmvp = qnorm(1 - Pmvp / 2)) %>%
  mutate(Zmvp = ifelse(flipped_mvp == 1  , -Zmvp*sign(Beta_mvp) , Zmvp))  

write.table(final, paste0(dataDir, "/lc_scc_three.txt"), append=F, quote=F, row.names=F, col.names=T)
