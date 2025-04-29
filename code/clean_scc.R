#clean scc
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/clean_scc.R")

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
rawUKB <- fread(paste0(dataDir,"/lung_cancer_summary_statistics/UKB/full_sccLung_maf00005_regenie.txt"))
rawUKB <- rawUKB %>% 
  select(CHROM,GENPOS,ALLELE0,ALLELE1,BETA,SE,LOG10P ) %>%
  na.omit() %>%
  mutate(chrpos = paste0(CHROM,":",GENPOS)) %>%
  mutate(pUKB = 10^-(LOG10P)) %>%
  rename(chrom = CHROM,pos = GENPOS,A1_UKB = ALLELE1, A2_UKB = ALLELE0 ,Beta_ukb = BETA, SE_ukb =SE)


#Read ILCCO 
rawILCCO <- fread(paste0(dataDir,"/lung_cancer_summary_statistics/ILCCO/LC_Squam_hg19_chr",1,".txt")) 


for (i in 2:22) {
  temp <- fread(paste0(dataDir,"/lung_cancer_summary_statistics/ILCCO/LC_Squam_hg19_chr",i,".txt"))
  rawILCCO <- rbind(rawILCCO,temp)
}

rawILCCO <- rawILCCO %>%
  select("chr:pos","RS","reference_allele","effect_allele","OR_fixed","StdError_fixed","P-value") %>%
  rename("chrpos" = "chr:pos","A1_ILCCO" = "effect_allele", "A2_ILCCO" = "reference_allele", "OR" = "OR_fixed","SE" = "StdError_fixed" ,"pILCCO" = "P-value" ) %>%
  mutate(Beta_ILCCO = log(OR)) %>%
  mutate(zILCCO = Beta_ILCCO/SE)


#merge ILCCO and UKB
allZ <- rawILCCO %>%
  inner_join(rawUKB, by = "chrpos") 


#Find the duplicated chrpos and retain the lists that have the same alt and ref alleles.
dup_z <- allZ[duplicated(allZ$chrpos) | duplicated(allZ$chrpos, fromLast = TRUE), ]

dup_z <- dup_z %>%
  mutate(merge_ilcco_allele = paste0(A1_ILCCO,',',A2_ILCCO)) %>%
  mutate(merge_ukb_allele_1 = paste0(A1_UKB,',',A2_UKB)) %>%
  mutate(merge_ukb_allele_2 = paste0(A2_UKB,',',A1_UKB)) 

uniq_z = dup_z[dup_z$merge_ilcco_allele == dup_z$merge_ukb_allele_1 | dup_z$merge_ilcco_allele == dup_z$merge_ukb_allele_2, ]
dup_z = dup_z[!(dup_z$merge_ilcco_allele == dup_z$merge_ukb_allele_1 | dup_z$merge_ilcco_allele == dup_z$merge_ukb_allele_2), ]

#Create a list of entries that should be removed.
rm_z = dup_z[,-c("merge_ilcco_allele","merge_ukb_allele_1","merge_ukb_allele_2")]


#remove them
allZ <- anti_join(allZ, rm_z)  


#check the direction of effect (ref allele) and calculate the Z for UKB
allZ <- allZ %>%
  mutate(same_ukb = ifelse(A1_ILCCO == A1_UKB & A2_ILCCO == A2_UKB ,1,0)) %>%
  mutate(flipped_ukb = ifelse(!(A1_ILCCO == A1_UKB) ,1,0)) %>%
  mutate(zUKB = Beta_ukb/SE_ukb) %>%
  mutate(zUKB = ifelse(flipped_ukb == 1  , -zUKB , zUKB))  



#read mvp
rawMVP <- fread(paste0(dataDir,"/lung_cancer_summary_statistics/MVP/LUSC-MVP-EA.tsv.gz"))
cleanMVP <- rawMVP %>%
  select(SNP_ID, Chromosome, Position,`Allele 1`,`Allele 2`, `P value`,`Estimate effect`,SE) %>%
  mutate(chrpos = paste0(Chromosome,":",Position)) %>%
  rename("A1_mvp" =`Allele 1`,"A2_mvp" =`Allele 2`, "Beta_mvp" = `Estimate effect`, "SE_mvp" = SE,"pMVP" = `P value`) %>%
  select(chrpos, pMVP, A1_mvp,  A2_mvp, SNP_ID, Beta_mvp,SE_mvp)



#merge all of them
final_z <- allZ %>%
  inner_join(cleanMVP,by="chrpos")


#Find the duplicated chrpos and retain the lists that have the same alt and ref alleles.
dup_z <- final_z[duplicated(final_z$chrpos) | duplicated(final_z$chrpos, fromLast = TRUE), ]
dup_z <- dup_z %>%
  mutate(merge_ilcco_allele = paste0(A1_ILCCO,A2_ILCCO)) %>%
  mutate(merge_mvp_allele_1 = paste0(A1_mvp,A2_mvp)) %>%
  mutate(merge_mvp_allele_2 = paste0(A2_mvp,A1_mvp)) 

uniq_z = dup_z[dup_z$merge_ilcco_allele == dup_z$merge_mvp_allele_1 | dup_z$merge_ilcco_allele == dup_z$merge_mvp_allele_2, ]
dup_z = dup_z[!(dup_z$merge_ilcco_allele == dup_z$merge_mvp_allele_1 | dup_z$merge_ilcco_allele == dup_z$merge_mvp_allele_2), ]

#Create a list of entries that should be removed.
rm_z = dup_z[,-c("merge_ilcco_allele","merge_mvp_allele_1","merge_mvp_allele_2")]


#remove them
final_z <- anti_join(final_z, rm_z)  



#check the direction of effect (ref allele) and calculate the Z for MVP
final <- final_z %>%
  mutate(same_mvp = ifelse(A1_mvp == A1_ILCCO & A2_mvp == A2_ILCCO ,1,0)) %>%
  mutate(flipped_mvp = ifelse(!(A1_mvp == A1_ILCCO) ,1,0)) %>%
  mutate(zMVP = Beta_mvp/SE_mvp) %>%
  mutate(zMVP = ifelse(flipped_mvp == 1  , -zMVP , zMVP)) 

write.table(final, paste0(dataDir, "/lc_scc_three.txt"), append=F, quote=F, row.names=F, col.names=T)
