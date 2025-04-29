# Make Figure 4: Manhattan Plots for Three GWAS Summary Statistics
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/Fig4.R")

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
library(ggpubr)
library(tidyr)


# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
outputDir_Fig <- here::here("Lung", "Fig_Tab")
dataDir <- here::here("Data")



fdrLimitNew <- 0.1

#replication analysis
three_Z <- fread(here::here(dataDir, "lc_overall_three.txt")) 


tempNew <- fread(paste0(outputDir ,"/three_rep_data_aID1_newlfdr.txt"), header=T, data.table=F)
tempDat <- three_Z %>% select(chrpos,zILCCO, zUKB, zMVP) %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)


#meta analysis
allZ_meta <- three_Z %>%  
  rename(P1 = pILCCO) %>%
  rename(Z1 = zILCCO) %>%
  rename(P2 = pUKB) %>%
  rename(Z2 = zUKB) %>%
  rename(P3 = pMVP) %>%
  rename(Z3 = zMVP) %>%
  select(P1,Z1,P2,Z2,P3,Z3,chrpos)


w1 = sqrt(85716) #overall
w2 = sqrt(332422) #ukb
w3 = sqrt(73106) #mvp

allZ_meta$tempZ1 = qnorm(allZ_meta$P1/ 2) * sign(allZ_meta$Z1)
allZ_meta$tempZ2 = qnorm(allZ_meta$P2/ 2) * sign(allZ_meta$Z2)
allZ_meta$tempZ3 = qnorm(allZ_meta$P3/ 2) * sign(allZ_meta$Z3)

allZ_meta$Z_meta = (allZ_meta$tempZ1*w1 + allZ_meta$tempZ2*w2 + allZ_meta$tempZ3*w3 )/sqrt(w1^2 + w2^2 + w3^2)
allZ_meta$P_meta =2*pnorm(abs(-allZ_meta$Z_meta), lower.tail=FALSE)
allZ_meta_sig <- allZ_meta %>%
  filter(P_meta < 1e-8)


###ILCCO
manDatILCCO <- three_Z %>%
  select("chrpos","chrom","pos","pILCCO") %>%
  filter(pILCCO < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pILCCO, method="BH"))) %>%
  mutate(chrom= paste0("chr", chrom)) 


manDatILCCO$cat <- forcats::fct_relevel(manDatILCCO$cat, c("Both" ,
                                                           "Replication Analysis",
                                                           "Meta Analysis",
                                                           "Neither"))



plotmanhattan = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  names(x_breaks)[c(9, 11, 13, 14,16, 17, 18,20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(FDR)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10),legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}
manPlotILCCO_overall <- plotmanhattan(manDatILCCO,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,3,5)], shapeValues=c(16,17,18,19),
                                      ylimits=c(3.9, 100), legName ="Method",title = "ILCCO" )


###UKB
manDatUKB<- three_Z %>%
  select("chrpos","chrom","pos","pUKB") %>%
  filter(pUKB < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pUKB, method="BH"))) %>%
  mutate(chrom= paste0("chr",chrom)) 


manDatUKB$cat <- forcats::fct_relevel(manDatUKB$cat, c("Both" ,
                                                       "Replication Analysis",
                                                       "Meta Analysis",
                                                       "Neither"))



manPlotUKB_overall <- plotmanhattan(manDatUKB,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,3,5)], shapeValues=c(16,17,18,19),
                                    ylimits=c(3.9, 6), legName ="Method",title = "UKB" )

###MVP
manDatMVP<- three_Z %>%
  select("chrpos","chrom","pos","pMVP") %>%
  filter(pMVP < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pMVP, method="BH"))) %>%
  mutate(chrom= paste0("chr",chrom)) 

manDatMVP$cat <- forcats::fct_relevel(manDatMVP$cat, c("Both" ,
                                                       "Replication Analysis",
                                                       "Meta Analysis",
                                                       "Neither"))

manPlotMVP_overall <- plotmanhattan(manDatMVP,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,3,5)], shapeValues=c(16,17,18,19),
                                    ylimits=c(3.9, 30), legName ="Method",title = "MVP" )






####### SCC #########


#replication analysis
three_Z <- fread(here::here(dataDir, "lc_scc_three.txt"))


tempNew <- fread(paste0(outputDir ,"/rep_three_scc_aID1_newlfdr.txt"), header=T, data.table=F)
tempDat <- three_Z %>% select("chrpos","zILCCO" ,"zUKB"  ,"zMVP") %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)

#meta analysis
allZ_meta <- three_Z %>%  
  rename(P1 = pILCCO) %>%
  rename(Z1 = zILCCO) %>%
  rename(P2 = pUKB) %>%
  rename(Z2 = zUKB) %>%
  rename(P3 = pMVP) %>%
  rename(Z3 = zMVP) %>%
  select(P1,Z1,P2,Z2,P3,Z3,chrpos)


w1 = sqrt(63053) #ilcco
w2 = sqrt(332422) #ukb
w3 = sqrt(10325) #mvp

allZ_meta$tempZ1 = qnorm(allZ_meta$P1/ 2) * sign(allZ_meta$Z1)
allZ_meta$tempZ2 = qnorm(allZ_meta$P2/ 2) * sign(allZ_meta$Z2)
allZ_meta$tempZ3 = qnorm(allZ_meta$P3/ 2) * sign(allZ_meta$Z3)

allZ_meta$Z_meta = (allZ_meta$tempZ1*w1 + allZ_meta$tempZ2*w2 + allZ_meta$tempZ3*w3 )/sqrt(w1^2 + w2^2 + w3^2)
allZ_meta$P_meta =2*pnorm(abs(-allZ_meta$Z_meta), lower.tail=FALSE)
allZ_meta_sig <- allZ_meta %>%
  filter(P_meta < 1e-8)


###ILCCO
manDatILCCO_scc <- three_Z %>%
  select("chrpos","chrom","pos","pILCCO") %>%
  filter(pILCCO < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pILCCO, method="BH"))) %>%
  mutate(chrom= paste0("chr", chrom)) 



manDatILCCO_scc$cat <- forcats::fct_relevel(manDatILCCO_scc$cat, c("Both" ,
                                                           "Meta Analysis",
                                                           "Replication Analysis",
                                                           "Neither"))



plotmanhattan_scc = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  names(x_breaks)[c(9, 11, 13,14, 16, 17, 19,20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(FDR)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10),legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}
manPlotILCCO_scc <- plotmanhattan_scc(manDatILCCO_scc,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,5,3)], shapeValues=c(16,17,15,18),
                                  ylimits=c(3.9, 40), legName ="Method",title = "ILCCO" )


###UKB
manDatUKB_scc<- three_Z %>%
  select("chrpos","chrom","pos","pUKB") %>%
  filter(pUKB < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pUKB, method="BH"))) %>%
  mutate(chrom= paste0("chr",chrom)) 

manDatUKB_scc$cat <- forcats::fct_relevel(manDatUKB_scc$cat, c("Both" ,
                                                       "Meta Analysis",
                                                       "Replication Analysis",
                                                       "Neither"))



manPlotUKB_scc <- plotmanhattan_scc(manDatUKB_scc,build = "hg19", colValues = wes_palette("Darjeeling1")[c(3)], shapeValues=c(18),
                                ylimits=c(3.9, 5), legName ="Method",title = "UKB" )

###MVP
manDatMVP_scc<- three_Z %>%
  select("chrpos","chrom","pos","pMVP") %>%
  filter(pMVP < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pMVP, method="BH"))) %>%
  mutate(chrom= paste0("chr",chrom)) 

manDatMVP_scc$cat <- forcats::fct_relevel(manDatMVP_scc$cat, c("Both" ,
                                                       "Meta Analysis",
                                                       "Replication Analysis",
                                                       "Neither"))

manPlotMVP_scc <- plotmanhattan_scc(manDatMVP_scc,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,5,3)], shapeValues=c(16,17,15,18),
                                ylimits=c(3.9, 20), legName ="Method",title = "MVP" )





##########Adeno#############
three_Z <- fread(here::here(dataDir, "lc_adeno_three.txt"))

tempNew <- fread(paste0(outputDir ,"/rep_three_adeno_aID1_newlfdr.txt"), header=T, data.table=F)
tempDat <- three_Z %>% select("chrpos","zILCCO" ,"zUKB"  ,"zMVP") %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)

#meta analysis
allZ_meta <- three_Z %>%  
  rename(P1 = pILCCO) %>%
  rename(Z1 = zILCCO) %>%
  rename(P2 = pUKB) %>%
  rename(Z2 = zUKB) %>%
  rename(P3 = pMVP) %>%
  rename(Z3 = zMVP) %>%
  select(P1,Z1,P2,Z2,P3,Z3,chrpos)


w1 = sqrt(61073) #ilcco
w2 = sqrt(332422) #ukb
w3 = sqrt(14132) #mvp

allZ_meta$tempZ1 = qnorm(allZ_meta$P1/ 2) * sign(allZ_meta$Z1)
allZ_meta$tempZ2 = qnorm(allZ_meta$P2/ 2) * sign(allZ_meta$Z2)
allZ_meta$tempZ3 = qnorm(allZ_meta$P3/ 2) * sign(allZ_meta$Z3)

allZ_meta$Z_meta = (allZ_meta$tempZ1*w1 + allZ_meta$tempZ2*w2 + allZ_meta$tempZ3*w3 )/sqrt(w1^2 + w2^2 + w3^2)
allZ_meta$P_meta =2*pnorm(abs(-allZ_meta$Z_meta), lower.tail=FALSE)
allZ_meta_sig <- allZ_meta %>%
  filter(P_meta < 1e-8)


###ILCCO
manDatILCCO_adeno <- three_Z %>%
  select("chrpos","chrom","pos","pILCCO") %>%
  filter(pILCCO < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pILCCO, method="BH"))) %>%
  mutate(chrom= paste0("chr", chrom)) 


manDatILCCO_adeno$cat <- forcats::fct_relevel(manDatILCCO_adeno$cat, c("Both" ,
                                                           "Meta Analysis",
                                                           "Replication Analysis",
                                                           "Neither"))



plotmanhattan_adeno = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  names(x_breaks)[c(9, 11, 13, 14,16, 17,18, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(FDR)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}
manPlotILCCO_adeno <- plotmanhattan(manDatILCCO_adeno,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,5,3)], shapeValues=c(16,17,15,18),
                                    ylimits=c(3.9, 45), legName ="Method",title = "ILCCO" )


###UKB
manDatUKB_adeno<- three_Z %>%
  select("chrpos","chrom","pos","pUKB") %>%
  filter(pUKB < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pUKB, method="BH"))) %>%
  mutate(chrom= paste0("chr",chrom)) 

manDatUKB_adeno$cat <- forcats::fct_relevel(manDatUKB_adeno$cat, c("Both" ,
                                                       "Meta Analysis",
                                                       "Replication Analysis",
                                                       "Neither"))



manPlotUKB_adeno <- plotmanhattan(manDatUKB_adeno,build = "hg19", colValues =wes_palette("Darjeeling1")[c(1,3)], shapeValues=c(16,18),
                                  ylimits=c(3.9, 6), legName ="Method",title = "UKB" )


###MVP
manDatMVP_adeno<- three_Z %>%
  select("chrpos","chrom","pos","pMVP") %>%
  filter(pMVP < 1e-4) %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(y = -log10(p.adjust(pMVP, method="BH"))) %>%
  mutate(chrom= paste0("chr",chrom)) 

manDatMVP_adeno$cat <- forcats::fct_relevel(manDatMVP_adeno$cat, c("Both" ,
                                                       "Meta Analysis",
                                                       "Replication Analysis",
                                                       "Neither"))

manPlotMVP_adeno<- plotmanhattan(manDatMVP_adeno,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,5,3)], shapeValues=c(16,15,18),
                                 ylimits=c(3.9, 20), legName ="Method",title = "MVP" )


mergePlot1 <- ggarrange(manPlotILCCO_overall, manPlotMVP_overall , manPlotUKB_overall, ncol=3, nrow=1)
mergePlot2 <- ggarrange(manPlotILCCO_scc, manPlotMVP_scc, manPlotUKB_scc, ncol=3, nrow=1)
mergePlot3 <- ggarrange(manPlotILCCO_adeno, manPlotMVP_adeno, manPlotUKB_adeno,ncol=3, nrow=1, common.legend = TRUE, legend="bottom")


.libPaths("/home/ychang11/R/ubuntu/4.3.1")
library(plotgardener)


pageCreate(width = 8, height = 11, default.units = "inches")
plotGG(
  plot = mergePlot1,
  x = 0.1, y = 0.6,
  width = 8, height = 3 , just = c("left", "top")
)
plotGG(
  plot = mergePlot2,
  x = 0.1, y = 4,
  width = 8, height = 3 , just = c("left", "top")
)
plotGG(
  plot = mergePlot3,
  x = 0.1, y = 7.3,
  width = 8, height = 3.1 , just = c("left", "top")
)
plotText(label = "A. Overall", x = 0.5, y = 0.3,
         fontsize = 15, fontface = "bold", just = "center",
         default.units = "inches")
plotText(label = "B. LUSC", x = 0.5, y = 3.7,
         fontsize = 15, fontface = "bold", just = "center",
         default.units = "inches")
plotText(label = "C. LUAD", x = 0.5, y = 7,
         fontsize = 15, fontface = "bold", just = "center",
         default.units = "inches")
pageGuideHide()