# Make Figure 5: Three-way Replication Analysis Results for ILCCO and MVP

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/Fig5.R")

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


#two-way replication analysis for ILCCO and MVP
#replication analysis
three_Z <- fread(here::here(dataDir, "lc_overall_three.txt"))
aID = 2

#parameter
para <- data.frame(
  p_1 = c("pILCCO","pILCCO","pUKB"),
  p_2 = c("pUKB","pMVP","pMVP"),
  Z_1 = c("zILCCO","zILCCO","zUKB"),
  Z_2 = c("zUKB","zMVP","zMVP"),
  ss_1 = c(85716,85716,332422),
  ss_2 = c(332422,73106,73106))

fdrLimitNew <- 0.1

tempNew <- fread(paste0(outputDir ,"/two_rep_overall_aID",aID,"_newlfdr.txt"), header=T, data.table=F)


tempDat <- three_Z %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)

#meta analysis
allZ_meta <- three_Z %>%  
  rename(P1 = para$p_1[aID]) %>%
  rename(Z1 = para$Z_1[aID]) %>%
  rename(P2 = para$p_2[aID]) %>%
  rename(Z2 = para$Z_2[aID]) %>%
  select(P1,Z1,P2,Z2,chrpos)


w1 = sqrt(para$ss_1[aID])
w2 = sqrt(para$ss_2[aID]) 


allZ_meta$tempZ1 = qnorm(allZ_meta$P1/ 2) * sign(allZ_meta$Z1)
allZ_meta$tempZ2 = qnorm(allZ_meta$P2/ 2) * sign(allZ_meta$Z2)

allZ_meta$Z_meta = (allZ_meta$Z1*w1 + allZ_meta$Z2*w2)/sqrt(w1^2 + w2^2)
allZ_meta$P_meta = 2*pnorm(abs(-allZ_meta$Z_meta), lower.tail=FALSE)
allZ_meta_sig <- allZ_meta %>%
  filter(P_meta < 1e-8) %>%
  mutate(method = "Meta Analysis") %>%
  select(chrpos,method)

#all replication
merge_results <- tempDat %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y,cat) 




merge_results$cat <- forcats::fct_relevel(merge_results$cat, c("Both" ,
                                                               "Meta Analysis",
                                                               "Replication Analysis",
                                                               "Neither"))

plotmanhattan = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
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
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold",size =12))
  
}


manPlot_overall <- plotmanhattan(merge_results,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,3)], shapeValues=c(16,17,18),
                                 ylimits=c(0, 70), legName ="Method",title = "Overall" )

ggsave(paste0(outputDir_Fig,"/two_rep_overall.png"), plot = manPlot_overall, width = 6,
       height = 3)




#########SCC###########
three_Z <- fread(here::here(dataDir, "lc_scc_three.txt"))


para <- data.frame(
  p_1 = c("pILCCO","pILCCO","pUKB"),
  p_2 = c("pUKB","pMVP","pMVP"),
  Z_1 = c("zILCCO","zILCCO","zUKB"),
  Z_2 = c("zUKB","zMVP","zMVP"),
  ss_1 = c(63053,63053,332422),
  ss_2 = c(332422,10325,10325))

fdrLimitNew <- 0.1

aID = 2

tempNew <- fread(paste0(outputDir ,"/two_rep_scc_aID",aID,"_newlfdr.txt"), header=T, data.table=F)


tempDat <- three_Z %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)

#meta analysis
allZ_meta <- three_Z %>%  
  rename(P1 = para$p_1[aID]) %>%
  rename(Z1 = para$Z_1[aID]) %>%
  rename(P2 = para$p_2[aID]) %>%
  rename(Z2 = para$Z_2[aID]) %>%
  select(P1,Z1,P2,Z2,chrpos)


w1 = sqrt(para$ss_1[aID])
w2 = sqrt(para$ss_2[aID]) 


allZ_meta$tempZ1 = qnorm(allZ_meta$P1/ 2) * sign(allZ_meta$Z1)
allZ_meta$tempZ2 = qnorm(allZ_meta$P2/ 2) * sign(allZ_meta$Z2)

allZ_meta$Z_meta = (allZ_meta$tempZ1*w1 + allZ_meta$tempZ2*w2)/sqrt(w1^2 + w2^2)
allZ_meta$P_meta = 2*pnorm(abs(-allZ_meta$Z_meta), lower.tail=FALSE)
allZ_meta_sig <- allZ_meta %>%
  filter(P_meta < 1e-8) %>%
  mutate(method = "Meta Analysis") %>%
  select(chrpos,method)

#all replication
merge_results_scc <- tempDat %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y,cat) 

merge_results_scc$cat <- forcats::fct_relevel(merge_results_scc$cat, c("Both" ,
                                                               "Meta Analysis",
                                                               "Replication Analysis",
                                                               "Neither"))

manPlot_scc <- plotmanhattan(merge_results_scc,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,5,3)], shapeValues=c(16,17,15,18),
                             ylimits=c(0, 20), legName ="Method",title = "LUSC" )



ggsave(paste0(outputDir_Fig,"/two_rep_scc.png"), plot = manPlot_scc , width = 6,
       height = 3)




########adeno#############
three_Z <- fread(here::here(dataDir, "lc_adeno_three.txt"))
para <- data.frame(
  p_1 = c("pILCCO","pILCCO","pUKB"),
  p_2 = c("pUKB","pMVP","pMVP"),
  Z_1 = c("zILCCO","zILCCO","zUKB"),
  Z_2 = c("zUKB","zMVP","zMVP"),
  ss_1 = c(61073,61073,332422),
  ss_2 = c(332422,14132,14132))

fdrLimitNew <- 0.1

aID = 2

tempNew <- fread(paste0(outputDir ,"/two_rep_adeno_aID",aID,"_newlfdr.txt"), header=T, data.table=F)


tempDat <- three_Z %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)

#meta analysis
allZ_meta <- three_Z %>%  
  rename(P1 = para$p_1[aID]) %>%
  rename(Z1 = para$Z_1[aID]) %>%
  rename(P2 = para$p_2[aID]) %>%
  rename(Z2 = para$Z_2[aID]) %>%
  select(P1,Z1,P2,Z2,chrpos)


w1 = sqrt(para$ss_1[aID])
w2 = sqrt(para$ss_2[aID]) 


allZ_meta$tempZ1 = qnorm(allZ_meta$P1/ 2) * sign(allZ_meta$Z1)
allZ_meta$tempZ2 = qnorm(allZ_meta$P2/ 2) * sign(allZ_meta$Z2)

allZ_meta$Z_meta = (allZ_meta$Z1*w1 + allZ_meta$Z2*w2)/sqrt(w1^2 + w2^2)
allZ_meta$P_meta = 2*pnorm(abs(-allZ_meta$Z_meta), lower.tail=FALSE)
allZ_meta_sig <- allZ_meta %>%
  filter(P_meta < 1e-8) %>%
  mutate(method = "Meta Analysis") %>%
  select(chrpos,method)

#all replication
merge_results_adeno <- tempDat %>%
  mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos,1,0)) %>%
  mutate(replication = ifelse(chrpos %in% rejectDat$chrpos,1,0)) %>% 
  mutate(cat = "") %>%  # Initialize cat column
  mutate(cat = ifelse(meta == 1 & replication == 1,"Both",cat)) %>%
  mutate(cat = ifelse(meta == 1 & replication == 0,"Meta Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 1,"Replication Analysis",cat)) %>%
  mutate(cat = ifelse(meta == 0 & replication == 0,"Neither",cat)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y,cat) 

merge_results_adeno$cat <- forcats::fct_relevel(merge_results_adeno$cat, c("Both" ,
                                                               "Meta Analysis",
                                                               "Replication Analysis",
                                                               "Neither"))




manPlot_adeno <- plotmanhattan(merge_results_adeno,build = "hg19", colValues = wes_palette("Darjeeling1")[c(1,2,3)], shapeValues=c(16,17,18),
                               ylimits=c(0, 25), legName ="Method",title = "LUAD" )

ggsave(paste0(outputDir_Fig,"/two_rep_adeno.png"), plot = manPlot_adeno, width = 6,
       height = 3)




# 
# mergePlot1 <- ggarrange(manPlot_overall, manPlot_scc, manPlot_adeno, ncol=3, nrow=1,common.legend = T,labels = c("A", "B", "C"))
# ggsave(paste0(outputDir_Fig,"/test.png"), plot = mergePlot1, width = 10,
#        height = 5)


