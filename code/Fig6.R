# Make Figure 6: Three-way Replication Analysis and Meta Analysis Results
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/Fig6.R")

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

######Overall#########
three_Z <- fread(here::here(dataDir, "lc_overall_three.txt")) 


tempNew <- fread(paste0(outputDir ,"/three_rep_data_aID1_newlfdr.txt"), header=T, data.table=F)
tempDat <- three_Z %>% select(chrpos,zILCCO, zUKB, zMVP) %>%
  mutate(origIdx = 1:nrow(.)) %>%   
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
  arrange(origIdx)#all replication

rejectDat <- tempDat %>% filter(rejNew == 1)


rep_results <- tempDat %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y)
manhattan_mod = function(gwas,build=c('hg18','hg19','hg38'),title,color1=wes_palette("Zissou1")[1],color2=wes_palette("Zissou1")[3]){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  data=add_color(data,color1 = color1,color2 = color2)
  data=add_shape(data,shape=16)
  data=add_fill(data)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  names(x_breaks)[c(9, 11, 13, 14,16, 17,18, 20, 21)] <- ""
  
  color_map=unique(data$color)
  names(color_map)=unique(data$color)
  
  ggplot2::ggplot(data,aes(x=cumulative_pos,y=y,color=color,shape=shape,fill=fill))+
    geom_point()+
    theme_classic()+
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome')+
    scale_y_continuous(limits=c(0,max(data$y)+1),expand=c(0.01,0),name=expression(paste("-",log[10],'(FDR)')))+
    scale_color_manual(values=color_map,guide='none')+
    scale_fill_manual(values = color_map, guide = 'none')+
    
    scale_shape_identity() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold",size =12))
}
plotDatRep_all <- manhattan_mod(rep_results,build='hg19',title = "Overall Lung Cancer Replication Analysis")
plotDatRep_all <- plotDatRep_all + geom_hline(yintercept = 1, linetype=2, color="grey") 
# ggsave(paste0(outputDir_Fig,"/overall_replication.png"), plot = plotDatRep_all , width = 6,
#        height = 3)

#all meta
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


plotDatMeta <- allZ_meta %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  separate(chrpos, c('Chr', 'BP'), sep=':') %>%
  mutate(pos = as.numeric(BP))  %>%
  mutate(chrom= paste0("chr",Chr))  %>%
  select(chrom,pos,y)

plotDatMeta_all <- manhattan_mod(plotDatMeta,build='hg19',title = "Overall Lung Cancer Meta Analysis")
plotDatMeta_all <- plotDatMeta_all + geom_hline(yintercept = 8, linetype=2, color="grey") 
# ggsave(paste0(outputDir_Fig,"/overall_meta.png"), plot = plotDatMeta_all, width = 6,
#          height = 3)



#########SCC##############
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

rep_results <- tempDat %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y)

plotDatRep_scc <- manhattan_mod(rep_results,build='hg19',title = "LUSC Replication Analysis")
plotDatRep_scc <- plotDatRep_scc + geom_hline(yintercept = 1, linetype=2, color="grey") 


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


plotDatMeta <- allZ_meta %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  separate(chrpos, c('Chr', 'BP'), sep=':') %>%
  mutate(pos = as.numeric(BP))  %>%
  mutate(chrom= paste0("chr",Chr))  %>%
  select(chrom,pos,y)

plotDatMeta_scc <- manhattan_mod(plotDatMeta,build='hg19',title = "LUSC Meta Analysis")
plotDatMeta_scc <- plotDatMeta_scc + geom_hline(yintercept = 8, linetype=2, color="grey") 
# ggsave(paste0(outputDir_Fig,"/scc_meta.png"), plot = plotDatMeta_scc, width = 6,
#         height = 3)


######Adeno########
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

#all replication
rep_results <- tempDat %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y)

plotDatRep_adeno <- manhattan_mod(rep_results,build='hg19',title = "LUAD Replication Analysis")
plotDatRep_adeno <- plotDatRep_adeno + geom_hline(yintercept = 1, linetype=2, color="grey") 

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
#all meta

plotDatMeta <- allZ_meta %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  separate(chrpos, c('Chr', 'BP'), sep=':') %>%
  mutate(pos = as.numeric(BP))  %>%
  mutate(chrom= paste0("chr",Chr))  %>%
  select(chrom,pos,y)

plotDatMeta_adeno <- manhattan_mod(plotDatMeta,build='hg19',title = "LUAD Meta Analysis")
plotDatMeta_adeno<- plotDatMeta_adeno + geom_hline(yintercept = 8, linetype=2, color="grey") 

# 
# mergePlot1 <- ggarrange(plotDatRep_all, plotDatMeta_all, ncol=2, nrow=1,labels = c("A","B"))
# mergePlot2 <- ggarrange(plotDatRep_scc, plotDatMeta_scc, ncol=2, nrow=1,labels = c("C","D"))
# mergePlot3 <- ggarrange(plotDatRep_adeno, plotDatMeta_adeno, ncol=2, nrow=1,labels = c("E","F"))

mergePlot <- ggarrange(plotDatRep_all, plotDatMeta_all, 
                       plotDatRep_scc, plotDatMeta_scc,
                       plotDatRep_adeno, plotDatMeta_adeno,
                       ncol=2, nrow=3,labels = c("A","B","C","D","E","F"))


ggsave(paste0(outputDir_Fig,"/Fig6_v2.png"), plot = mergePlot , width = 6,
       height = 9)
