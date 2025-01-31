# Make Figure and Table for replication analysis (Lung)

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/plot_replication.R")

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


# 4 is replication ILCCO overall and UKB lc

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
outputDir_Fig <- here::here("Lung", "Fig_Tab")
dataDir <- here::here("Data")

# for colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#----------------------------------------------------------------------------------------#


### replication plot for Lung

# add position information to data
s4 <- fread(here::here(outputDir, "reject_bmi_with_overall_neg5_reject_S_1_aID4.txt"))
s4new <- s4 %>% filter(rejNew == 1) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))

# for plotting axes
allZ <- fread(here::here(dataDir, "replication_with_lcoverall.txt")) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP))
chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- allZ %>% filter(Chr == chr_it)
  maxPos <- max(tempDat$BP)
  chrCounts[chr_it + 1] <- maxPos
}

#Significant SNPs from McKay 2017 Table 2
toDoSNP <- fread(here::here(dataDir, "toDoSNP.txt")) 

# data for manhattan plot - replication
manDataRep <- s4new %>% 
  select(Chr, BP, chrpos, newLfdr, cumNew) %>% 
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  mutate(RepLC = ifelse(chrpos %in% toDoSNP$chrpos, 1, 0)) %>%
  left_join(toDoSNP, by = c("chrpos","Chr","BP")) %>%
  mutate(cat = case_when(
    Cancer == "Overall" ~ "Lung",
    Cancer == "Adeno" ~ "Adeno",
    Cancer == "Squam" ~ "SQC",
    TRUE ~ "Not in McKay's Paper"
  )) %>%
  select(Chr, BP, chrpos, newLfdr,cat,cumNew) %>%
  # distinct() keeps the first one, so we arrange first - want to show Pleio lfdr if in pleio
  distinct(., chrpos, .keep_all = TRUE)  




# plot manhattan function
# plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
#   # arrange data by chromosome
#   plotRes <- plotRes %>% arrange(Chr)
#   uniqueChrs <- sort(unique(plotRes$Chr))
#   
#   # add true positions
#   truePos <- rep(NA, nrow(plotRes))
#   counter <- 1
#   for (chr_it in 1:length(uniqueChrs)) {
#     tempChr <- uniqueChrs[chr_it]
#     tempDat <- plotRes %>% filter(Chr == tempChr)
#     truePos[counter:(counter + nrow(tempDat) - 1)] <- rep(sum(chrCounts[1:tempChr]), nrow(tempDat)) + tempDat$BP
#     counter <- counter + nrow(tempDat)
#   }
#   
#   # plot
#   xBreaks <- cumsum(chrCounts[-1])
#   xBreaksLabs <- 1:22
#   xBreaksLabs[c(9, 11, 13, 15, 16, 17, 19, 20, 21)] <- ""
#   
#   plotDat <- plotRes %>% mutate(truePos = truePos)
#   
#   returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr), color=as.factor(cat), shape=as.factor(cat))) +
#     geom_point() +
#     xlab("Chromosome") + ylab("-log10(lfdr)") +
#     #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
#     scale_color_manual(name=legName, values=colValues) +
#     scale_shape_manual(name=legName, values=shapeValues) +
#     scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
#     ylim(ylimits) +
#     theme_cowplot() +
#     theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
#           legend.title=element_text(size=18), legend.text=element_text(size=12)) +
#     guides(colour = guide_legend(override.aes = list(size=4)))
#   
#   
# }
# 
# 
# manPlotRep <- plotManhattan(plotRes = manDataRep, chrCounts,
#                             colValues=wes_palette("GrandBudapest2"), shapeValues=c(17, 18, 16, 15),
#                             ylimits=c(0, 12.5), legName="Lung Cancer")
# manPlotRep
# 






#make manhattan plot (Significant csmGmm only (lfdr))
plotDat_lfdr <- manDataRep %>%
  mutate(y = -log10(newLfdr)) %>%
  mutate(chrom= paste0("chr",Chr))
colnames(plotDat_lfdr)[2] = "pos"

plotmanhattan = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
 
  names(x_breaks)[c(9, 11, 13, 16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y)) +
    geom_point(color = colValues[1]) +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(lfdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    #scale_color_manual(name=legName, values=colValues) +
    #scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}


manPlotRep <- plotmanhattan(plotDat_lfdr,build = "hg19", colValues=wes_palette("Zissou1"), shapeValues=c(17),
               ylimits=c(0, 12.5), legName ="",title = "Significant Variants of Replication Analysis" )

manPlotRep
ggsave(paste0(outputDir_Fig,"/Sig_Rep_Manhatton_lfdr.png"), plot = manPlotRep, width = 6,
       height = 3 )


#make manhattan plot (Significant csmGmm only (Fdr))
plotDat_fdr <- manDataRep %>%
  mutate(y = -log10(cumNew)) %>%
  mutate(chrom= paste0("chr",Chr))
colnames(plotDat_fdr)[2] = "pos"

plotmanhattan = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  names(x_breaks)[c(9, 11, 13, 16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y)) +
    geom_point(color = colValues[1]) +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    #scale_color_manual(name=legName, values=colValues) +
    #scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}


manPlotRep <- plotmanhattan(plotDat_fdr,build = "hg19", colValues=wes_palette("Zissou1"), shapeValues=c(17),
                            ylimits=c(0, 12.5), legName ="",title = "Significant Variants of Replication Analysis" )

manPlotRep
ggsave(paste0(outputDir_Fig,"/Sig_Rep_Manhatton_fdr.png"), plot = manPlotRep, width = 6,
       height = 3 )

#----------------------------------------------------------------------------------------#


### replication table for Lung
#duplicated with McKay's
tab1 <- manDataRep %>%
  filter(cat %in% c("Lung", "Adeno", "SQC"))  %>%
  left_join(toDoSNP[,c("RS", "Gene", "chrpos")], by = "chrpos") %>%
  left_join(allZ[,c("chrpos", "Zoverall", "Zlcukb", "pOverall", "pLCukb")], by = "chrpos")

# top 5 lowest new LFDR
tab2 <- manDataRep %>%
  left_join(allZ[,c("chrpos", "Zoverall", "Zlcukb", "pOverall", "pLCukb")], by = "chrpos") %>%
  arrange(newLfdr) %>%
  slice(1:5) %>%
  select("chrpos", "newLfdr", "cat", "Zoverall", "Zlcukb", "pOverall", "pLCukb")

write.table(tab1, paste0(outputDir, "/tab1_replication.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(tab2, paste0(outputDir, "/tab2_replication.txt"), append=F, quote=F, row.names=F, col.names=T)
#----------------------------------------------------------------------------------------#


### scatter plot for P and Z

Fig2_data <- manDataRep %>%
  left_join(allZ[,c("chrpos", "Zoverall", "Zlcukb", "pOverall", "pLCukb")], by = "chrpos")

ggplot(Fig2_data, aes(x=pOverall, y=pLCukb, color=cat)) + 
  geom_point() 


ggplot(Fig2_data, aes(x=Zlcukb, y=Zoverall, color=cat)) + 
  geom_point() 

sig_Z = allZ %>%
  filter(pOverall < 5e-5 )

ggplot(sig_Z, aes(x=Zoverall, y=Zlcukb)) + 
  geom_point() 



#----------------------------------------------------------------------------------------#
#manhattan plot for meta-analysis
meta_analysis <- function(data,n1,n2){
  
  w1 = sqrt(n1) 
  w2 = sqrt(n2) 
  
  data$tempZ1 = qnorm(data$P1/ 2) * sign(data$Z1)
  data$tempZ2 = qnorm(data$P2/ 2) * sign(data$Z2)
  data$Z_meta = (data$tempZ1*w1 + data$tempZ2*w2)/sqrt(w1^2 + w1^2)
  data$P_meta =2*pnorm(abs(-data$Z_meta), lower.tail=FALSE)
  return(data)
}

allZ_meta <- allZ %>%
  rename(P1 = pLCukb) %>%
  rename(Z1 = Zlcukb) %>%
  rename(P2 = pOverall) %>%
  rename(Z2 = Zoverall)

meta_sig <- meta_analysis(allZ_meta,332000,85716) 
meta_sig <- meta_sig %>% filter(P_meta < 1e-8)


plotDatMeta <- meta_sig %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  mutate(chrom= paste0("chr",Chr))  %>%
  rename(pos = BP)


plotmanhattan_meta = function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  names(x_breaks)[c(9, 11, 13, 16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y = y)) +
    geom_point(color = colValues[1]) +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    #scale_color_manual(name=legName, values=colValues) +
    #scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}

manPlotmeta <- plotmanhattan_meta(plotDatMeta,build = "hg19",colValues=wes_palette("GrandBudapest2"), shapeValues=c(17),
                              ylimits=c(0, 25), legName = "" ,title = "Significant Variants of Meta Analysis")

manPlotmeta
ggsave(paste0(outputDir_Fig,"/Sig_Meta_Manhatton_fdr.png"), plot = manPlotmeta, width = 6,
       height = 3 )





#----------------------------------------------------------------------------------------#
#manhattan plot - replication and meta
plotDat_fdr <- manDataRep %>%
  mutate(y = -log10(cumNew)) %>%
  mutate(chrom= paste0("chr",Chr))
colnames(plotDat_fdr)[2] = "pos"

manDataMetaRep <- rbind(
  plotDatMeta %>% select(chrpos,chrom,pos,y) %>%   
    mutate(method = "Meta Analysis"),
  plotDat_fdr %>% select(chrpos,chrom,pos,y) %>%
  mutate(method = "Replication Analysis")) %>%
  #distinct() keeps the first one, so we arrange first, and show Meta Analysis first
  arrange(chrpos) %>% 
  distinct(., chrpos, .keep_all = TRUE)


plotmanhattanMetaRep <- function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  names(x_breaks)[c(9, 11, 13, 16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=as.factor(method), shape=as.factor(method))) +
    geom_point() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}

manPlotMetaRep <- plotmanhattanMetaRep(manDataMetaRep,build = "hg19",colValues= c(wes_palette("GrandBudapest2")[1],wes_palette("Zissou1")[1]), shapeValues=c(16,16),
                              ylimits=c(0, 25), legName="Methods", title = "Replication Analysis vs Meta Analysis")

manPlotMetaRep
ggsave(paste0(outputDir_Fig,"/Sig_Rep_Meta_Manhattan_fdr.png"), plot = manPlotMetaRep, width = 6,
       height = 3 )


#----------------------------------------------------------------------------------------#
#FAVOR manhatton plot Chr15
manDataFavor <- readxl::read_xlsx(here::here(outputDir, "FAVOR_results.xlsx")) %>%
  rename(Chr = Chromosome, BP = Position)

manDataFavor <- manDataRep %>%
  left_join(manDataFavor,by = c("Chr", "BP")) %>% 
  mutate(AC_rank = case_when(
  is.na(ApcConservationV2) ~ "NoInfor",
  ApcConservationV2 > 25 ~ "high (>25)",
  ApcConservationV2 > 10 & ApcConservationV2 <= 25 ~ "medium (10 - 25)",
  TRUE ~ "low (< 10)"
)) %>% 
  mutate(AE_rank = case_when(
    is.na(ApcEpigeneticsActive) ~ "NoInfor",
    ApcEpigeneticsActive > 25 ~ "high (>25)",
    ApcEpigeneticsActive > 10 & ApcEpigeneticsActive <= 25 ~ "medium (10 - 25)",
    TRUE ~ "low (< 10)"
  ))
  

#Chromosome 15 (fdr) AE
plotDatFavor <- manDataFavor %>% select(Chr, BP, newLfdr,cumNew, AC_rank, AE_rank) %>%
  mutate(y = -log10(cumNew)) %>%
  rename(chrom = Chr, pos = BP, y = y) %>%
  mutate(chrom = paste0("chr",chrom)) %>%
  filter(chrom == "chr15")
# plotDatFavor$AC_rank <- factor(x = plotDatFavor$AC_rank, levels = c("high (>25)" ,
#                                                                     "medium (10 - 25)",
#                                                                     "low (< 10)",
#                                                                     "NoInfor"))
# plotDatFavor$AE_rank <- factor(x = plotDatFavor$AE_rank, levels = c("high (>25)" ,
#                                                                     "medium (10 - 25)",
#                                                                     "low (< 10)",
#                                                                     "NoInfor"))

plotDatFavor$AC_rank <- forcats::fct_relevel(plotDatFavor$AC_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "NoInfor"))

plotDatFavor$AE_rank <- forcats::fct_relevel(plotDatFavor$AE_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "NoInfor"))
plotDatFavor <- plotDatFavor[order(plotDatFavor$AE_rank, decreasing=TRUE), ]
plotDatFavor <- plotDatFavor[order(plotDatFavor$AC_rank, decreasing=TRUE), ]


plotmanhattanFavor_AE <- function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AE_rank, shape=AE_rank)) +
    geom_point() +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(2307285719,2500171864),expand=c(0.001,0),breaks=x_breaks[c(14,15,16)],
                       labels=names(x_breaks[c(14,15,16)]),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}

manPlotFavor_AE<- plotmanhattanFavor_AE(plotDatFavor,build = "hg19",colValues=wes_palette("Darjeeling1")[c(2,3,4)], shapeValues=c(16,16,16,16),
                                       ylimits=c(1, 6), legName="AE score", title = "Manhattan plot with AE score")

manPlotFavor_AE
ggsave(paste0(outputDir_Fig,"/Sig_chr15_FAVOR_AE.png"), plot = manPlotFavor_AE, width = 3,
       height = 3 )




#Chromosome 15 (fdr) AC

plotmanhattanFavor_AC <- function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=as.factor(AC_rank), shape=as.factor(AC_rank))) +
    geom_point() +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(2307285719,2500171864),expand=c(0.001,0),breaks=x_breaks[c(14,15,16)],
                       labels=names(x_breaks[c(14,15,16)]),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}

manPlotFavor_AC<- plotmanhattanFavor_AC(plotDatFavor,build = "hg19",colValues=wes_palette("Darjeeling1"), shapeValues=c(16,16,16,16),
                                  ylimits=c(1, 6), legName="AC score", title = "Manhattan plot with AC score")

manPlotFavor_AC
ggsave(paste0(outputDir_Fig,"/Sig_chr15_FAVOR_AC.png"), plot = manPlotFavor_AC, width = 3,
       height = 3 )








#FAVOR chr 5
plotDatFavor <- manDataFavor %>% select(Chr, BP, newLfdr,cumNew, AC_rank, AE_rank) %>%
  mutate(y = -log10(cumNew)) %>%
  rename(chrom = Chr, pos = BP, y = y) %>%
  mutate(chrom = paste0("chr",chrom)) %>%
  filter(chrom == "chr5")
plotDatFavor$AC_rank <- forcats::fct_relevel(plotDatFavor$AC_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "NoInfor"))

plotDatFavor$AE_rank <- forcats::fct_relevel(plotDatFavor$AE_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "NoInfor"))
plotDatFavor <- plotDatFavor[order(plotDatFavor$AE_rank, decreasing=TRUE), ]
plotDatFavor <- plotDatFavor[order(plotDatFavor$AC_rank, decreasing=TRUE), ]

plotmanhattanFavor_AE <- function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AE_rank, shape=AE_rank)) +
    geom_point() +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(786049562,972084330),expand=c(0.001,0),breaks=x_breaks[c(4,5)],
                       labels=names(x_breaks[c(4,5)]),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}

manPlotFavor_AE<- plotmanhattanFavor_AE(plotDatFavor,build = "hg19",colValues=wes_palette("Darjeeling1"), shapeValues=c(16,16,16,16),
                                        ylimits=c(1, 6), legName="AE score", title = "Manhattan plot with AE score")



#CHr5 AC

plotmanhattanFavor_AC <- function(gwas,build=c('hg18','hg19','hg38'),colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AC_rank, shape=AC_rank)) +
    geom_point() +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(786049562,972084330),expand=c(0.001,0),breaks=x_breaks[c(4,5)],
                       labels=names(x_breaks[c(4,5)]),name='Chromosome') +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ggtitle(title)
  
}

manPlotFavor_AC<- plotmanhattanFavor_AC(plotDatFavor,build = "hg19",colValues=wes_palette("Darjeeling1")[c(2,3,4)], shapeValues=c(16,16,16,16),
                                        ylimits=c(1, 6), legName="AC score", title = "Manhattan plot with AE score")

manPlotFavor_AC
ggsave(paste0(outputDir_Fig,"/Sig_chr5_FAVOR_AC.png"), plot = manPlotFavor_AC, width = 3,
       height = 3 )

#----------------------------------------------------------------------------------------#
#Table for simulation
final_result <- data.frame(nCausal =NA, propCausal = NA ,nRejNew = NA, nRejMeta = NA, nRejThreshold = NA,
                           powerNew = NA, powerMeta = NA, powerThreshold = NA,
                           fdpNew = NA, fdpMeta = NA,fdpThreshold = NA)

sim_results <- fread(paste0(outputDir, "/rep_sim_aID", 1, ".txt")) %>%
  select("nCausal","propCausal", "nRejNew", "nRejMeta", "nRejThreshold", 
         "powerNew","powerMeta","powerThreshold", 
         "fdpNew", "fdpMeta" ,"fdpThreshold" )

for (aID in 2:100) {
  tempResults <- fread(paste0(outputDir, "/rep_sim_aID", aID, ".txt")) %>%
    select("nCausal","propCausal", "nRejNew", "nRejMeta", "nRejThreshold", 
           "powerNew","powerMeta","powerThreshold", 
           "fdpNew", "fdpMeta" ,"fdpThreshold" )
  sim_results <- rbind(sim_results,tempResults)
}

pCausal1 <- subset(sim_results,sim_results$propCausal == 2e-04)
pCausal1 <- sapply(pCausal1,mean)
pCausal2 <- subset(sim_results,sim_results$propCausal == 2e-03)
pCausal2 <- sapply(pCausal2,mean)
pCausal3 <- subset(sim_results,sim_results$propCausal == 1e-02)
pCausal3 <- sapply(pCausal3,mean)

final_result[1,] = pCausal1
final_result[2,] = pCausal2
final_result[3,] = pCausal3
write.table(final_result, paste0(outputDir_Fig, "/Tab_avg_sim_results.txt"), append=F, quote=F, row.names=F, col.names=T)

#----------------------------------------------------------------------------------------#
#Fig for simulation
#Box plot for power
plotDatPower <- sim_results %>%
  select(nCausal,powerNew,powerMeta,powerThreshold) %>%
  mutate(nCausal = as.character(nCausal)) %>%
  melt(id = c("nCausal")) %>% 
  ggplot(aes(x=nCausal, y=value)) + geom_boxplot(aes(fill=variable)) + ylim(0,1)
plotDatPower
ggsave(paste0(outputDir_Fig,"/boxplot_power.png"), plot = plotDatPower, width = 6,
       height = 3 )

#Box plot for fdp
plotDatFDP <- sim_results %>%
  select(nCausal,fdpNew,fdpMeta,fdpThreshold) %>%
  mutate(nCausal = as.character(nCausal)) %>%
  melt(id = c("nCausal")) %>% 
  ggplot(aes(x=nCausal, y=value)) + geom_boxplot(aes(fill=variable)) 
plotDatFDP
ggsave(paste0(outputDir_Fig,"/boxplot_fdp.png"), plot = plotDatFDP, width = 6,
       height = 3 )



#----------------------------------------------------------------------------------------#
#manhattan plot for all SNPs
#prepare all meta results
allZ_meta <- allZ %>%
  rename(P1 = pLCukb) %>%
  rename(Z1 = Zlcukb) %>%
  rename(P2 = pOverall) %>%
  rename(Z2 = Zoverall)

meta_results <- meta_analysis(allZ_meta,332000,85716) 


plotDatMeta <- meta_results %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  mutate(chrom= paste0("chr",Chr))  %>%
  rename(pos = BP) %>%
  select(chrom,pos,y)


#prepare all replication results

rep_results <- tempDat %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(chrom = paste0("chr",Chr), pos = as.numeric(BP), y = -log10(cumNew)) %>%
  select(chrom,pos,y)

#manhattan plot for all Meta
manhattan_mod = function(gwas,build=c('hg18','hg19','hg38'),color1=wes_palette("Zissou1")[1],color2=wes_palette("Zissou1")[2]){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  data=add_color(data,color1 = color1,color2 = color2)
  data=add_shape(data,shape=16)
  data=add_fill(data)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  color_map=unique(data$color)
  names(color_map)=unique(data$color)
  
  ggplot2::ggplot(data,aes(x=cumulative_pos,y=y,color=color,shape=shape,fill=fill))+
    geom_point()+
    theme_classic()+
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome')+
    scale_y_continuous(limits=c(0,max(data$y)+1),expand=c(0.01,0),name=expression(paste("-",log[10],'(Fdr)')))+
    scale_color_manual(values=color_map,guide='none')+
    scale_fill_manual(values = color_map, guide = 'none')+
    scale_shape_identity()
}
plotDatMeta_all <- manhattan_mod(plotDatMeta,build='hg19')
ggsave(paste0(outputDir_Fig,"/man_Meta_all.png"), plot = plotDatMeta_all, width = 6,
       height = 3)
#manhattan plot for all Meta
plotDatRep_all <- manhattan_mod(rep_results,build='hg19')
ggsave(paste0(outputDir_Fig,"/man_Rep_all.png"), plot = plotDatRep_all, width = 6,
       height = 3)
