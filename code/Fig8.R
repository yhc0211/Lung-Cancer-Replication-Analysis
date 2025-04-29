# Make Figure 8: Three-way Meta Analysis Functional Annotation
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/Fig8.R")

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
library(plotgardener)


# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
outputDir_Fig <- here::here("Lung", "Fig_Tab")
dataDir <- here::here("Data")


### meta table for Lung
three_Z <- fread(here::here(dataDir, "lc_overall_three.txt"))
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

results_hg38 <- fread(here::here(outputDir, "hg38_overall_meta.bed"),header=F)

hg_ref <- cbind(allZ_meta_sig,results_hg38) %>%
  select(chrpos,V1) %>%
  separate(V1,c("chrom","BP_hg38"),sep = "-")



manDataFavor <- readxl::read_xlsx(here::here(outputDir, "FAVOR_overall_meta.xlsx")) %>%
  rename(Chr = Chromosome, BP = Position) %>%
  cbind(hg_ref) %>%
  select("chrpos","ApcConservationV2","ApcEpigeneticsActive")


manDataFavor <- allZ_meta_sig %>%
  left_join(manDataFavor,by = "chrpos")


manDataFavor <- manDataFavor %>% 
  mutate(AC_rank = case_when(
    is.na(ApcConservationV2) ~ "Not Significant",
    ApcConservationV2 > 25 ~ "high (>25)",
    ApcConservationV2 > 10 & ApcConservationV2 <= 25 ~ "medium (10 - 25)",
    TRUE ~ "low (< 10)"
  )) %>% 
  mutate(AE_rank = case_when(
    is.na(ApcEpigeneticsActive) ~ "Not Significant",
    ApcEpigeneticsActive > 25 ~ "high (>25)",
    ApcEpigeneticsActive > 10 & ApcEpigeneticsActive <= 25 ~ "medium (10 - 25)",
    TRUE ~ "low (< 10)"
  )) 


manDataFavor$AC_rank <- forcats::fct_relevel(manDataFavor$AC_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "Not Significant"))

manDataFavor$AE_rank <- forcats::fct_relevel(manDataFavor$AE_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "Not Significant"))
manDataFavor <- manDataFavor[order(manDataFavor$AE_rank, decreasing=TRUE), ]
manDataFavor <- manDataFavor[order(manDataFavor$AC_rank, decreasing=TRUE), ]

#Chromosome 15  AE
plotDatFavor <- manDataFavor %>% 
  separate(chrpos, c("Chr","pos"),sep=":") %>% 
  mutate(pos = as.numeric(pos)) %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  mutate(chrom = paste0("chr",Chr)) %>%
  filter(chrom == "chr15")



plotmanhattanFavor_AE <- function(gwas,build=c('hg18','hg19','hg38'), colValues, shapeValues, ylimits, legName){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  

  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AE_rank, shape=AE_rank)) +
    geom_point(size = 2) +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(2385997522- 20000, max(data$cumulative_pos) + 20000),expand=c(0.001,0),breaks=x_breaks[c(14,15,16)],
                       labels=names(x_breaks[c(14,15,16)]),name='Chromosome 15') +
    # scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
    #                    labels=names(x_breaks),name='Chromosome')+
    ylim(ylimits) +
    #xlim(c(min(data$cumulative_pos),max(data$cumulative_pos)))+
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) 
  
}

manPlotFavor_AE<- plotmanhattanFavor_AE(plotDatFavor,build = "hg19",colValues=wes_palette("FantasticFox1")[c(5,3,4)],shapeValues=c(16,17,18),
                                        ylimits=c(0, 70), legName="Epigenetics Score")

# ggsave(paste0(outputDir_Fig,"/Fig7_1.png"), plot = manPlotFavor_AE, width = 6,
#         height = 6)


#Chromosome 15  AC
plotmanhattanFavor_AC <- function(gwas,build=c('hg18','hg19','hg38'), colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AC_rank, shape=AC_rank)) +
    geom_point(size = 2) +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(2385997522- 20000, max(data$cumulative_pos) + 20000),expand=c(0.001,0),breaks=x_breaks[c(14,15,16)],
                       labels=names(x_breaks[c(14,15,16)]),name='Chromosome 15') +
    # scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
    #                    labels=names(x_breaks),name='Chromosome')+
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) 
  
}

manPlotFavor_AC<- plotmanhattanFavor_AC(plotDatFavor,build = "hg19",colValues=wes_palette("FantasticFox1")[c(5,3,4)],shapeValues=c(16,17,18),
                                        ylimits=c(0, 70), legName="Conservation Score")

# ggsave(paste0(outputDir_Fig,"/Fig7_2.png"), plot = manPlotFavor_AC, width = 6,
#         height = 6)





### replication table for SCC Lung
three_Z <- fread(here::here(dataDir, "lc_scc_three.txt"))
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

results_hg38 <- fread(here::here(outputDir, "hg38_scc_meta.bed"),header=F)

hg_ref <- cbind(allZ_meta_sig,results_hg38) %>%
  select(chrpos,V1) %>%
  separate(V1,c("chrom","BP_hg38"),sep = "-")


manDataFavor <- readxl::read_xlsx(here::here(outputDir, "FAVOR_scc_meta.xlsx")) %>%
  rename(Chr = Chromosome, BP = Position) %>%
  cbind(hg_ref) %>%
  select("chrpos","ApcConservationV2","ApcEpigeneticsActive")


manDataFavor <- allZ_meta_sig %>%
  left_join(manDataFavor,by = "chrpos")


manDataFavor <- manDataFavor %>%
  mutate(AC_rank = case_when(
    is.na(ApcConservationV2) ~ "Not Significant",
    ApcConservationV2 > 25 ~ "high (>25)",
    ApcConservationV2 > 10 & ApcConservationV2 <= 25 ~ "medium (10 - 25)",
    TRUE ~ "low (< 10)"
  )) %>%
  mutate(AE_rank = case_when(
    is.na(ApcEpigeneticsActive) ~ "Not Significant",
    ApcEpigeneticsActive > 25 ~ "high (>25)",
    ApcEpigeneticsActive > 10 & ApcEpigeneticsActive <= 25 ~ "medium (10 - 25)",
    TRUE ~ "low (< 10)"
  ))


manDataFavor$AC_rank <- forcats::fct_relevel(manDataFavor$AC_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "Not Significant"))

manDataFavor$AE_rank <- forcats::fct_relevel(manDataFavor$AE_rank, c("high (>25)" ,
                                                                     "medium (10 - 25)",
                                                                     "low (< 10)",
                                                                     "Not Significant"))
manDataFavor <- manDataFavor[order(manDataFavor$AE_rank, decreasing=TRUE), ]
manDataFavor <- manDataFavor[order(manDataFavor$AC_rank, decreasing=TRUE), ]

#Chromosome 15  AE
plotDatFavor <- manDataFavor %>% 
  separate(chrpos, c("Chr","pos"),sep=":") %>% 
  mutate(pos = as.numeric(pos)) %>%
  mutate(y = -log10(p.adjust(P_meta, method="BH"))) %>%
  mutate(chrom = paste0("chr",Chr)) %>%
  filter(chrom == "chr15")



plotmanhattanFavor_AE <- function(gwas,build=c('hg18','hg19','hg38'), colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AE_rank, shape=AE_rank)) +
    geom_point(size = 2) +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(min(data$cumulative_pos)- 20000, max(data$cumulative_pos) + 20000),expand=c(0.001,0),breaks=x_breaks[c(14,15,16)],
                       labels=names(x_breaks[c(14,15,16)]),name='Chromosome 15') +
    # scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
    #                    labels=names(x_breaks),name='Chromosome')+
    ylim(ylimits) +
    #xlim(c(min(data$cumulative_pos),max(data$cumulative_pos)))+
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) 
  
}

manPlotFavor_AE_SCC<- plotmanhattanFavor_AE(plotDatFavor,build = "hg19",colValues=wes_palette("FantasticFox1")[c(5,3,4)],shapeValues=c(16,17,18),
                                            ylimits=c(0, 20), legName="Epigenetics Score")

# ggsave(paste0(outputDir_Fig,"/Fig7_3.png"), plot = manPlotFavor_AE_SCC, width = 6,
#         height = 6)


#Chromosome 15  AC
plotmanhattanFavor_AC <- function(gwas,build=c('hg18','hg19','hg38'), colValues, shapeValues, ylimits, legName, title){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  chrom_lengths=get_chrom_lengths(build)
  #chrom_lengths =  chrom_lengths[c("chr14","chr15","chr16")]
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  #names(x_breaks)[c(9, 11, 13,  16, 17, 20, 21)] <- ""
  
  returnPlot <- ggplot(data, aes(x=cumulative_pos, y=y, color=AC_rank, shape=AC_rank)) +
    geom_point(size = 2) +
    #geom_jitter() +
    xlab("Chromosome") + ylab(expression(paste("-",log[10],'(Fdr)'))) +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(limits=c(min(data$cumulative_pos) - 20000, max(data$cumulative_pos) + 20000), expand=c(0.001,0),breaks=x_breaks[c(14,15,16)],
                       labels=names(x_breaks[c(14,15,16)]),name='Chromosome 15') +
    # scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
    #                    labels=names(x_breaks),name='Chromosome')+
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.title=element_text(size=12), legend.text=element_text(size=12)) +
    guides(colour = guide_legend(override.aes = list(size=4))) 
  
}

manPlotFavor_AC_SCC<- plotmanhattanFavor_AC(plotDatFavor,build = "hg19",colValues=wes_palette("FantasticFox1")[c(5,3,4)], shapeValues=c(16,17,18),
                                            ylimits=c(0, 20), legName="Conservation Score")

# ggsave(paste0(outputDir_Fig,"/Fig6_4.png"), plot = manPlotFavor_AC_SCC, width = 6,
#        height = 6)
merge_plot1 <- ggarrange(manPlotFavor_AC,manPlotFavor_AC_SCC,ncol=2, nrow=1,common.legend = TRUE  ,legend="bottom")
merge_plot2 <- ggarrange(manPlotFavor_AE,manPlotFavor_AE_SCC,ncol=2, nrow=1,common.legend = TRUE,  legend="bottom")
merge_plot <- ggarrange(merge_plot1,merge_plot2,ncol=1, nrow=2,common.legend = TRUE  ,legend="bottom")


pageCreate(width = 12, height = 6.5, default.units = "inches")

plotGG(
  plot = merge_plot,
  x = 0.1, y = 0.25,
  width = 12, height = 6, just = c("left", "top")
)

plotText(label = "A", x = 0.1, y = 0.3,
         fontsize = 10, fontface = "bold", just = "center",
         default.units = "inches")
plotText(label = "B", x = 6, y = 0.3,
         fontsize = 10, fontface = "bold", just = "center",
         default.units = "inches")
plotText(label = "C", x = 0.1, y = 3,
         fontsize = 10, fontface = "bold", just = "center",
         default.units = "inches")
plotText(label = "D", x = 6, y = 3,
         fontsize = 10, fontface = "bold", just = "center",
         default.units = "inches")
pageGuideHide()


ggsave(paste0(outputDir_Fig,"/Fig8.png"), plot = merge_plot, width = 12,
       height = 6)
