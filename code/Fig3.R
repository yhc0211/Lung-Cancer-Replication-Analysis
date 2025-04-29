# Make Figure 3: Three-way Replication Analysis Simulation of the different mean and proportion of causal variants

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/Fig3.R")

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




# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
outputDir_Fig <- here::here("Lung", "Fig_Tab")
dataDir <- here::here("Data")

# define colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"


final_result <- data.frame(nCausal =NA, propCausal = NA, mu = NA, nRejNew = NA, nRejMeta = NA, nRejThreshold1 = NA, nRejThreshold2 = NA, nRejThreshold3 = NA,
                           powerNew = NA, powerMeta = NA, powerThreshold1 = NA, powerThreshold2 = NA, powerThreshold3 = NA,
                           fdpNew = NA, fdpMeta = NA,fdpThreshold1 = NA,fdpThreshold2 = NA,fdpThreshold3 = NA)



sim_results <- fread(paste0(outputDir, "/rep_three_sim_Snum",1,"_aID", 1, ".txt")) %>%
  select("nCausal","propCausal","mu", "nRejNew", "nRejMeta", "nRejThreshold1", "nRejThreshold2", "nRejThreshold3", 
         "powerNew","powerMeta","powerThreshold1", "powerThreshold2","powerThreshold3",
         "fdpNew", "fdpMeta" ,"fdpThreshold1" ,"fdpThreshold2","fdpThreshold3")



for (Snum in 1:15) {
  for (aID in 1:100) {
    tempResults <- fread(paste0(outputDir, "/rep_three_sim_Snum",Snum,"_aID", aID, ".txt")) %>%
      select("nCausal","propCausal","mu", "nRejNew", "nRejMeta", "nRejThreshold1", "nRejThreshold2", "nRejThreshold3", 
             "powerNew","powerMeta","powerThreshold1", "powerThreshold2","powerThreshold3",
             "fdpNew", "fdpMeta" ,"fdpThreshold1" ,"fdpThreshold2","fdpThreshold3")
    sim_results <- rbind(sim_results,tempResults)
  }
}

sim_results <- sim_results[-1,]

#----------------------------------------------------------------------------------------#
#power plot for different mu 

plot_func <- function(data, pCausal, type){
  
  sim_results = data
  #filter proportion of causal
  sim_results = sim_results %>% filter(sim_results$propCausal == pCausal)
  
  #calculate mean for each mu
  mu1 <- subset(sim_results,sim_results$mu == 2)
  mu1 <- sapply(mu1,mean)
  mu2 <- subset(sim_results,sim_results$mu == 3)
  mu2 <- sapply(mu2,mean)
  mu3 <- subset(sim_results,sim_results$mu == 4)
  mu3 <- sapply(mu3,mean)
  mu4 <- subset(sim_results,sim_results$mu == 5)
  mu4 <- sapply(mu4,mean)
  mu5 <- subset(sim_results,sim_results$mu == 6)
  mu5 <- sapply(mu5,mean)
  
  final_result[1,] = mu1
  final_result[2,] = mu2
  final_result[3,] = mu3
  final_result[4,] = mu4
  final_result[5,] = mu5
  
  
  
  
  if (type == "Power"){
    plotDatPowerMu <- final_result %>%
      select("mu", "powerNew", "powerMeta", "powerThreshold1", "powerThreshold2","powerThreshold3") %>%
      rename( "Replication Analysis" = "powerNew" , "Meta Analysis" = "powerMeta" , "Threshold 10-8" = "powerThreshold1" , "Threshold 10-5" = "powerThreshold2",  "Threshold 10-4" = "powerThreshold3") %>%
      melt(id.vars = "mu") 
    colnames(plotDatPowerMu) <- c("Mean","Method","Power")
    
    plot1 <- ggplot(data=plotDatPowerMu, aes(x=Mean, y=Power, group=Method)) +
      geom_line(aes(color=Method,linetype=Method),lwd=1.2) +
      # scale_color_discrete(labels = c(expression(paste("Replication Analysis")),
      #                                 expression(paste("Meta Analysis")),
      #                                 expression(paste("P <",10^-8)),
      #                                 expression(paste("P <",10^-6)),
      #                                 expression(paste("P <",10^-5)))) + 
      ylab("Power") +
      xlab("Mean") +
      ylim(c(0, 1)) + xlim(c(2, 6)) +
      #geom_hline(yintercept = 0.1, linetype=2, color="grey") +
      scale_linetype_manual(values=1:6,labels = c(expression(paste("Replication Analysis")),
                                                  expression(paste("Meta Analysis")),
                                                  expression(paste("P <",10^-8)),
                                                  expression(paste("P <",10^-6)),
                                                  expression(paste("P <",10^-5)))) +
      scale_color_manual(values=mycols,labels = c(expression(paste("Replication Analysis")),
                                                  expression(paste("Meta Analysis")),
                                                  expression(paste("P <",10^-8)),
                                                  expression(paste("P <",10^-6)),
                                                  expression(paste("P <",10^-5)))) +
      theme_cowplot() +
      theme(axis.title = element_text(size=14), axis.text = element_text(size=14)) +
      theme(legend.title = element_text(size=14), legend.text = element_text(size=14)) +
      ggtitle(paste("Proportion of Causal Variants is", format(pCausal, scientific=F))) 
    
  } else if (type == "FDP") {
    
    #FDP plot for different mean
    plotDatFDPMu <- final_result %>%
      select("mu", "fdpNew", "fdpMeta", "fdpThreshold1", "fdpThreshold2","fdpThreshold3") %>%
      rename( "Replication Analysis" = "fdpNew" , "Meta Analysis" = "fdpMeta" , "Threshold 10-8" = "fdpThreshold1" , "Threshold 10-5" = "fdpThreshold2",  "Threshold 10-4" = "fdpThreshold3") %>%
      melt(id.vars = "mu") %>%
      replace(is.na(.), 0)
    colnames(plotDatFDPMu) <- c("Mean","Method","FDP")
    
    plot1 <- ggplot(data=plotDatFDPMu, aes(x=Mean, y=FDP, group=Method)) +
      geom_line(aes(color=Method,linetype=Method),lwd=1.2) +
      # scale_color_discrete(labels = c(expression(paste("Replication Analysis")),
      #                                 expression(paste("Meta Analysis")),
      #                                 expression(paste("P <",10^-8)),
      #                                 expression(paste("P <",10^-6)),
      #                                 expression(paste("P <",10^-5)))) + 
      ylab("FDP") +
      xlab("Mean") +
      ylim(c(0, 1)) + xlim(c(2, 6)) +
      geom_hline(yintercept = 0.1, linetype=2, color="grey") +
      scale_linetype_manual(values=1:6,labels = c(expression(paste("Replication Analysis")),
                                                  expression(paste("Meta Analysis")),
                                                  expression(paste("P <",10^-8)),
                                                  expression(paste("P <",10^-6)),
                                                  expression(paste("P <",10^-5)))) +
      scale_color_manual(values=mycols,labels = c(expression(paste("Replication Analysis")),
                                                  expression(paste("Meta Analysis")),
                                                  expression(paste("P <",10^-8)),
                                                  expression(paste("P <",10^-6)),
                                                  expression(paste("P <",10^-5)))) +
      theme_cowplot() +
      theme(axis.title = element_text(size=14), axis.text = element_text(size=14)) +
      theme(legend.title = element_text(size=14), legend.text = element_text(size=14)) +
      theme(legend.key.size = unit(3,"line")) +
      ggtitle(paste("Proportion of Causal Variants is", format(pCausal, scientific=F))) 
    
    
  } 
  return(plot1)
}



fig1a <- plot_func(data = sim_results,pCausal = 0.0002, type = "Power")
fig1b <- plot_func(data = sim_results,pCausal = 0.0002, type = "FDP")

fig1c <- plot_func(data = sim_results,pCausal = 0.002, type = "Power")
fig1d <- plot_func(data = sim_results,pCausal = 0.002, type = "FDP")

fig1e <- plot_func(data = sim_results,pCausal = 0.01, type = "Power")
fig1f <- plot_func(data = sim_results,pCausal = 0.01, type = "FDP")

mergePlot <- ggarrange(fig1a, fig1b,fig1c,fig1d,fig1e,fig1f , ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
ggsave(paste0(outputDir_Fig,"/Fig3.png"), plot = mergePlot, width = 13,
       height = 9)
