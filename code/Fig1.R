# Make Fig1 - power and fdp plot for different mean and number of causal

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/Fig1.R")

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
outputDir <- here::here("Lung", "output_revised")
outputDir_Fig <- here::here("Lung", "Fig_Tab")
dataDir <- here::here("Data")


# define colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

mycols <- gg_color_hue(8)
mycols[4] <- "green"
mycols[5] <- "blue"
mycols[7] <- "black"
mycols[8] <- "purple"

sim_results <- fread(paste0(outputDir, "/rep_two_sim_Snum", 1, "_aID", 1, ".txt")) %>%
  dplyr::select(
    nCausal, propCausal, mu,
    powerNew, powerMeta, powerMetaBH, powerRepfdr, powerMamba,
    powerThreshold1, powerThreshold2, powerThreshold3,
    fdpNew, fdpMeta, fdpMetaBH, fdpRepfdr, fdpMamba,
    fdpThreshold1, fdpThreshold2, fdpThreshold3,
    inconNew, inconRepfdr
  )

for (Snum in 1:15) {
  for (aID in 1:100) {
    tempResults <- fread(paste0(outputDir, "/rep_two_sim_Snum", Snum, "_aID", aID, ".txt")) %>%
      dplyr::select(
        nCausal, propCausal, mu,
        powerNew, powerMeta, powerMetaBH, powerRepfdr, powerMamba,
        powerThreshold1, powerThreshold2, powerThreshold3,
        fdpNew, fdpMeta, fdpMetaBH, fdpRepfdr, fdpMamba,
        fdpThreshold1, fdpThreshold2, fdpThreshold3,
        inconNew, inconRepfdr
      )
    sim_results <- rbind(sim_results, tempResults)
  }
}

sim_results <- sim_results[-1, ]

#----------------------------------------------------------------------------------------#
#power plot for different mu 

plot_func <- function(data, pCausal, type) {
  
  sim_results <- data %>%
    dplyr::filter(propCausal == pCausal)
  
  final_result <- sim_results %>%
    dplyr::group_by(mu) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  method_order <- c(
    "Replication Analysis",
    "Meta-analysis",
    "Meta-analysis (BH)",
    "P-value threshold < 1e-8",
    "P-value threshold < 1e-6",
    "P-value threshold < 1e-5",
    "repfdr"
  )
  
  color_map <- c(
    "Replication Analysis"      = mycols[1],
    "Meta-analysis"             = mycols[2],
    "Meta-analysis (BH)"        = mycols[3],
    "P-value threshold < 1e-8"  = mycols[4],
    "P-value threshold < 1e-6"  = mycols[5],
    "P-value threshold < 1e-5"  = mycols[6],
    "repfdr"                    = mycols[7]
  )
  
  linetype_map <- c(
    "Replication Analysis"      = 1,
    "Meta-analysis"             = 2,
    "Meta-analysis (BH)"        = 3,
    "P-value threshold < 1e-8"  = 4,
    "P-value threshold < 1e-6"  = 5,
    "P-value threshold < 1e-5"  = 6,
    "repfdr"                    = 7
  )
  
  legend_labels <- c(
    "Replication Analysis",
    "Meta-analysis",
    "Meta-analysis (BH)",
    expression("P-value threshold < " * 10^{-8}),
    expression("P-value threshold < " * 10^{-6}),
    expression("P-value threshold < " * 10^{-5}),
    "repfdr"
  )
  
  if (type == "Power") {
    
    plotDat <- final_result %>%
      dplyr::select(
        mu, powerNew, powerMeta, powerMetaBH,
        powerThreshold1, powerThreshold2, powerThreshold3, powerRepfdr
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = powerNew,
        "Meta-analysis"             = powerMeta,
        "Meta-analysis (BH)"        = powerMetaBH,
        "P-value threshold < 1e-8"  = powerThreshold1,
        "P-value threshold < 1e-6"  = powerThreshold2,
        "P-value threshold < 1e-5"  = powerThreshold3,
        "repfdr"                    = powerRepfdr
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      )
    
    ylab_text <- "Power"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 0.5) +
      ylab(ylab_text) +
      xlab("Mean Effect Size") +
      ylim(0, 1) +
      xlim(2, 6) +
      scale_color_manual(
        values = color_map,
        breaks = method_order,
        labels = legend_labels
      ) +
      scale_linetype_manual(
        values = linetype_map,
        breaks = method_order,
        labels = legend_labels
      ) +
      theme_cowplot() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "line")
      ) +
      ggtitle(
        paste(
          "Proportion of Causal Variants is",
          format(pCausal * 100, scientific = FALSE),
          "%"
        )
      )
    
  } else if (type == "FDP") {
    
    plotDat <- final_result %>%
      dplyr::select(
        mu, fdpNew, fdpMeta, fdpMetaBH,
        fdpThreshold1, fdpThreshold2, fdpThreshold3, fdpRepfdr
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = fdpNew,
        "Meta-analysis"             = fdpMeta,
        "Meta-analysis (BH)"        = fdpMetaBH,
        "P-value threshold < 1e-8"  = fdpThreshold1,
        "P-value threshold < 1e-6"  = fdpThreshold2,
        "P-value threshold < 1e-5"  = fdpThreshold3,
        "repfdr"                    = fdpRepfdr
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      ) %>%
      tidyr::replace_na(list(Value = 0))
    
    ylab_text <- "FDR"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 0.5) +
      ylab(ylab_text) +
      xlab("Mean Effect Size") +
      ylim(0, 1) +
      xlim(2, 6) +
      geom_hline(yintercept = 0.1, linetype = 2, color = "grey") +
      scale_color_manual(
        values = color_map,
        breaks = method_order,
        labels = legend_labels
      ) +
      scale_linetype_manual(
        values = linetype_map,
        breaks = method_order,
        labels = legend_labels
      ) +
      theme_cowplot() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "line")
      ) +
      ggtitle(
        paste(
          "Proportion of Causal Variants is",
          format(pCausal * 100, scientific = FALSE),
          "%"
        )
      )
    
  } else {
    stop("type must be either 'Power' or 'FDP'")
  }
  
  return(plot1)
}
fig1a <- plot_func(data = sim_results,pCausal = 0.0002, type = "FDP")
fig1b <- plot_func(data = sim_results,pCausal = 0.0002, type = "Power")

fig1c <- plot_func(data = sim_results,pCausal = 0.002, type = "FDP")
fig1d <- plot_func(data = sim_results,pCausal = 0.002, type = "Power")

fig1e <- plot_func(data = sim_results,pCausal = 0.01, type = "FDP")
fig1f <- plot_func(data = sim_results,pCausal = 0.01, type = "Power")


mergePlot <- ggarrange(fig1a, fig1b,fig1c,fig1d,fig1e,fig1f , ncol=2, nrow=3, common.legend = TRUE, legend="bottom",labels = c("A","B","C","D","E","F"))
ggsave(paste0(outputDir_Fig,"/Fig1_revised.pdf"), plot = mergePlot, width = 14,
       height = 9)

plot_incon_func <- function(data, pCausal) {
  
  sim_results <- data %>%
    dplyr::filter(propCausal == pCausal)
  
  final_result <- sim_results %>%
    dplyr::group_by(mu) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  method_order <- c(
    "Replication Analysis",
    "Meta-analysis",
    "Meta-analysis (BH)",
    "P-value threshold < 1e-8",
    "P-value threshold < 1e-6",
    "P-value threshold < 1e-5",
    "repfdr"
  )
  
  color_map <- c(
    "Replication Analysis"      = mycols[1],
    "Meta-analysis"             = mycols[2],
    "Meta-analysis (BH)"        = mycols[3],
    "P-value threshold < 1e-8"  = mycols[4],
    "P-value threshold < 1e-6"  = mycols[5],
    "P-value threshold < 1e-5"  = mycols[6],
    "repfdr"                    = mycols[7]
  )
  
  linetype_map <- c(
    "Replication Analysis"      = 1,
    "Meta-analysis"             = 2,
    "Meta-analysis (BH)"        = 3,
    "P-value threshold < 1e-8"  = 4,
    "P-value threshold < 1e-6"  = 5,
    "P-value threshold < 1e-5"  = 6,
    "repfdr"                    = 7
  )
  
  legend_labels <- c(
    "Replication Analysis",
    "Meta-analysis",
    "Meta-analysis (BH)",
    expression("P-value threshold < " * 10^{-8}),
    expression("P-value threshold < " * 10^{-6}),
    expression("P-value threshold < " * 10^{-5}),
    "repfdr"
  )
  

    plotDat <- final_result %>%
      dplyr::select(
        mu, inconNew, inconRepfdr
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = inconNew,
        "repfdr"                    = inconRepfdr
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      )
    
    ylab_text <- "Num Incongruous"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 0.5) +
      ylab(ylab_text) +
      xlab("Mean Effect Size") +
      ylim(0, 10) +
      xlim(2, 6) +
      scale_color_manual(
        values = color_map,
        breaks = method_order,
        labels = legend_labels
      ) +
      scale_linetype_manual(
        values = linetype_map,
        breaks = method_order,
        labels = legend_labels
      ) +
      theme_cowplot() +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "line")
      ) +
      ggtitle(
        paste(
          "Proportion of Causal Variants is",
          format(pCausal * 100, scientific = FALSE),
          "%"
        )
      )
    
 
}
figa <- plot_incon_func(data = sim_results,pCausal = 0.0002)
figb <- plot_incon_func(data = sim_results,pCausal = 0.002)
figc <- plot_incon_func(data = sim_results,pCausal = 0.01)
mergePlot_incon <- ggarrange(figa, figb,figc, ncol=1, nrow=3, common.legend = TRUE, legend="bottom",labels = c("A","B","C"))
ggsave(paste0(outputDir_Fig,"/Fig1_incon_revised.pdf"), plot = mergePlot_incon, width = 14,
       height = 9)
