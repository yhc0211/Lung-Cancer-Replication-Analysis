# Make all simulation plots

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/plot_simulation.R")

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
library(scales)



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
mycols[7] <- "purple"


sim_results_two <- NULL

for (Snum in 1:15) {
  for (aID in 1:100) {
    
    file_path <- paste0(outputDir, "/rep_two_sim_Snum", Snum, "_aID", aID, ".txt")
    
    if (!file.exists(file_path)) {
      message("Skipping missing file: ", file_path)
      next
    }
    
    tempResults <- fread(file_path) %>%
      dplyr::select(
        nCausal, propCausal, mu,
        powerNew, powerMeta, powerMetaBH, powerRepfdr, powerMamba,
        powerThreshold1, powerThreshold2, powerThreshold3,
        fdpNew, fdpMeta, fdpMetaBH, fdpRepfdr, fdpMamba,
        fdpThreshold1, fdpThreshold2, fdpThreshold3,
        inconNew, inconRepfdr,
        fdpKernel, powerKernel, nRejKernel, inconKernel
      )
    
    sim_results_two <- rbind(sim_results_two, tempResults)
  }
}


sim_results_three <- NULL

for (Snum in 1:15) {
  for (aID in 1:100) {
    
    file_path <- paste0(outputDir, "/rep_three_sim_Snum", Snum, "_aID", aID, ".txt")
    
    if (!file.exists(file_path)) {
      message("Skipping missing file: ", file_path)
      next
    }
    
    tempResults <- fread(file_path) %>%
      dplyr::select(
        nCausal, propCausal, mu,
        powerNew, powerMeta, powerMetaBH, powerRepfdr, powerMamba,
        powerThreshold1, powerThreshold2, powerThreshold3,
        fdpNew, fdpMeta, fdpMetaBH, fdpRepfdr, fdpMamba,
        fdpThreshold1, fdpThreshold2, fdpThreshold3,
        inconNew, inconRepfdr,
        fdpKernel, powerKernel, nRejKernel, inconKernel
      )
    
    sim_results_three <- rbind(sim_results_three, tempResults)
  }
}





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
    "P-value threshold < 1e-8",
    "P-value threshold < 1e-6",
    "P-value threshold < 1e-5",
    "repfdr",
    "kernel"
  )
  
  color_map <- c(
    "Replication Analysis"      = mycols[1],
    "Meta-analysis"             = mycols[2],
    "P-value threshold < 1e-8"  = mycols[3],
    "P-value threshold < 1e-6"  = mycols[4],
    "P-value threshold < 1e-5"  = mycols[5],
    "repfdr"                    = mycols[6],
    "kernel"                    = mycols[7]
  )
  
  linetype_map <- c(
    "Replication Analysis"      = 1,
    "Meta-analysis"             = 2,
    "P-value threshold < 1e-8"  = 3,
    "P-value threshold < 1e-6"  = 4,
    "P-value threshold < 1e-5"  = 5,
    "repfdr"                    = 6,
    "kernel"                    = 7
  )
  
  legend_labels <- c(
    "Replication Analysis",
    "Meta-analysis",
    expression("P-value threshold < " * 10^{-8}),
    expression("P-value threshold < " * 10^{-6}),
    expression("P-value threshold < " * 10^{-5}),
    "repfdr",
    "kernel"
  )
  
  if (type == "Power") {
    
    plotDat <- final_result %>%
      dplyr::select(
        mu, powerNew, powerNeff, 
        powerThreshold1, powerThreshold2, powerThreshold3, powerRepfdr, powerKernel
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = powerNew,
        "Meta-analysis"             = powerMeta,
        "P-value threshold < 1e-8"  = powerThreshold1,
        "P-value threshold < 1e-6"  = powerThreshold2,
        "P-value threshold < 1e-5"  = powerThreshold3,
        "repfdr"                    = powerRepfdr,
        "kernel"                    = powerKernel
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      )
    
    ylab_text <- "Power"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 1) +
      ylab(ylab_text) +
      xlab("Mean Z-score under alternative") +
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
      ) 
    
  } else if (type == "FDP") {
    
    plotDat <- final_result %>%
      dplyr::select(
        mu, fdpNew, fdpNeff, 
        fdpThreshold1, fdpThreshold2, fdpThreshold3, fdpRepfdr, fdpKernel
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = fdpNew,
        "Meta-analysis"             = fdpMeta,
        "P-value threshold < 1e-8"  = fdpThreshold1,
        "P-value threshold < 1e-6"  = fdpThreshold2,
        "P-value threshold < 1e-5"  = fdpThreshold3,
        "repfdr"                    = fdpRepfdr,
        "kernel"                    = fdpKernel
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      ) %>%
      tidyr::replace_na(list(Value = 0))
    
    ylab_text <- "FDR"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 1) +
      ylab(ylab_text) +
      xlab("Mean Z-score under alternative") +
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
      ) 
    
  } else {
    stop("type must be either 'Power' or 'FDP'")
  }
  
  return(plot1)
}

fig1a <- plot_func(data = sim_results_two,pCausal = 0.002, type = "FDP")
fig1b <- plot_func(data = sim_results_two,pCausal = 0.002, type = "Power")
fig1c <- plot_func(data = sim_results_three,pCausal = 0.002, type = "FDP")
fig1d <- plot_func(data = sim_results_three,pCausal = 0.002, type = "Power")

mergePlot <- ggarrange(fig1a, fig1b,fig1c,fig1d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom",labels = c("A","B","C","D"))
ggsave(paste0(outputDir_Fig,"/simulation_main.pdf"), plot = mergePlot, width = 11, height = 9)

supple_plot_func <- function(data, pCausal, type) {
  
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
    "P-value threshold < 1e-8",
    "P-value threshold < 1e-6",
    "P-value threshold < 1e-5",
    "repfdr",
    "kernel"
  )
  
  color_map <- c(
    "Replication Analysis"      = mycols[1],
    "Meta-analysis"             = mycols[2],
    "P-value threshold < 1e-8"  = mycols[3],
    "P-value threshold < 1e-6"  = mycols[4],
    "P-value threshold < 1e-5"  = mycols[5],
    "repfdr"                    = mycols[6],
    "kernel"                    = mycols[7]
  )
  
  linetype_map <- c(
    "Replication Analysis"      = 1,
    "Meta-analysis"             = 2,
    "P-value threshold < 1e-8"  = 3,
    "P-value threshold < 1e-6"  = 4,
    "P-value threshold < 1e-5"  = 5,
    "repfdr"                    = 6,
    "kernel"                    = 7
  )
  
  legend_labels <- c(
    "Replication Analysis",
    "Meta-analysis",
    expression("P-value threshold < " * 10^{-8}),
    expression("P-value threshold < " * 10^{-6}),
    expression("P-value threshold < " * 10^{-5}),
    "repfdr",
    "kernel"
  )
  
  if (type == "Power") {
    
    plotDat <- final_result %>%
      dplyr::select(
        mu, powerNew, powerNeff, 
        powerThreshold1, powerThreshold2, powerThreshold3, powerRepfdr, powerKernel
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = powerNew,
        "Meta-analysis"             = powerMeta,
        "P-value threshold < 1e-8"  = powerThreshold1,
        "P-value threshold < 1e-6"  = powerThreshold2,
        "P-value threshold < 1e-5"  = powerThreshold3,
        "repfdr"                    = powerRepfdr,
        "kernel"                    = powerKernel
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      )
    
    ylab_text <- "Power"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 1) +
      ylab(ylab_text) +
      xlab("Mean Z-score under alternative") +
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
      ) 
  } else if (type == "FDP") {
    
    plotDat <- final_result %>%
      dplyr::select(
        mu, fdpNew, fdpNeff, 
        fdpThreshold1, fdpThreshold2, fdpThreshold3, fdpRepfdr, fdpKernel
      ) %>%
      dplyr::rename(
        "Replication Analysis"      = fdpNew,
        "Meta-analysis"             = fdpMeta,
        "P-value threshold < 1e-8"  = fdpThreshold1,
        "P-value threshold < 1e-6"  = fdpThreshold2,
        "P-value threshold < 1e-5"  = fdpThreshold3,
        "repfdr"                    = fdpRepfdr,
        "kernel"                    = fdpKernel
      ) %>%
      tidyr::pivot_longer(
        cols = -mu,
        names_to = "Method",
        values_to = "Value"
      ) %>%
      tidyr::replace_na(list(Value = 0))
    
    ylab_text <- "FDR"
    
    plot1 <- ggplot(plotDat, aes(x = mu, y = Value, color = Method, linetype = Method)) +
      geom_line(linewidth = 1) +
      ylab(ylab_text) +
      xlab("Mean Z-score under alternative") +
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
          "Proportion of Causal Variants:",
          format(pCausal * 100, scientific = FALSE),
          "%"
        )
      )
    
  } else {
    stop("type must be either 'Power' or 'FDP'")
  }
  
  return(plot1)
}

supple_fig1a <- supple_plot_func(data = sim_results_two,pCausal = 0.0002, type = "FDP")
supple_fig1b <- supple_plot_func(data = sim_results_two,pCausal = 0.0002, type = "Power")
supple_fig1c <- supple_plot_func(data = sim_results_two,pCausal = 0.01, type = "FDP")
supple_fig1d <- supple_plot_func(data = sim_results_two,pCausal = 0.01, type = "Power")
supple_mergePlot <- ggarrange(supple_fig1a, supple_fig1b,supple_fig1c,supple_fig1d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom",labels = c("A","B","C","D"))
ggsave(paste0(outputDir_Fig,"/simulation_supple_two_cohort.pdf"), plot = supple_mergePlot, width = 11, height = 9)


supple_fig1a <- supple_plot_func(data = sim_results_three,pCausal = 0.0002, type = "FDP")
supple_fig1b <- supple_plot_func(data = sim_results_three,pCausal = 0.0002, type = "Power")
supple_fig1c <- supple_plot_func(data = sim_results_three,pCausal = 0.01, type = "FDP")
supple_fig1d <- supple_plot_func(data = sim_results_three,pCausal = 0.01, type = "Power")
supple_mergePlot <- ggarrange(supple_fig1a, supple_fig1b,supple_fig1c,supple_fig1d, ncol=2, nrow=2, common.legend = TRUE, legend="bottom",labels = c("A","B","C","D"))
ggsave(paste0(outputDir_Fig,"/simulation_supple_three_cohort.pdf"), plot = supple_mergePlot, width = 11, height = 9)






sim_results_two_mean <- sim_results_two %>%
  
  group_by(propCausal, mu) %>%
  
  summarise(
    
    across(
      
      .cols = where(is.numeric),
      
      .fns  = ~ mean(.x, na.rm = TRUE)
      
    ),
    
    .groups = "drop"
    
  )




sim_results_three_mean <- sim_results_three %>%
  
  group_by(propCausal, mu) %>%
  
  summarise(
    
    across(
      
      .cols = where(is.numeric),
      
      .fns  = ~ mean(.x, na.rm = TRUE)
      
    ),
    
    .groups = "drop"
    
  )

write.csv(
  
  sim_results_two_mean,
  
  file = paste0(outputDir_Fig,"/sim_results_two_mean.csv"),
  
  row.names = FALSE
  
)


write.csv(
  
  sim_results_three_mean,
  
  file = paste0(outputDir_Fig,"/sim_results_three_mean.csv"),
  
  row.names = FALSE
  
)



# ============================================================
# Create incongruous-ranking plot data
# ============================================================

# ============================================================
# Build plot_dat_two / plot_dat_three for incongruous plot
# ============================================================

make_incon_plot_dat <- function(sim_results) {
  sim_results %>%
    dplyr::group_by(mu, propCausal) %>%
    dplyr::summarise(
      `Replication Analysis` = mean(inconNew,    na.rm = TRUE),
      `repfdr`               = mean(inconRepfdr, na.rm = TRUE),
      `kernel`               = mean(inconKernel, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols      = c(`Replication Analysis`, `repfdr`, `kernel`),
      names_to  = "method",
      values_to = "mean_incon"
    ) %>%
    dplyr::mutate(
      method = factor(
        method,
        levels = c("kernel", "repfdr", "Replication Analysis")  # y-axis bottom -> top
      ),
      propCausal_lab = factor(
        paste0(
          "Proportion of causal variants: ",
          format(propCausal * 100, scientific = FALSE, trim = TRUE),
          "%"
        ),
        levels = paste0(
          "Proportion of causal variants: ",
          format(sort(unique(propCausal)) * 100, scientific = FALSE, trim = TRUE),
          "%"
        )
      )
    )
}

plot_dat_two   <- make_incon_plot_dat(sim_results_two)
plot_dat_three <- make_incon_plot_dat(sim_results_three)


write.csv(
  
  plot_dat_two,
  
  file = paste0(outputDir_Fig,"/sim_results_two_incon.csv"),
  
  row.names = FALSE
  
)

write.csv(
  
  plot_dat_three,
  
  file = paste0(outputDir_Fig,"/sim_results_three_incon.csv"),
  
  row.names = FALSE
  
)



plot_dat_two <- plot_dat_two %>%
  mutate(
    log_mean_incon = log10(mean_incon + 1)
  )


p_incon_two <- ggplot(
  plot_dat_two,
  aes(
    x = mu,
    y = method,
    size = log_mean_incon
  )
) +
  geom_point(alpha = 0.85, color = "#2E86AB") +
  facet_wrap(
    ~ propCausal_lab,
    nrow = 1,
    strip.position = "top"
  ) +
  scale_size_continuous(
    range = c(1.8, 14),
    breaks = log10(c(0, 10, 100, 1000) + 1),
    labels = c("0", "10", "100", "1000"),
    name = "Mean incongruous"
  ) +
  labs(
    x = "Mean Z-score under alternative",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 14),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  )
p_incon_two

plot_dat_three <- plot_dat_three %>%
  mutate(
    log_mean_incon = log10(mean_incon + 1)
  )

p_incon_three <- ggplot(
  plot_dat_three,
  aes(
    x = mu,
    y = method,
    size = log_mean_incon
  )
) +
  geom_point(alpha = 0.85, color = "#2E86AB") +
  facet_wrap(
    ~ propCausal_lab,
    nrow = 1,
    strip.position = "top"
  ) +
  scale_size_continuous(
    range = c(1.8, 14),
    breaks = log10(c(0, 10, 100, 1000) + 1),
    labels = c("0", "10", "100", "1000"),
    name = "Mean incongruous"
  ) +
  labs(
    x = "Mean Z-score under alternative",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 18),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 14),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  )
p_incon_three

incon_mergePlot <- ggarrange(
  p_incon_two,
  p_incon_three,
  ncol = 1,
  nrow = 2,
  common.legend = FALSE,
  legend = "right",
  labels = c("B", "C")
)




library(ggplot2)
library(dplyr)
library(tibble)
library(grid)
library(patchwork)



# ============================================================
# A. Scatter plot: Z-score evidence
# ============================================================

library(tidyverse)
library(grid)

scatter_dat <- tibble(
  snp  = c("rs397640", "rs1265053"),
  z1   = c(-5.83, -5.09),
  z2   = c(-3.07, -2.70),
  lfdr = c(0.13, 0.08)
)

p_scatter <- ggplot(scatter_dat, aes(x = z1, y = z2)) +
  geom_point(
    aes(shape = snp, fill = snp),
    size = 6,
    stroke = 1.1,
    color = "black"
  ) +
  
  # Arrow from rs1265053 to rs397640
  geom_segment(
    x = -5.09, y = -2.70,
    xend = -5.83, yend = -3.07,
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.9,
    color = "black"
  ) +
  
  # Annotation for rs1265053
  annotate(
    "text",
    x = -4.95, y = -2.55,
    label = "SNP 2\nrs1265053\n(-5.09, -2.70)",
    hjust = 0,
    vjust = 0,
    size = 6
  ) +
  
  # Annotation for rs397640
  annotate(
    "text",
    x = -5.95, y = -3.20,
    label = "SNP 1\nrs397640\n(-5.83, -3.07)",
    hjust = 1,
    vjust = 1,
    size = 6
  ) +
  
  scale_shape_manual(
    values = c("rs397640" = 21, "rs1265053" = 23)
  ) +
  scale_fill_manual(
    values = c("rs397640" = "#2C7FB8", "rs1265053" = "#F28E2B")
  ) +
  labs(
    x = "ILCCO Z-score",
    y = "MVP Z-score"
  ) +
  coord_cartesian(
    xlim = c(-6.20, -4.50),
    ylim = c(-3.40, -2.40),
    clip = "off"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.margin = margin(15, 80, 15, 60)
  )

p_scatter
# ============================================================
# B. Rank crossing plot
# ============================================================
rank_dat <- tibble(
  snp = c("rs397640", "rs1265053"),
  x1  = 1,
  x2  = 2,
  y1  = c(1, 2),   # rank by Z-score evidence
  y2  = c(2, 1),   # rank by lfdr
  z1  = c(-5.83, -5.09),
  z2  = c(-3.07, -2.70),
  lfdr = c(0.13, 0.08)
)

p_rank <- ggplot(rank_dat) +
  geom_segment(
    aes(x = x1, xend = x2, y = y1, yend = y2, color = snp),
    linewidth = 1.3
  ) +
  geom_point(
    aes(x = x1, y = y1, shape = snp, fill = snp, color = snp),
    size = 5.5,
    stroke = 1.1
  ) +
  geom_point(
    aes(x = x2, y = y2, shape = snp, fill = snp, color = snp),
    size = 5.5,
    stroke = 1.1
  ) +
  
  # Left side: Z-score evidence ranking
  annotate(
    "text",
    x = 0.82, y = 1,
    label = "rank 1\nrs397640\n(-5.83, -3.07)",
    hjust = 1,
    fontface = "bold",
    size = 7,
    color = "#2C7FB8"
  ) +
  annotate(
    "text",
    x = 0.82, y = 2,
    label = "rank 2\nrs1265053\n(-5.09, -2.70)",
    hjust = 1,
    fontface = "bold",
    size = 7,
    color = "#F28E2B"
  ) +
  
  # Right side: lfdr ranking
  annotate(
    "text",
    x = 2.18, y = 1,
    label = "rank 1\nrs1265053\nlfdr = 0.08",
    hjust = 0,
    fontface = "bold",
    size = 7,
    color = "#F28E2B"
  ) +
  annotate(
    "text",
    x = 2.18, y = 2,
    label = "rank 2\nrs397640\nlfdr = 0.13",
    hjust = 0,
    fontface = "bold",
    size = 7,
    color = "#2C7FB8"
  ) +
  
  annotate(
    "label",
    x = 1.5, y = 1.5,
    label = "Bayesian lfdr methods",
    size = 7,
    label.size = 0.3,
    fill = "white"
  ) +
  
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Z-score evidence\nrank", "lfdr\nrank"),
    position = "top",
    limits = c(0.45, 2.55)
  ) +
  scale_y_reverse(
    breaks = c(1, 2),
    labels = c("rank 1", "rank 2"),
    limits = c(2.35, 0.65)
  ) +
  scale_shape_manual(
    values = c("rs397640" = 21, "rs1265053" = 23)
  ) +
  scale_fill_manual(
    values = c("rs397640" = "#2C7FB8", "rs1265053" = "#F28E2B")
  ) +
  scale_color_manual(
    values = c("rs397640" = "#2C7FB8", "rs1265053" = "#F28E2B")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x.top = element_text(
      size = 18,
      face = "bold",
      margin = margin(b = 8)
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 20, 5, 20)
  )

p_rank
# ============================================================
# C. Combine A and B with controlled width and height
# ============================================================

top_row <- p_scatter + p_rank +
  plot_layout(
    widths = c(1.0, 1.45)   # adjust width of A and B here
  )

# ============================================================
# D. Combine with p_incon_two and p_incon_three
# ============================================================
# Assumes you already created:
# p_incon_two
# p_incon_three

figure2 <- 
  top_row /
  p_incon_two /
  p_incon_three +
  plot_layout(
    heights = c(1.8, 1, 1)  # adjust height of top A/B row here
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold")
    )
  )

figure2

ggsave(paste0(outputDir_Fig,"/simulation_incon_main.pdf"), plot = figure2, width = 25, height = 16)
