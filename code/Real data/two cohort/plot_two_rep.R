# Make Two-way Replication Analysis Results
#
# Pairs (from make_para):
#   aID = 1 -> ILCCO_UKB
#   aID = 2 -> ILCCO_MVP
#   aID = 3 -> UKB_MVP
#


# ============================================================
# Setup
# ============================================================

setwd("/rsrch8/home/epi/ychang11/csmGmm_reproduce/Lung")
here::i_am("Lung/Fig5.R")

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(xtable)
library(devtools)
library(wesanderson)
library(here)
library(manhattan)
library(ggpubr)
library(tidyr)
library(forcats)
library(purrr)
library(ggrastr)

# ============================================================
# Read aID from command line / array job
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

aID <- as.integer(args[1])


# Pair label for messages / filenames
pair_labels <- c("ILCCO_UKB", "ILCCO_MVP", "UKB_MVP")
pair_name <- pair_labels[aID]
cat("Running aID =", aID, " (pair:", pair_name, ")\n")

# ============================================================
# Source supporting code
# ============================================================

codePath <- here::here("SupportingCode")
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# ============================================================
# Directories
# ============================================================

outputDir <- here::here("Lung", "output")
outputDir_Fig <- here::here("Lung", "Fig_Tab")
dataDir <- here::here("Data")

dir.create(outputDir_Fig, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper functions
# ============================================================

calc_neff <- function(n_case, n_ctrl) {
  4 / (1 / n_case + 1 / n_ctrl)
}

make_para <- function(analysis = c("overall", "LUSC", "LUAD")) {
  
  analysis <- match.arg(analysis)
  
  if (analysis == "overall") {
    N_ILCCO <- calc_neff(29266, 56450)
    N_MVP   <- calc_neff(10398, 62708)
    N_UKB   <- calc_neff(2404, 330018)
  }
  
  if (analysis == "LUSC") {
    N_ILCCO <- calc_neff(7426, 56450)
    N_MVP   <- calc_neff(1475, 62708)
    N_UKB   <- calc_neff(489, 330018)
  }
  
  if (analysis == "LUAD") {
    N_ILCCO <- calc_neff(11273, 56450)
    N_MVP   <- calc_neff(2019, 62708)
    N_UKB   <- calc_neff(862, 330018)
  }
  
  para <- data.frame(
    pair   = c("ILCCO_UKB", "ILCCO_MVP", "UKB_MVP"),
    p_1    = c("pILCCO", "pILCCO", "pUKB"),
    p_2    = c("pUKB",   "pMVP",   "pMVP"),
    Z_1    = c("zILCCO", "zILCCO", "zUKB"),
    Z_2    = c("zUKB",   "zMVP",   "zMVP"),
    Neff_1 = c(N_ILCCO, N_ILCCO, N_UKB),
    Neff_2 = c(N_UKB,   N_MVP,   N_MVP)
  )
  
  return(para)
}

run_meta_pair <- function(three_Z, para, aID, p_cutoff = 5e-8) {
  
  p1_col <- para$p_1[aID]
  p2_col <- para$p_2[aID]
  z1_col <- para$Z_1[aID]
  z2_col <- para$Z_2[aID]
  
  w1 <- sqrt(para$Neff_1[aID])
  w2 <- sqrt(para$Neff_2[aID])
  
  out <- three_Z %>%
    transmute(
      chrpos = chrpos,
      RS = RS,
      P1 = .data[[p1_col]],
      Z1 = .data[[z1_col]],
      P2 = .data[[p2_col]],
      Z2 = .data[[z2_col]]
    ) %>%
    filter(!is.na(Z1), !is.na(Z2)) %>%
    mutate(
      Z_meta = (Z1 * w1 + Z2 * w2) / sqrt(w1^2 + w2^2),
      P_meta = 2 * pnorm(abs(Z_meta), lower.tail = FALSE),
      method = "Meta-analysis",
      pair = para$pair[aID]
    )
  
  out_sig <- out %>%
    filter(P_meta < p_cutoff) %>%
    select(chrpos, RS, pair, method, Z_meta, P_meta)
  
  return(list(
    all = out,
    sig = out_sig
  ))
}

plotmanhattan <- function(gwas,
                          build = c("hg18", "hg19", "hg38"),
                          colValues,
                          shapeValues,
                          ylimits,
                          legName,
                          title) {
  
  data <- gwas
  build <- match.arg(build)
  
  cat_levels <- c("Both", "Meta-analysis", "Replication Analysis", "Neither")
  data$cat <- factor(data$cat, levels = cat_levels)
  
  data <- add_cumulative_pos(data, build)
  chrom_lengths <- get_chrom_lengths(build)
  xmax <- get_total_length(chrom_lengths)
  x_breaks <- get_x_breaks(chrom_lengths)
  
  names(x_breaks)[c(11, 13, 14, 16, 17, 18, 20, 21)] <- ""
  
  returnPlot <- ggplot(
    data,
    aes(x = cumulative_pos, y = y, color = cat, shape = cat)
  ) +
    ggrastr::rasterise(
      geom_point(size = 0.8, alpha = 0.85),
      dpi = 300
    ) +
    xlab("Chromosome") +
    ylab(expression(paste("-", log[10], "(FDR)"))) +
    scale_color_manual(
      name = legName, values = colValues,
      limits = cat_levels, breaks = cat_levels, drop = FALSE
    ) +
    scale_shape_manual(
      name = legName, values = shapeValues,
      limits = cat_levels, breaks = cat_levels, drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, xmax),
      expand = c(0.01, 0),
      breaks = x_breaks,
      labels = names(x_breaks),
      name = "Chromosome"
    ) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(
      axis.text  = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 12)
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 4)),
      shape  = guide_legend(override.aes = list(size = 4))
    ) +
    ggtitle(title)
  
  return(returnPlot)
}

run_one_analysis <- function(
    gwas_file,
    lfdr_file,
    analysis,
    title,
    ylimits,
    fig_name,
    aID,
    fdrLimitNew  = 0.1,
    p_cutoff_meta = 5e-8
) {
  
  three_Z <- fread(file.path(dataDir, gwas_file))
  tempNew <- fread(file.path(outputDir, lfdr_file), header = TRUE, data.table = FALSE)
  
  tempDat <- three_Z %>%
    mutate(origIdx = 1:nrow(.)) %>%
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr)) %>%
    mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
    arrange(origIdx)
  
  rejectDat <- tempDat %>% filter(rejNew == 1)
  
  para <- make_para(analysis)
  
  meta_res <- run_meta_pair(
    three_Z = three_Z,
    para    = para,
    aID     = aID,
    p_cutoff = p_cutoff_meta
  )
  
  allZ_meta_sig <- meta_res$sig %>%
    mutate(method = "Meta-analysis") %>%
    select(chrpos, method)
  
  merge_results <- tempDat %>%
    mutate(meta = ifelse(chrpos %in% allZ_meta_sig$chrpos, 1, 0)) %>%
    mutate(replication = ifelse(chrpos %in% rejectDat$chrpos, 1, 0)) %>%
    mutate(cat = case_when(
      meta == 1 & replication == 1 ~ "Both",
      meta == 1 & replication == 0 ~ "Meta-analysis",
      meta == 0 & replication == 1 ~ "Replication Analysis",
      TRUE                          ~ "Neither"
    )) %>%
    mutate(
      chars    = nchar(chrpos),
      colonPos = as.numeric(gregexpr(":", chrpos)),
      Chr      = substr(chrpos, 1, colonPos - 1),
      BP       = substr(chrpos, colonPos + 1, chars),
      chrom    = paste0("chr", Chr),
      pos      = as.numeric(BP),
      y        = -log10(cumNew)
    ) %>%
    select(chrom, pos, y, cat, meta, replication)
  
  merge_results$cat <- factor(
    merge_results$cat,
    levels = c("Both", "Meta-analysis", "Replication Analysis", "Neither")
  )
  
  n_meta             <- sum(merge_results$meta == 1, na.rm = TRUE)
  n_replication      <- sum(merge_results$replication == 1, na.rm = TRUE)
  n_both             <- sum(merge_results$meta == 1 & merge_results$replication == 1, na.rm = TRUE)
  n_meta_only        <- sum(merge_results$meta == 1 & merge_results$replication == 0, na.rm = TRUE)
  n_replication_only <- sum(merge_results$meta == 0 & merge_results$replication == 1, na.rm = TRUE)
  
  cat("\n==============================\n")
  cat("Analysis:", analysis, "\n")
  cat("Title:",    title,    "\n")
  cat("aID:",      aID,      "\n")
  cat("Pair:",     para$pair[aID], "\n")
  cat("Meta-analysis significant SNPs:", n_meta, "\n")
  cat("Replication discoveries:",        n_replication, "\n")
  cat("Both discoveries:",               n_both, "\n")
  cat("Meta-analysis only:",             n_meta_only, "\n")
  cat("Replication analysis only:",      n_replication_only, "\n")
  cat("==============================\n")
  
  colValues <- c(
    "Both"                 = wes_palette("Darjeeling1")[1],
    "Meta-analysis"        = wes_palette("Darjeeling1")[2],
    "Replication Analysis" = wes_palette("Darjeeling1")[5],
    "Neither"              = wes_palette("Darjeeling1")[3]
  )
  
  shapeValues <- c(
    "Both"                 = 16,
    "Meta-analysis"        = 17,
    "Replication Analysis" = 15,
    "Neither"              = 18
  )
  
  manPlot <- plotmanhattan(
    merge_results,
    build       = "hg19",
    colValues   = colValues,
    shapeValues = shapeValues,
    ylimits     = ylimits,
    legName     = "Method",
    title       = title
  )
  
  ggsave(
    filename = file.path(outputDir_Fig, fig_name),
    plot     = manPlot,
    width    = 6,
    height   = 3,
    device   = cairo_pdf
  )
  
  return(invisible(list(
    three_Z       = three_Z,
    tempDat       = tempDat,
    rejectDat     = rejectDat,
    meta          = meta_res,
    meta_sig      = allZ_meta_sig,
    merge_results = merge_results,
    plot          = manPlot,
    fig_path      = file.path(outputDir_Fig, fig_name),
    counts = list(
      n_meta             = n_meta,
      n_replication      = n_replication,
      n_both             = n_both,
      n_meta_only        = n_meta_only,
      n_replication_only = n_replication_only
    )
  )))
}

save_combined_plot <- function(res_overall, res_scc, res_adeno, outputDir_Fig, aID) {
  
  # Robust legend extractor (works on old and new cowplot)
  get_legend_robust <- function(plot) {
    leg <- tryCatch(
      cowplot::get_plot_component(plot, "guide-box-bottom", return_all = TRUE),
      error = function(e) NULL
    )
    if (is.null(leg) || inherits(leg, "zeroGrob")) {
      leg <- tryCatch(cowplot::get_legend(plot), error = function(e) NULL)
    }
    if (is.null(leg) || inherits(leg, "zeroGrob")) {
      leg <- tryCatch(ggpubr::get_legend(plot), error = function(e) NULL)
    }
    leg
  }
  
  # Strip each panel's own legend so they don't clutter the row
  p_overall <- res_overall$plot + theme(legend.position = "none")
  p_adeno   <- res_adeno$plot   + theme(legend.position = "none")
  p_scc     <- res_scc$plot     + theme(legend.position = "none")
  
  # Dummy plot supplies ONE shared legend with all four categories
  cat_levels <- c("Both", "Meta-analysis", "Replication Analysis", "Neither")
  
  colValues <- c(
    "Both"                 = wes_palette("Darjeeling1")[1],
    "Meta-analysis"        = wes_palette("Darjeeling1")[2],
    "Replication Analysis" = wes_palette("Darjeeling1")[5],
    "Neither"              = wes_palette("Darjeeling1")[3]
  )
  
  shapeValues <- c(
    "Both"                 = 16,
    "Meta-analysis"        = 17,
    "Replication Analysis" = 15,
    "Neither"              = 18
  )
  
  dummy_legend_data <- data.frame(
    x   = seq_along(cat_levels),
    y   = seq_along(cat_levels),
    cat = factor(cat_levels, levels = cat_levels)
  )
  
  legend_src <- ggplot(
    dummy_legend_data,
    aes(x = x, y = y, color = cat, shape = cat)
  ) +
    geom_point(size = 4) +
    scale_color_manual(
      name = "Method", values = colValues,
      limits = cat_levels, breaks = cat_levels, drop = FALSE
    ) +
    scale_shape_manual(
      name = "Method", values = shapeValues,
      limits = cat_levels, breaks = cat_levels, drop = FALSE
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 4)),
      shape  = guide_legend(override.aes = list(size = 4))
    ) +
    theme_cowplot() +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 11)
    )
  
  common_legend <- get_legend_robust(legend_src)
  
  # Top row: three panels with labels (order: Overall, LUAD, LUSC)
  top_row <- cowplot::plot_grid(
    p_overall,
    p_adeno,
    p_scc,
    ncol       = 3,
    labels     = c("A", "B", "C"),
    label_size = 14
  )
  
  final_plot <- cowplot::plot_grid(
    top_row,
    common_legend,
    ncol        = 1,
    rel_heights = c(1, 0.12)
  )
  
  pdf_file <- file.path(
    outputDir_Fig,
    paste0("two_rep_combined_aID", aID, ".pdf")
  )
  ggsave(
    filename = pdf_file,
    plot     = final_plot,
    width    = 12,
    height   = 4.5,
    device   = cairo_pdf
  )
  
  cat("\nSaved combined two-way figure:\n")
  cat(pdf_file, "\n")
  
  invisible(NULL)
}

# ============================================================
# Run analyses (uses aID from command line)
# ============================================================

fdrLimitNew   <- 0.1
p_cutoff_meta <- 5e-8

# -----------------------------
# Overall lung cancer
# -----------------------------
res_overall <- run_one_analysis(
  gwas_file = "lc_overall_three.txt",
  lfdr_file = paste0("two_rep_overall_aID", aID, "_newlfdr.txt"),
  analysis  = "overall",
  title     = "Overall",
  ylimits   = c(0, 70),
  fig_name  = paste0("two_rep_overall_aID", aID, ".png"),
  aID       = aID,
  fdrLimitNew   = fdrLimitNew,
  p_cutoff_meta = p_cutoff_meta
)

# -----------------------------
# LUSC / SCC
# -----------------------------
res_scc <- run_one_analysis(
  gwas_file = "lc_scc_three.txt",
  lfdr_file = paste0("two_rep_scc_aID", aID, "_newlfdr.txt"),
  analysis  = "LUSC",
  title     = "LUSC",
  ylimits   = c(0, 20),
  fig_name  = paste0("two_rep_scc_aID", aID, ".png"),
  aID       = aID,
  fdrLimitNew   = fdrLimitNew,
  p_cutoff_meta = p_cutoff_meta
)

# -----------------------------
# LUAD / Adenocarcinoma
# -----------------------------
res_adeno <- run_one_analysis(
  gwas_file = "lc_adeno_three.txt",
  lfdr_file = paste0("two_rep_adeno_aID", aID, "_newlfdr.txt"),
  analysis  = "LUAD",
  title     = "LUAD",
  ylimits   = c(0, 25),
  fig_name  = paste0("two_rep_adeno_aID", aID, ".png"),
  aID       = aID,
  fdrLimitNew   = fdrLimitNew,
  p_cutoff_meta = p_cutoff_meta
)

# -----------------------------
# Save combined plot (Order: Overall, LUAD, LUSC)
# -----------------------------
save_combined_plot(
  res_overall   = res_overall,
  res_scc       = res_scc,
  res_adeno     = res_adeno,
  outputDir_Fig = outputDir_Fig,
  aID           = aID
)