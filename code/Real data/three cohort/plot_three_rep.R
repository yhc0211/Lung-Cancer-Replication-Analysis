# Make Three-way Replication Analysis Results Manhattan plot

# ============================================================
# Setup
# ============================================================

setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/plot_three_rep.R")

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

# Source the .R scripts from the SupportingCode/ folder
codePath <- here::here("SupportingCode")
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# Directories
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

make_neff_three <- function(analysis = c("overall", "LUAD", "LUSC")) {
  
  analysis <- match.arg(analysis)
  
  if (analysis == "overall") {
    N_ILCCO <- calc_neff(29266, 56450)
    N_MVP   <- calc_neff(10398, 62708)
    N_UKB   <- calc_neff(2404, 330018)
  }
  
  if (analysis == "LUAD") {
    N_ILCCO <- calc_neff(11273, 56450)
    N_MVP   <- calc_neff(2019, 62708)
    N_UKB   <- calc_neff(862, 330018)
  }
  
  if (analysis == "LUSC") {
    N_ILCCO <- calc_neff(7426, 56450)
    N_MVP   <- calc_neff(1475, 62708)
    N_UKB   <- calc_neff(489, 330018)
  }
  
  list(
    ILCCO = N_ILCCO,
    MVP = N_MVP,
    UKB = N_UKB
  )
}

run_meta_three <- function(three_Z,
                           analysis = c("overall", "LUAD", "LUSC"),
                           p_cutoff = 5e-8) {
  
  analysis <- match.arg(analysis)
  neff <- make_neff_three(analysis)
  
  w_ILCCO <- sqrt(neff$ILCCO)
  w_UKB   <- sqrt(neff$UKB)
  w_MVP   <- sqrt(neff$MVP)
  
  out <- three_Z %>%
    transmute(
      chrpos = chrpos,
      RS = if ("RS" %in% colnames(three_Z)) RS else NA,
      zILCCO = zILCCO,
      zUKB = zUKB,
      zMVP = zMVP,
      pILCCO = pILCCO,
      pUKB = pUKB,
      pMVP = pMVP
    ) %>%
    filter(
      !is.na(zILCCO),
      !is.na(zUKB),
      !is.na(zMVP)
    ) %>%
    mutate(
      Z_meta = (
        zILCCO * w_ILCCO +
          zUKB   * w_UKB +
          zMVP   * w_MVP
      ) / sqrt(w_ILCCO^2 + w_UKB^2 + w_MVP^2),
      P_meta = 2 * pnorm(abs(Z_meta), lower.tail = FALSE),
      method = "Meta-analysis"
    )
  
  out_sig <- out %>%
    filter(P_meta < p_cutoff)
  
  list(
    all = out,
    sig = out_sig
  )
}

make_replication_result <- function(three_Z,
                                    lfdr_file,
                                    fdrLimitNew = 0.1) {
  
  tempNew <- fread(
    file.path(outputDir, lfdr_file),
    header = TRUE,
    data.table = FALSE
  )
  
  tempDat <- three_Z %>%
    dplyr::select(chrpos, zILCCO, zUKB, zMVP) %>%
    mutate(origIdx = 1:nrow(.)) %>%
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr)) %>%
    mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
    arrange(origIdx)
  
  rejectDat <- tempDat %>%
    filter(rejNew == 1)
  
  list(
    tempDat = tempDat,
    rejectDat = rejectDat
  )
}

make_manhattan_data <- function(three_Z,
                                study = c("ILCCO", "MVP", "UKB"),
                                meta_sig_chrpos,
                                replication_chrpos,
                                p_filter = 1e-4) {
  
  study <- match.arg(study)
  
  if (study == "ILCCO") {
    p_col <- "pILCCO"
  }
  
  if (study == "MVP") {
    p_col <- "pMVP"
  }
  
  if (study == "UKB") {
    p_col <- "pUKB"
  }
  
  out <- three_Z %>%
    transmute(
      chrpos = chrpos,
      chrom = chrom,
      pos = pos,
      p = as.numeric(.data[[p_col]])
    ) %>%
    filter(
      !is.na(p),
      p > 0,
      p < p_filter
    ) %>%
    mutate(
      meta = ifelse(chrpos %in% meta_sig_chrpos, 1, 0),
      replication = ifelse(chrpos %in% replication_chrpos, 1, 0),
      cat = case_when(
        meta == 1 & replication == 1 ~ "Both",
        meta == 1 & replication == 0 ~ "Meta-analysis",
        meta == 0 & replication == 1 ~ "Replication Analysis",
        TRUE ~ "Neither"
      ),
      y = -log10(p),
      chrom = paste0("chr", chrom)
    ) %>%
    filter(!is.na(y), is.finite(y))
  
  out$cat <- factor(
    out$cat,
    levels = c("Both", "Meta-analysis", "Replication Analysis", "Neither")
  )
  
  out
}

plotmanhattan_three <- function(gwas,
                                build = c("hg18", "hg19", "hg38"),
                                ylimits,
                                title,
                                show_legend = FALSE) {
  
  data <- gwas
  build <- match.arg(build)
  
  cat_levels <- c("Both", "Meta-analysis", "Replication Analysis", "Neither")
  data$cat <- factor(data$cat, levels = cat_levels)
  
  colValues <- c(
    "Both" = wes_palette("Darjeeling1")[1],
    "Meta-analysis" = wes_palette("Darjeeling1")[2],
    "Replication Analysis" = wes_palette("Darjeeling1")[5],
    "Neither" = wes_palette("Darjeeling1")[3]
  )
  
  shapeValues <- c(
    "Both" = 16,
    "Meta-analysis" = 17,
    "Replication Analysis" = 15,
    "Neither" = 18
  )
  
  data <- add_cumulative_pos(data, build)
  chrom_lengths <- get_chrom_lengths(build)
  xmax <- get_total_length(chrom_lengths)
  x_breaks <- get_x_breaks(chrom_lengths)
  
  names(x_breaks)[c(9, 11, 13, 14, 16, 17, 18, 20, 21)] <- ""
  
  p <- ggplot(
    data,
    aes(
      x = cumulative_pos,
      y = y,
      color = cat,
      shape = cat
    )
  ) +
    geom_point(size = 0.75, alpha = 0.85) +
    xlab("Chromosome") +
    ylab(expression(paste("-", log[10], "(P)"))) +
    scale_color_manual(
      name = "Method",
      values = colValues,
      limits = cat_levels,
      breaks = cat_levels,
      drop = FALSE
    ) +
    scale_shape_manual(
      name = "Method",
      values = shapeValues,
      limits = cat_levels,
      breaks = cat_levels,
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, xmax),
      expand = c(0.01, 0),
      breaks = x_breaks,
      labels = names(x_breaks),
      name = "Chromosome"
    ) +
    coord_cartesian(ylim = ylimits) +
    theme_cowplot() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 11)
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 4)),
      shape = guide_legend(override.aes = list(size = 4))
    ) +
    ggtitle(title)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "bottom")
  }
  
  p
}

safe_ylim <- function(dat, lower = 3.9, upper_pad = 1) {
  upper <- max(dat$y, na.rm = TRUE) + upper_pad
  if (!is.finite(upper) || upper <= lower) {
    upper <- lower + upper_pad
  }
  c(lower, upper)
}

run_one_threeway_analysis <- function(gwas_file,
                                      lfdr_file,
                                      analysis = c("overall", "LUAD", "LUSC"),
                                      p_cutoff_meta = 5e-8,
                                      fdrLimitNew = 0.1,
                                      show_legend_row = FALSE) {
  
  analysis <- match.arg(analysis)
  
  three_Z <- fread(file.path(dataDir, gwas_file))
  
  rep_res <- make_replication_result(
    three_Z = three_Z,
    lfdr_file = lfdr_file,
    fdrLimitNew = fdrLimitNew
  )
  
  meta_res <- run_meta_three(
    three_Z = three_Z,
    analysis = analysis,
    p_cutoff = p_cutoff_meta
  )
  
  meta_sig_chrpos <- unique(meta_res$sig$chrpos)
  replication_chrpos <- unique(rep_res$rejectDat$chrpos)
  
  n_meta <- length(meta_sig_chrpos)
  n_replication <- length(replication_chrpos)
  n_both <- length(intersect(meta_sig_chrpos, replication_chrpos))
  n_meta_only <- length(setdiff(meta_sig_chrpos, replication_chrpos))
  n_replication_only <- length(setdiff(replication_chrpos, meta_sig_chrpos))
  
  cat("\n==============================\n")
  cat("Three-way analysis:", analysis, "\n")
  cat("Meta-analysis significant SNPs:", n_meta, "\n")
  cat("Replication discoveries:", n_replication, "\n")
  cat("Both discoveries:", n_both, "\n")
  cat("Meta-analysis only:", n_meta_only, "\n")
  cat("Replication analysis only:", n_replication_only, "\n")
  cat("==============================\n")
  
  manDatILCCO <- make_manhattan_data(
    three_Z = three_Z,
    study = "ILCCO",
    meta_sig_chrpos = meta_sig_chrpos,
    replication_chrpos = replication_chrpos,
    p_filter = 1e-4
  )
  
  manDatMVP <- make_manhattan_data(
    three_Z = three_Z,
    study = "MVP",
    meta_sig_chrpos = meta_sig_chrpos,
    replication_chrpos = replication_chrpos,
    p_filter = 1e-4
  )
  
  manDatUKB <- make_manhattan_data(
    three_Z = three_Z,
    study = "UKB",
    meta_sig_chrpos = meta_sig_chrpos,
    replication_chrpos = replication_chrpos,
    p_filter = 1e-4
  )
  
  plotILCCO <- plotmanhattan_three(
    manDatILCCO,
    build = "hg19",
    ylimits = safe_ylim(manDatILCCO),
    title = "ILCCO",
    show_legend = FALSE
  )
  
  plotMVP <- plotmanhattan_three(
    manDatMVP,
    build = "hg19",
    ylimits = safe_ylim(manDatMVP),
    title = "MVP",
    show_legend = FALSE
  )
  
  plotUKB <- plotmanhattan_three(
    manDatUKB,
    build = "hg19",
    ylimits = safe_ylim(manDatUKB),
    title = "UKB",
    show_legend = show_legend_row
  )
  
  row_plot <- ggarrange(
    plotILCCO,
    plotMVP,
    plotUKB,
    ncol = 3,
    nrow = 1,
    legend = "none"
  )
  
  invisible(list(
    three_Z = three_Z,
    replication = rep_res,
    meta = meta_res,
    manDatILCCO = manDatILCCO,
    manDatMVP = manDatMVP,
    manDatUKB = manDatUKB,
    plotILCCO = plotILCCO,
    plotMVP = plotMVP,
    plotUKB = plotUKB,
    row_plot = row_plot,
    counts = list(
      n_meta = n_meta,
      n_replication = n_replication,
      n_both = n_both,
      n_meta_only = n_meta_only,
      n_replication_only = n_replication_only
    )
  ))
}

save_threeway_combined_plot <- function(res_overall,
                                        res_luad,
                                        res_lusc,
                                        output_file = "three_rep_combined.pdf") {
  
  # ------------------------------------------------------------
  # Robust legend extractor (works on old and new cowplot)
  # ------------------------------------------------------------
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
  
  # ------------------------------------------------------------
  # Helper: add left-aligned row label
  # ------------------------------------------------------------
  add_left_label <- function(plot, label) {
    ggdraw() +
      draw_label(
        label,
        x = 0.01, y = 0.98,
        hjust = 0, vjust = 1,
        fontface = "bold", size = 15
      ) +
      draw_plot(plot, x = 0, y = 0, width = 1, height = 0.92)
  }
  
  row_overall <- add_left_label(res_overall$row_plot, "A. Overall")
  row_luad    <- add_left_label(res_luad$row_plot,    "B. LUAD")
  row_lusc    <- add_left_label(res_lusc$row_plot,    "C. LUSC")
  
  # ------------------------------------------------------------
  # Dummy plot whose only job is to supply ONE legend with all
  # four categories (Both, Meta-analysis, Replication Analysis,
  # Neither) regardless of what appears in the real panels.
  # ------------------------------------------------------------
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
      legend.title    = element_text(size = 11),
      legend.text     = element_text(size = 10)
    )
  
  common_legend <- get_legend_robust(legend_src)
  
  # ------------------------------------------------------------
  # Stack the three rows (each row has legend = "none" already)
  # ------------------------------------------------------------
  main_plot <- ggarrange(
    row_overall,
    row_luad,
    row_lusc,
    ncol    = 1,
    nrow    = 3,
    heights = c(1, 1, 1),
    legend  = "none"
  )
  
  # ------------------------------------------------------------
  # Final figure = main panels on top, ONE legend on the bottom.
  # plot_grid handles grobs more reliably than ggarrange here.
  # ------------------------------------------------------------
  final_plot <- cowplot::plot_grid(
    main_plot,
    common_legend,
    ncol        = 1,
    rel_heights = c(1, 0.08)
  )
  
  pdf_file <- file.path(outputDir_Fig, output_file)
  
  ggsave(
    filename = pdf_file,
    plot     = final_plot,
    width    = 8,
    height   = 11,
    device   = "pdf"
  )
  
  cat("\nSaved combined three-way figure:\n")
  cat(pdf_file, "\n")
  
  invisible(NULL)
}# ============================================================
# Run analyses
# ============================================================

fdrLimitNew <- 0.1
p_cutoff_meta <- 5e-8

res_overall <- run_one_threeway_analysis(
  gwas_file = "lc_overall_three.txt",
  lfdr_file = "three_rep_data_aID1_newlfdr.txt",
  analysis = "overall",
  p_cutoff_meta = p_cutoff_meta,
  fdrLimitNew = fdrLimitNew,
  show_legend_row = FALSE
)

res_luad <- run_one_threeway_analysis(
  gwas_file = "lc_adeno_three.txt",
  lfdr_file = "rep_three_adeno_aID1_newlfdr.txt",
  analysis = "LUAD",
  p_cutoff_meta = p_cutoff_meta,
  fdrLimitNew = fdrLimitNew,
  show_legend_row = FALSE
)

res_lusc <- run_one_threeway_analysis(
  gwas_file = "lc_scc_three.txt",
  lfdr_file = "rep_three_scc_aID1_newlfdr.txt",
  analysis = "LUSC",
  p_cutoff_meta = p_cutoff_meta,
  fdrLimitNew = fdrLimitNew,
  show_legend_row = FALSE
)

# Save final combined PDF
save_threeway_combined_plot(
  res_overall = res_overall,
  res_luad = res_luad,
  res_lusc = res_lusc,
  output_file = "three_rep_combined.pdf"
)




#FAVOR Part
three_Z <- fread(here::here(dataDir, "lc_overall_three.txt"))

meta_favor_overall <- res_overall$meta$sig

allZ_meta_sig <- meta_favor_overall %>%
  separate(chrpos,c("Chr","BP"),sep = ":") %>%
  dplyr::select(Chr,BP) %>%
  mutate(final = paste0("chr",Chr,":",BP,"-",BP))

#write.table(allZ_meta_sig$final, paste0(outputDir, "/for_liftover_overall_meta_revised.txt"), append=F, quote=F, row.names=F, col.names=F)

###for FAVOR 
# results_hg38 <- fread(here::here(outputDir, "hg38_overall_meta_revised.bed"),header=F)
# 
# FAVOR_ref <- cbind(allZ_meta_sig,results_hg38)
# FAVOR_ref <- FAVOR_ref %>%
#   rename(hg19 = final,hg38 = V1)
# 
# for_FAVOR <- FAVOR_ref %>%
#   mutate(chrpos = paste0(Chr,":",BP)) %>%
#   left_join(three_Z[,c("chrpos","A1_ILCCO","A2_ILCCO")],by = c("chrpos")) %>%
#   separate(hg38, c('tempChr', 'FAVOR'), sep='-') %>%
#   mutate(FAVOR1 = paste0(Chr,"-",FAVOR,"-",A1_ILCCO,"-",A2_ILCCO)) %>%
#   mutate(FAVOR2 = paste0(Chr,"-",FAVOR,"-",A2_ILCCO,"-",A1_ILCCO))   
# 
# final <- append(for_FAVOR$FAVOR1,for_FAVOR$FAVOR2)#1-based
# #write.table(final, paste0(outputDir, "/for_FAVOR_meta_three_way_overall_revised.txt"), append=F, quote=F, row.names=F, col.names=F)
# 
# 
# 
# 
# 
# 
# three_Z <- fread(here::here(dataDir, "lc_scc_three.txt"))
# meta_favor_scc <- res_lusc$meta$sig
# 
# allZ_meta_sig <- meta_favor_scc %>%
#   separate(chrpos,c("Chr","BP"),sep = ":") %>%
#   dplyr::select(Chr,BP) %>%
#   mutate(final = paste0("chr",Chr,":",BP,"-",BP))
# 
# #write.table(allZ_meta_sig$final, paste0(outputDir, "/for_liftover_lusc_meta_revised.txt"), append=F, quote=F, row.names=F, col.names=F)
# 
# ###for FAVOR 
# results_hg38 <- fread(here::here(outputDir, "hg38_lusc_meta_revised.bed"),header=F)
# 
# FAVOR_ref <- cbind(allZ_meta_sig,results_hg38)
# FAVOR_ref <- FAVOR_ref %>%
#   rename(hg19 = final,hg38 = V1)
# 
# for_FAVOR <- FAVOR_ref %>%
#   mutate(chrpos = paste0(Chr,":",BP)) %>%
#   left_join(three_Z[,c("chrpos","A1_ILCCO","A2_ILCCO")],by = c("chrpos")) %>%
#   separate(hg38, c('tempChr', 'FAVOR'), sep='-') %>%
#   mutate(FAVOR1 = paste0(Chr,"-",FAVOR,"-",A1_ILCCO,"-",A2_ILCCO)) %>%
#   mutate(FAVOR2 = paste0(Chr,"-",FAVOR,"-",A2_ILCCO,"-",A1_ILCCO))   
# 
# final <- append(for_FAVOR$FAVOR1,for_FAVOR$FAVOR2)#1-based
# #write.table(final, paste0(outputDir, "/for_FAVOR_meta_three_way_lusc_revised.txt"), append=F, quote=F, row.names=F, col.names=F)
# 
