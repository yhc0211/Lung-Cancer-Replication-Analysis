# Make Fig6 - manhattan plot for FAVOR (SCC / LUSC)
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


# ---------------------------------------------------------------------------
# Zoom window for the chr15q25 locus (shared by replication + meta panels)
# Adjust if you want a wider/narrower view.
# ---------------------------------------------------------------------------
CHR15_BP_MIN <- 78.5e6
CHR15_BP_MAX <- 79.2e6


# ---------------------------------------------------------------------------
# Shared plotting function: chr15 Manhattan with Mb-labeled x-axis
# ---------------------------------------------------------------------------
# `rank_col` is the column name used for color/shape (AC_rank or AE_rank).
plotmanhattanFavor <- function(gwas, build = c('hg18','hg19','hg38'),
                               rank_col, colValues, shapeValues,
                               ylimits, legName, title) {
  data <- gwas
  build <- match.arg(build)
  data <- add_cumulative_pos(data, build)
  
  # window bounds (small padding on each side)
  xmin <- min(data$cumulative_pos) - 20000
  xmax <- max(data$cumulative_pos) + 20000
  brks <- pretty(c(xmin, xmax), n = 4)
  
  # convert cumulative_pos breaks back to chr15 BP (Mb) for axis labels
  chr15_offset <- min(data$cumulative_pos) - min(data$BP)
  lbls <- paste0(round((brks - chr15_offset) / 1e6, 1), " Mb")
  
  ggplot(data, aes(x = cumulative_pos, y = y,
                   color = .data[[rank_col]],
                   shape = .data[[rank_col]])) +
    geom_point(size = 2) +
    labs(
      title = title,
      x = "Chromosome 15 (Mb)",
      y = expression(paste("-", log[10], "(FDR)"))
    ) +
    scale_color_manual(name = legName, values = colValues) +
    scale_shape_manual(name = legName, values = shapeValues) +
    scale_x_continuous(
      limits = c(xmin, xmax),
      expand = c(0.01, 0),
      breaks = brks,
      labels = lbls,
      name   = "Chromosome 15 (Mb)"
    ) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(
      axis.text    = element_text(size = 12),
      axis.title   = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12),
      axis.line.x  = element_line(color = "black", linewidth = 0.5),
      axis.line.y  = element_line(color = "black", linewidth = 0.5),
      axis.ticks   = element_line(color = "black")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4)))
}


# ===========================================================================
# Replication analysis (SCC lung)
# ===========================================================================
three_Z <- fread(here::here(dataDir, "lc_scc_three.txt"))
tempNew <- fread(paste0(outputDir, "/rep_three_scc_aID1_newlfdr.txt"),
                 header = TRUE, data.table = FALSE)

tempDat <- three_Z %>%
  mutate(origIdx = 1:nrow(.)) %>%
  mutate(newLfdr = tempNew$x) %>%
  arrange(newLfdr) %>%
  mutate(cumNew = cummean(newLfdr)) %>%
  mutate(rejNew = ifelse(cumNew < 0.1, 1, 0)) %>%
  arrange(origIdx)
rejectDat <- tempDat %>% filter(rejNew == 1)

results_hg38 <- fread(here::here(outputDir, "hg38_scc.bed"), header = FALSE)

hg_ref <- cbind(rejectDat, results_hg38) %>%
  select(chrpos, V1) %>%
  separate(V1, c("chrom", "BP_hg38"), sep = "-")

manDataFavor <- readxl::read_xlsx(here::here(outputDir, "FAVOR_results_scc.xlsx")) %>%
  rename(Chr = Chromosome, BP = Position) %>%
  cbind(hg_ref) %>%
  select("chrpos", "ApcConservationV2", "ApcEpigeneticsActive")

manDataFavor <- rejectDat %>%
  left_join(manDataFavor, by = "chrpos")

manDataFavor <- manDataFavor %>%
  mutate(AC_rank = case_when(
    is.na(ApcConservationV2)                              ~ "Not Significant",
    ApcConservationV2 > 25                                ~ "high (>25)",
    ApcConservationV2 > 10 & ApcConservationV2 <= 25      ~ "medium (10 - 25)",
    TRUE                                                  ~ "low (< 10)"
  )) %>%
  mutate(AE_rank = case_when(
    is.na(ApcEpigeneticsActive)                                  ~ "Not Significant",
    ApcEpigeneticsActive > 25                                    ~ "high (>25)",
    ApcEpigeneticsActive > 10 & ApcEpigeneticsActive <= 25       ~ "medium (10 - 25)",
    TRUE                                                         ~ "low (< 10)"
  ))

manDataFavor$AC_rank <- forcats::fct_relevel(manDataFavor$AC_rank,
                                             c("high (>25)", "medium (10 - 25)", "low (< 10)", "Not Significant"))
manDataFavor$AE_rank <- forcats::fct_relevel(manDataFavor$AE_rank,
                                             c("high (>25)", "medium (10 - 25)", "low (< 10)", "Not Significant"))

manDataFavor <- manDataFavor[order(manDataFavor$AE_rank, decreasing = TRUE), ]
manDataFavor <- manDataFavor[order(manDataFavor$AC_rank, decreasing = TRUE), ]

# Chromosome 15 replication data
plotDatFavor <- manDataFavor %>%
  mutate(BP = as.numeric(pos)) %>%
  mutate(y = -log10(cumNew)) %>%
  mutate(chrom = paste0("chr", chrom)) %>%
  filter(chrom == "chr15")

manPlotFavor_AE_SCC <- plotmanhattanFavor(
  plotDatFavor, build = "hg19", rank_col = "AE_rank",
  colValues  = wes_palette("FantasticFox1")[c(5, 3, 4)],
  shapeValues = c(16, 17, 18),
  ylimits = c(0, 10),
  legName = "Epigenetics Score",
  title   = "Replication analysis"
)

manPlotFavor_AC_SCC <- plotmanhattanFavor(
  plotDatFavor, build = "hg19", rank_col = "AC_rank",
  colValues  = wes_palette("FantasticFox1")[c(5, 3, 4)],
  shapeValues = c(16, 17, 18),
  ylimits = c(0, 10),
  legName = "Conservation Score",
  title   = "Replication analysis"
)


# ===========================================================================
# Meta-analysis (SCC / LUSC)
# ===========================================================================
allZ_meta <- three_Z %>%
  rename(P1 = pILCCO, Z1 = zILCCO,
         P2 = pUKB,   Z2 = zUKB,
         P3 = pMVP,   Z3 = zMVP) %>%
  select(P1, Z1, P2, Z2, P3, Z3, chrpos)

calc_neff <- function(n_case, n_ctrl) {
  4 / (1 / n_case + 1 / n_ctrl)
}

N_ILCCO <- calc_neff(7426, 56450)
N_UKB   <- calc_neff(489,  330018)
N_MVP   <- calc_neff(1475, 62708)
w1 <- sqrt(N_ILCCO)
w2 <- sqrt(N_UKB)
w3 <- sqrt(N_MVP)

allZ_meta_sig <- allZ_meta %>%
  mutate(
    Z_meta = (Z1 * w1 + Z2 * w2 + Z3 * w3) / sqrt(w1^2 + w2^2 + w3^2),
    P_meta = 2 * pnorm(abs(Z_meta), lower.tail = FALSE)
  ) %>%
  filter(P_meta < 5e-8)

results_hg38 <- fread(here::here(outputDir, "hg38_lusc_meta_revised.bed"),
                      header = FALSE)

hg_ref <- cbind(allZ_meta_sig, results_hg38) %>%
  select(chrpos, V1) %>%
  separate(V1, c("chrom", "BP_hg38"), sep = "-")

manDataFavor <- readxl::read_xlsx(
  here::here(outputDir, "FAVOR_lusc_meta_revised.xlsx")
) %>%
  rename(Chr = Chromosome, BP = Position) %>%
  mutate(Chr = as.character(Chr), BP = as.numeric(BP))

hg_ref2 <- hg_ref %>%
  mutate(Chr = sub(":.*$", "", chrpos),
         BP  = as.numeric(BP_hg38))

manDataFavor <- manDataFavor %>%
  inner_join(hg_ref2, by = c("Chr", "BP")) %>%
  select(chrpos, ApcConservationV2, ApcEpigeneticsActive)

manDataFavor <- allZ_meta_sig %>%
  left_join(manDataFavor, by = "chrpos")

manDataFavor <- manDataFavor %>%
  mutate(AC_rank = case_when(
    is.na(ApcConservationV2)                              ~ "Not Significant",
    ApcConservationV2 > 25                                ~ "high (>25)",
    ApcConservationV2 > 10 & ApcConservationV2 <= 25      ~ "medium (10 - 25)",
    TRUE                                                  ~ "low (< 10)"
  )) %>%
  mutate(AE_rank = case_when(
    is.na(ApcEpigeneticsActive)                                  ~ "Not Significant",
    ApcEpigeneticsActive > 25                                    ~ "high (>25)",
    ApcEpigeneticsActive > 10 & ApcEpigeneticsActive <= 25       ~ "medium (10 - 25)",
    TRUE                                                         ~ "low (< 10)"
  ))

manDataFavor$AC_rank <- forcats::fct_relevel(manDataFavor$AC_rank,
                                             c("high (>25)", "medium (10 - 25)", "low (< 10)", "Not Significant"))
manDataFavor$AE_rank <- forcats::fct_relevel(manDataFavor$AE_rank,
                                             c("high (>25)", "medium (10 - 25)", "low (< 10)", "Not Significant"))

manDataFavor <- manDataFavor[order(manDataFavor$AE_rank, decreasing = TRUE), ]
manDataFavor <- manDataFavor[order(manDataFavor$AC_rank, decreasing = TRUE), ]

# Chromosome 15 meta data
# - add BP = pos so the shared plot function can compute the Mb offset
# - zoom into the chr15q25 locus (same window as replication panels)
plotDatFavor <- manDataFavor %>%
  separate(chrpos, c("Chr", "pos"), sep = ":") %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(BP  = pos) %>%
  mutate(y   = -log10(p.adjust(P_meta, method = "BH"))) %>%
  mutate(chrom = paste0("chr", Chr)) %>%
  filter(chrom == "chr15") %>%
  filter(BP >= CHR15_BP_MIN, BP <= CHR15_BP_MAX)        # <-- zoom

manPlotFavor_AE_SCC_meta <- plotmanhattanFavor(
  plotDatFavor, build = "hg19", rank_col = "AE_rank",
  colValues  = wes_palette("FantasticFox1")[c(5, 3, 4)],
  shapeValues = c(16, 17, 18),
  ylimits = c(0, 60),
  legName = "Epigenetics Score",
  title   = "Meta-analysis"
)

manPlotFavor_AC_SCC_meta <- plotmanhattanFavor(
  plotDatFavor, build = "hg19", rank_col = "AC_rank",
  colValues  = wes_palette("FantasticFox1")[c(5, 3, 4)],
  shapeValues = c(16, 17, 18),
  ylimits = c(0, 60),
  legName = "Conservation Score",
  title   = "Meta-analysis"
)


# ===========================================================================
# Assemble and save
# ===========================================================================
merge_plot1 <- ggarrange(manPlotFavor_AC_SCC, manPlotFavor_AC_SCC_meta,
                         ncol = 2, nrow = 1,
                         common.legend = TRUE, legend = "bottom",
                         labels = c("A", "B"))
merge_plot2 <- ggarrange(manPlotFavor_AE_SCC, manPlotFavor_AE_SCC_meta,
                         ncol = 2, nrow = 1,
                         common.legend = TRUE, legend = "bottom",
                         labels = c("C", "D"))
merge_plot  <- ggarrange(merge_plot1, merge_plot2,
                         ncol = 1, nrow = 2,
                         common.legend = TRUE, legend = "bottom")

ggsave(paste0(outputDir_Fig, "/FAVOR_LUSC_revised.pdf"),
       plot = merge_plot, width = 15, height = 10)
