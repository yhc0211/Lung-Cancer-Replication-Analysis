# calculate_eqtl_stats_synthetic.R

# For each of the 18 sentinel SNPs identified by ILCCO, test if it is
# associated with the expressions of a group of genes at the other 17 loci in lung.

# We cannot share the original GTEx data due to privacy concerns. Thus we have prepared
# a synthetic substitute, synthDat.txt. Using the synthetic data will not reproduce the 
# results in the paper. However, if the user possesses the original GTEx data, using it
# will reproduce our findings.

# We have placed the original output (summary statistics) of this step in the Data/ folder. Using that output
# will reproduce our findings.

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/calculate_eqtl_stats_synthetic.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("csmGmm_reproduce/Lung")
here::i_am("Lung/calculate_eqtl_stats.R")

library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)
library(data.table)
library(here)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output and data directories
outputDir <- here::here("Fig4", "output")
dataDir <- here::here("Data")

# the TWAS part of the mediation data (see DATA folder)
twasLoc <- here::here(dataDir, "scc_lung_addchr1.csv") 

# the gene location information (see DATA folder)
geneInfoLoc <- here::here(dataDir, "ensembl_refgene_hg19_20180109.rda")  

# where is the synthetic data expression data (see DATA folder)
synthLoc <- here::here(dataDir, "synthDat.txt")  

# original=1 means we have the real data
original <- 0

# some parameters for analysis
sentinelWindow <- 5000
cisLength <- 2500000

# read gds format data
read_gds <- function(gdsName) {
  gds <- GdsGenotypeReader(gdsName, YchromCode=24L, XYchromCode=25L)
  
  # all information packed into the file by GTEx 
  scanID <- getScanID(gds)
  family_info <- getVariable(gds, "sample.annot/family")
  father_info <- getVariable(gds, "sample.annot/father")
  mother_info <- getVariable(gds, "sample.annot/mother")
  # sex must be recorded as M/F/NA - it appears to be all "" in this dataset
  sex <- getVariable(gds, "sample.annot/sex")
  sex[sex == ""] <- NA
  phenotype <- getVariable(gds, "sample.annot/phenotype")
  
  # have to put the annotation data in a special structure, that special structure will go inside
  # another special structure, all this to satisfy the file container requirements
  scanAnnot <- data.frame(scanID, father_info, mother_info, sex, phenotype, stringsAsFactors = FALSE) %>%
    ScanAnnotationDataFrame()
  
  # some information about the SNPs, positions + alleles
  snpID <- getSnpID(gds)
  chromosome <- getChromosome(gds)
  position <- getPosition(gds)
  alleleA <- getAlleleA(gds)
  alleleB <- getAlleleB(gds)
  rsID <- getVariable(gds, "snp.rs.id")
  snpAnnot <- data.frame(snpID, chromosome, position, rsID, alleleA, alleleB,
                         stringsAsFactors = FALSE) %>%
    SnpAnnotationDataFrame(YchromCode=24L, XYchromCode=25L)
  
  # final object to be returned
  genoData <- GenotypeData(gds, scanAnnot = scanAnnot, snpAnnot = snpAnnot)
  # return
  return(list(genoData = genoData, gds = gds))
}

# clean and subset the genotype and gene expression data
clean_and_subset <- function(exprGeneNames, tissue, genoData, snpsToGet,
                             covarDF, exprMeta, exprDF) {
  
  # need the snpIDs of the SNPs that we specified
  snpDF <- c()
  uniqueChrs <- unique(snpsToGet$Chr)
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempToGet <- snpsToGet %>% filter(Chr == tempChr)
    
    # specify SNPs precisely by position 
    tempPos <- tempToGet %>% select(BP) %>% unlist(.)
    tempDF <- snpPosDF %>% filter(chromosome == tempChr & position %in% tempPos)
    snpDF <- rbind(snpDF, tempDF)
  }
  
  # pull the genotypes from the gds 
  gMat <- getGenotypeSelection(genoData, snpID = snpDF$snpID)
  
  # check names
  if (class(gMat)[1] == "matrix") { 
    stopifnot(rownames(gMat) == as.character(snpDF$snpID))
    if (nrow(gMat) == 0) {return(-1)} 
  } else {
    tempNames <- names(gMat)
    gMat <- matrix(data=gMat, nrow=1)
    colnames(gMat) <- tempNames
  }
  # use rsIDs
  rownames(gMat) <- snpDF$rsID
  
  # extract expression relevant tissue samples
  tissueExprSamps <- exprMeta %>% filter(Note == tissue)
  
  # overlap genotype and expression subjects
  overlapIDs <- intersect(tissueExprSamps$subj_id, colnames(gMat))
  
  # get the overlap data, pick a single tissue sample if there are multiple
  tissueExprOverlap <- tissueExprSamps %>% filter(subj_id %in% overlapIDs) %>%
    filter(Sample %in% colnames(exprDF)) %>%
    unique(by = "subj_id") 
  
  # after this, exprSelected is p*(n+2) where rows are genes and columns are name/description of gene + sample ID
  # gMatSelected is n*q where rows are subjects and columns are snpID
  exprSelected <- exprDF %>% as.data.frame(.) %>%
    select(Name, Description, all_of(tissueExprOverlap$Sample)) %>%
    filter(Description %in% exprGeneNames) %>%
    # sometimes two columns for the same gene
    distinct_at("Description", .keep_all = TRUE) %>% 
    column_to_rownames(var = "Description") %>% 
    select(-Name) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "Sample") %>%
    merge(., tissueExprOverlap %>% select(Sample, subj_id), by="Sample") %>%
    select(-Sample)
  
  # GTEx did inverse normal transform for each gene across samples
  n <- nrow(exprSelected)
  # there are actually no NAs in the entire gene tpm file, so we don't have to worry about that
  saveID <- exprSelected$subj_id
  exprINT <- exprSelected %>% select(-subj_id)
  for (gene_it in 1:ncol(exprINT)) {
    tempCol <- exprINT[, gene_it]
    # use the blom offset of k=3/8
    tempTrans <- qnorm( (rank(tempCol) - 3/8) / (n + 0.25) )
    exprINT[, gene_it] <- tempTrans
  }
  # add back subject ID
  exprINT <- exprINT %>% mutate(subj_id = saveID)
  
  # select and transpose genotype matrix
  genotypeSubToSelect <- colnames(gMat)[which(colnames(gMat) %in% tissueExprOverlap$subj_id)]
  gMatSelected <- as.data.frame(gMat) %>% select(all_of(genotypeSubToSelect)) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "subj_id") %>%
    # can't have numbers be columns for lm()
    set_colnames(paste0("rs_", colnames(.)))
  IDcol <- which(colnames(gMatSelected) == "rs_subj_id")
  colnames(gMatSelected)[IDcol] <- "subj_id"
  
  # covariates
  covarSelected <- covarDF %>% column_to_rownames(var="ID") %>%
    t(.) %>%
    as.data.frame(.) %>% 
    rownames_to_column(var="subj_id")
  # there's a gene named C2
  cIdx1 <- which(colnames(covarSelected) == "C2")
  colnames(covarSelected)[cIdx1] <- "Cov2"
  # there's also a gene named C3
  cIdx2 <- which(colnames(covarSelected) == "C3")
  colnames(covarSelected)[cIdx2] <- "Cov3"
  
  # merge
  allDat <- merge(gMatSelected, exprINT, by="subj_id") %>%
    merge(covarSelected, by="subj_id")
  
  return(list(allDat = allDat, snpDF = snpDF, exprINT = exprINT, gMatSelected = gMatSelected))
}

# do the linear regression
MPreg <- function(covarNames, genoNames, exprNames, allDat) {
  
  nGeno <- length(genoNames)
  
  # string for regression formula
  formulaRoot <- paste0(" ~ ", paste(covarNames, collapse=" + "))
  
  # two loops, outer is expression, inner is SNP
  resultsDF <- c() 
  setbasedDF <- c() 
  for (snp_it in 1:nGeno) {
    
    tempSNP <- genoNames[snp_it]
    tempGeno <- allDat %>% select(all_of(tempSNP)) %>% unlist(.)
    tempMAF <- mean(tempGeno, na.rm = TRUE) / 2
    
    # check for NA 
    NAidx <- which(is.na(tempGeno))
    if (length(NAidx) > 0) {
      tempGeno[NAidx] <- rbinom(n=length(NAidx), size=2, prob=tempMAF)
    }
    
    # all outcomes
    rawY <- allDat %>% select(all_of(exprNames)) %>%
      as.matrix(.)
    # filter out ones that don't vary
    lengthUnique <- function(x) {length(unique(x))}
    numUnique <- apply(rawY, 2, lengthUnique)
    yMat <- rawY[, which(numUnique > 10)]
    
    # design matrix
    xMat <- cbind(1, allDat %>% select(all_of(covarNames))) %>%
      as.matrix(.)
    projX <- xMat %*% solve(t(xMat) %*% xMat) %*% t(xMat)
    # residuals 
    fittedYMat <- projX %*% yMat
    residMat <- yMat - fittedYMat
    # estimated variance 
    sigmaSqYHat <- apply(residMat^2, 2, sum) / (nrow(xMat) - ncol(xMat))
    # variance-type quantity involving projection matrix 
    Ppart <- as.numeric(sum(tempGeno^2) - t(tempGeno) %*% projX %*% tempGeno)
    
    # one SNP against all the outcomes
    outcomeScoreStats <- as.numeric(t(tempGeno) %*% residMat)
    # can check against function checkMP()
    outcomeStandStats <- outcomeScoreStats / (sqrt(rep(Ppart, length(sigmaSqYHat)) * sigmaSqYHat))
    
    # record 
    tempResults <- data.frame(Gene=rep(NA, length(outcomeStandStats)), SNP=NA, MAF=NA, testStat=NA, pval=NA)
    tempResults$Gene <- names(outcomeStandStats)
    tempResults$SNP <- tempSNP
    tempResults$MAF <- tempMAF
    tempResults$testStat <- outcomeStandStats
    tempResults$pval <- 1 - pchisq(outcomeStandStats^2, df=1)
    
    # append
    resultsDF <- rbind(resultsDF, tempResults)
  }
  # return 
  return(list(individualDF = resultsDF, setbasedDF = setbasedDF))
}

#----------------------------------------------------_#
# start analysis

# the original GTEx data comes in some very complicated formats, we 
# keep the synthetic data much more clean
if (original) {
  # read from gds on disk
  retGds <- read_gds(genoLoc)
  genoData <- retGds$genoData
  
  # make SNP position df
  snpPosDF <- getSnpAnnotation(genoData) %>%
    pData() 
  
  # read expression metadata
  exprMeta <- fread(metaLoc)
  # read expression dataset - lung
  exprDF <- fread(exprLoc)
  # read covariates
  covarDF <- fread(covarLoc)
  
  # convert GTEx IDs to subject IDs
  subjIDs <- strsplit(exprMeta$Sample, "-") %>%
    sapply(function(x) {paste0(x[1], "-", x[2])})
  exprMeta$subj_id <- subjIDs
}

# read gene location information
load(file=geneInfoLoc)
geneRanges <-  data.table(ensembl_refgene_hg19_20180109)

# loop through the 18 loci given in Table 2 of McKay et al.
tab2DF <- data.frame(RS = c("rs71658797", "rs6920364", "rs11780471", "rs11571833", 
                            "rs66759488", "rs55781567", "rs56113850", "rs13080835",
                            "rs7705526", "rs4236709", "rs885518", "rs11591710", 
                            "rs1056562", "rs77468143", "rs41309931", "rs116822326",
                            "rs7953330", "rs17879961"),
                     Gene = c("FUBP1", "RNASET2", "CHRNA2", "BRCA2", "SEMA6D", 
                              "CHRNA5", "CYP2A6", "TP63", "TERT", "NRG1", "CDKN2A",
                              "OBFC1", "AMICA1", "SECISBP2L", "RTEL1", "HCP5", "RAD52", "CHEK2"),
                     BP = c(77967507, 167376466, 27344719, 32972626, 47577451, 78857986,
                            41353107, 189357199, 1285974, 32410110, 21830157, 105687632,
                            118125625, 49376624, 62326579, 31434111, 998819, 29121087),
                     Chr = c(1, 6, 8, 13, 15, 15, 19, 3, 5, 8, 9, 10, 11, 15, 20, 6, 12, 22))
individualDF <- c()
setbasedDF <- c()
locusToDo <- tab2DF$Gene
# outer loop through each SNP locus
for (snp_locus_it in 1:length(locusToDo)) {
  
  # one locus at a time 
  tempSNPlocus <- as.character(locusToDo[snp_locus_it])
  snpsToGet <- tab2DF %>% filter(Gene == tempSNPlocus) %>%
    mutate(start = BP - sentinelWindow, end = BP + sentinelWindow)
  
  # clean and merge expression, genotypes, covariates
  if (original) { 
    
    # gene expressions we want
    exprNames <- read.csv(twasLoc) %>%
      select(gene_name) %>%
      unlist(.)
    
    # clean data 
    lungCleaned <- clean_and_subset(exprGeneNames = exprNames, tissue = "Lung", genoData = genoData, 
                                    snpsToGet = snpsToGet, covarDF = covarDF, 
                                    exprMeta = exprMeta, exprDF = exprDF) 
    # found the SNP?
    if (class(lungCleaned) == "numeric") {next}
    
    # names of genes and SNPs
    exprNames <- colnames(lungCleaned$exprINT %>% select(-subj_id))
    genoNames <- colnames(lungCleaned$gMatSelected %>% select(-subj_id))
    # names of covariates
    covarNames <- unlist(covarDF$ID)
    # there is a gene named C2 and C3
    covarNames[which(covarNames == "C2")] <- "Cov2"
    covarNames[which(covarNames == "C3")] <- "Cov3"
    allDat <- lungCleaned$allDat
  } else {
    allDat <- fread(synthLoc)
    genoNames <- colnames(allDat)[2]
    exprNames <- colnames(allDat)[3:13571]
    covarNames <- colnames(allDat)[13572:13636]
  } 
  
  # checkpoint
  cat("Starting number ", snp_locus_it, "\n")
  cat("Number of genes: ", length(exprNames), "\n")
  cat("Number of SNPs: ", length(genoNames), "\n")
  
  # run eQTL analysis
  mpOutput <- MPreg(covarNames = covarNames, genoNames = genoNames, 
                    exprNames = exprNames, allDat = allDat)
  
  tempFname <- here::here(outputDir, paste0(snpsToGet$RS[1], "_", snpsToGet$Gene[1], "_", snpsToGet$Chr[1], "_", snpsToGet$BP[1], "_analysis.txt"))
  write.table(mpOutput$individualDF, tempFname, append=F, quote=F, row.names=F, col.names=T, sep='\t')
  
}

