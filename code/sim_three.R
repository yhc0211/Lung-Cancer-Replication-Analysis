# 3D simulation proportion of causal is 0.0002

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/Fig3D_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.

#threshold1: 1e-8
#threshold2: 1e-5
#threshold3: 1e-4
here::i_am("Lung/sim_three.R")

# load libraries
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
library(ks)
library(csmGmm)
library(here)


# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Lung", "output")
outName <- paste0(outputDir, "/rep_three_sim_Snum",Snum,"_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- FALSE
# these names are for if saveData <- TRUE
testStatsName <- here::here(outputDir, "rep_sim_allZ")


# simulation parameters
qvalue <- 0.1
case = Snum*27
nSNPs <- 6744757
K <- 3



# determines the proportion of four scenario 

sProp1 = c(0.9002, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 
          0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0055, 0.0001, 0.0001, 
          0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001)
sProp2 = c(0.794, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 
          0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.001, 0.001, 
          0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
sProp3 = c(0.6, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
          0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.005, 0.005, 0.005, 
          0.005, 0.005, 0.005, 0.005, 0.005)

sProp <- c(sProp1,sProp2,sProp3)
#make a Hmat
Hmat <- expand.grid(rep(list(-1:1), K))

# attach the bl
blVec <- rep(0, nrow(Hmat))
slVec <- rep(0, nrow(Hmat))
for (k_it in K:1) {
  blVec <- blVec + 2^(K - k_it) * abs(Hmat[, k_it])
  slVec <- slVec + abs(Hmat[, k_it])
}
# symmetric alternative
symAltVec <- ifelse(apply(Hmat[, 1:K], 1, sum) == K | apply(Hmat[, 1:K], 1, sum) == -K, 1, 0)
# sort Hmat
Hmat <- Hmat %>% mutate(bl = blVec) %>%
  mutate(sl = slVec) %>%
  mutate(symAlt = symAltVec) %>%
  arrange(blVec, Var1, Var2, Var3) %>%
  mutate(l = 0:(nrow(.) - 1)) %>%
  relocate(l, .before = symAlt) 

Hmat <- Hmat %>% 
  bind_rows(Hmat, Hmat) %>%
  cbind(sProp) %>%
  select(Var1, Var2, Var3, sProp)




mu1 <- c(Hmat$Var1 * 2,Hmat$Var1 * 3,Hmat$Var1 * 4, Hmat$Var1 * 5,Hmat$Var1 * 6)
mu2 <- c(Hmat$Var2 * 2,Hmat$Var2 * 3,Hmat$Var2 * 4, Hmat$Var2 * 5,Hmat$Var2 * 6)
mu3 <- c(Hmat$Var3 * 2,Hmat$Var3 * 3,Hmat$Var3 * 4, Hmat$Var3 * 5,Hmat$Var3 * 6)

Hmat <- Hmat %>% 
  bind_rows(Hmat,Hmat,Hmat,Hmat) %>%
  bind_cols(mu1, mu2, mu3) %>%
  rename("mu1" = "...5", "mu2" = "...6", "mu3" = "...7")

Hmat <- Hmat %>% 
  mutate(number = as.integer(nSNPs * Hmat$sProp)) 




Hmat <- Hmat[(case-26):case,]
  

# record results here
powerRes <- data.frame(nCausal = rep(0,1), 
                       propCausal = rep(0,1), 
                       mu = rep(0,1), 
                       seed=NA, 
                       nRejNew=NA, nRejMeta=NA, nRejThreshold1=NA, nRejThreshold2=NA, nRejThreshold3=NA,
                       powerNew=NA, powerMeta=NA, powerThreshold1=NA, powerThreshold2=NA, powerThreshold3=NA,
                       fdpNew=NA, fdpMeta=NA, fdpThreshold1=NA, fdpThreshold2=NA, fdpThreshold3=NA)



  
  # set the seed 
  set.seed(aID * 10^3)
  powerRes$seed <- aID * 10^3 
  powerRes$nCausal <- Hmat$number[nrow(Hmat)]*2
  powerRes$propCausal <- Hmat$sProp[nrow(Hmat)]*2
  powerRes$mu <- Hmat$mu1[nrow(Hmat)]
  
  

    
  allZ = c()
  for (i in 1:nrow(Hmat)) {
      tempZ = mvrnorm(n = Hmat[,"number"][i], mu = c(Hmat[,"mu1"][i],Hmat[,"mu2"][i],Hmat[,"mu3"][i]), Sigma = diag(K))
      allZ = rbind(allZ,tempZ)
    }
    
    

  
  # adjustment to not get p-values of 0 needed for DACT and HDMT 
  for (col_it in 1:ncol(allZ)) {
    tooBig <- which(allZ[, col_it] > 8.1)
    tooSmall <- which(allZ[, col_it] < -8.1)
    if (length(tooBig) > 0) {
      allZ[tooBig, col_it] <- 8.1
    }
    if (length(tooSmall) > 0) {
      allZ[tooSmall, col_it] <- -8.1
    }
  }
  
  colnames(allZ) <- c("Z1","Z2","Z3")
  
  ###hold the results
  causalVec = c(rep(0,sum(Hmat[,'number'][c(1:19)])),
                rep(1,Hmat[,'number'][20]), 
                rep(0,sum(Hmat[,'number'][c(21:26)])),
                rep(1,Hmat[,'number'][27])
  )
  
  # p-value matrix
  allP <- 1- pchisq(as.matrix(allZ)^2, df=1)
  totOut <- data.frame(Z1 = allZ[, 1], Z2 = allZ[, 2], Z3 = allZ[,3], P1 = allP[,1], P2 = allP[,2], P3 = allP[,3], origIdx = 1:nrow(allP), causal=causalVec) 
  
  
  
  ### teseting part
  
  #First: threshold method 
  
  same_dir <- subset(totOut, (totOut$Z1>0 & totOut$Z2>0 & totOut$Z3>0) | (totOut$Z1<0 & totOut$Z2<0 & totOut$Z3<0))
  rejThres1 <- subset(same_dir, same_dir$P1 < 1e-8 & same_dir$P2 < 1e-8 & same_dir$P3 < 1e-8)
  rejThres2 <- subset(same_dir, same_dir$P1 < 1e-6 & same_dir$P2 < 1e-6 & same_dir$P3 < 1e-6)
  rejThres3 <- subset(same_dir, same_dir$P1 < 1e-5 & same_dir$P2 < 1e-5 & same_dir$P3 < 1e-5)
  
  totOut <- totOut %>%
    mutate(threshold1 = ifelse(origIdx %in% rejThres1$origIdx, 1, 0 )) %>%
    mutate(threshold2 = ifelse(origIdx %in% rejThres2$origIdx, 1, 0 )) %>%
    mutate(threshold3 = ifelse(origIdx %in% rejThres3$origIdx, 1, 0 ))  
  
  powerRes$fdpThreshold1 <- length(which(totOut$threshold1 == 1  & totOut$causal == 0)) /  length(which(totOut$threshold1 == 1))
  powerRes$powerThreshold1 <- length(which(totOut$threshold1 == 1  & totOut$causal == 1)) / sum(causalVec)
  powerRes$nRejThreshold1 <- nrow(rejThres1)
  
  powerRes$fdpThreshold2 <- length(which(totOut$threshold2 == 1  & totOut$causal == 0)) /  length(which(totOut$threshold2 == 1))
  powerRes$powerThreshold2 <- length(which(totOut$threshold2 == 1  & totOut$causal == 1)) / sum(causalVec)
  powerRes$nRejThreshold2 <- nrow(rejThres2)
  
  powerRes$fdpThreshold3 <- length(which(totOut$threshold3 == 1  & totOut$causal == 0)) /  length(which(totOut$threshold3 == 1))
  powerRes$powerThreshold3 <- length(which(totOut$threshold3 == 1  & totOut$causal == 1)) / sum(causalVec)
  powerRes$nRejThreshold3 <- nrow(rejThres3)
  
  
  #Second: new method
  initPiList <- list(c(0.82))
  for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
  initPiList[[8]] <- c(0.1)
  tempH <- expand.grid(c(0, 1), c(0, 1), c(0, 1)) %>%
    mutate(s = Var1 + Var2 + Var3) %>%
    arrange(s) %>%
    select(-s) %>%
    as.matrix(.)
  initMuList <- list(matrix(data=0, nrow=3, ncol=1))
  for (i in 2:7) {
    initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
  }
  initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
  newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, sameDirAlt=TRUE, eps=10^(-5))


  totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
    arrange(newLfdr) %>%
    mutate(newAvg = cummean(newLfdr)) %>%
    arrange(origIdx)

  powerRes$fdpNew <- length(which(totOut$newAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$newAvg < 0.1))
  powerRes$powerNew <- length(which(totOut$newAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
  powerRes$nRejNew <- length(which(totOut$newAvg < 0.1))



  #Third: Meta Analysis (Sample size based)
  w1 = sqrt(85716) #overall
  w2 = sqrt(332422) #ukb
  w3 = sqrt(73106) #mvp
  
  
  totOut$Z_meta = (totOut$Z1*w1 + totOut$Z2*w2 + totOut$Z3*w3)/sqrt(w1^2 + w2^2 + w3^2)
  totOut$P_meta =2*pnorm(abs(-totOut$Z_meta), lower.tail=FALSE)
  
  rejMeta <- subset(totOut, totOut$P_meta < 1e-8)
  
  totOut <- totOut %>%
    mutate(meta = ifelse(origIdx %in% rejMeta$origIdx, 1, 0))
  
  powerRes$fdpMeta <- length(which(totOut$meta == 1 & totOut$causal == 0)) /  length(which(totOut$meta == 1 ))
  powerRes$powerMeta <- length(which(totOut$meta == 1 & totOut$causal == 1)) / sum(causalVec)
  powerRes$nRejMeta <- nrow(rejMeta)
  
  



write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')
