# Figure 3D

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/Fig3D_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Lung/replication_simulation.R")

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
outName <- paste0(outputDir, "/rep_sim_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- FALSE
# these names are for if saveData <- TRUE
testStatsName <- here::here(outputDir, "rep_sim_allZ")


# simulation parameters
qvalue <- 0.1
nSNPs <- 6744757
K <- 2
nSims <- 100


# determines the proportion of four scenario 

sProp1 = c(0.915, 0.003, 0.003, 0.0393, 0.0393, 0.0001, 0.0001, 0.0001, 0.0001)
sProp2 = c(0.8, 0.03085, 0.03085, 0.06715, 0.06715, 0.001, 0.001, 0.001, 0.001)
sProp3 = c(0.6,0.07685,0.07685,0.11315,0.11315,0.005,0.005,0.005,0.005)

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
  arrange(blVec, Var1, Var2) %>%
  mutate(l = 0:(nrow(.) - 1)) %>%
  relocate(l, .before = symAlt) %>%
  cbind(sProp1) %>%
  cbind(sProp2) %>%
  cbind(sProp3) %>%
  select(Var1, Var2, sProp1, sProp2, sProp3)

number1 <- c()
number2 <- c()
number3 <- c()
for (s_it in 0:nrow(Hmat)) {
  number1 <- append(number1, as.integer(Hmat$sProp1[s_it] * nSNPs ))
  number2 <- append(number2, as.integer(Hmat$sProp2[s_it] * nSNPs )) 
  number3 <- append(number3, as.integer(Hmat$sProp3[s_it] * nSNPs )) 
}
Hmat <- Hmat %>% 
  mutate(number1 = number1) %>% 
  mutate(number2 = number2) %>% 
  mutate(number3 = number3)

muinfor <- fread(paste0(outputDir,"/rep_data_aID4_muInfo.txt")) %>%
  as.matrix()
Hmat$mu1 <- append(muinfor[1,], rep(muinfor[1,ncol(muinfor)],K^2-1))
Hmat$mu2 <- append(muinfor[2,], rep(muinfor[2,ncol(muinfor)],K^2-1))

for (case in 1:nrow(Hmat)) {
  if (Hmat$Var1[case] < 0){
    Hmat$mu1[case] = -Hmat$mu1[case]
  } 
  
  if(Hmat$Var2[case] < 0){
    Hmat$mu2[case] = -Hmat$mu2[case]
  } 
}


# record results here
powerRes <- data.frame(nCausal = rep(0,3), 
                       propCausal = rep(0,3), 
                       seed=NA, 
                       nRejNew=NA, nRejMeta=NA, nRejThreshold=NA,
                       powerNew=NA, powerMeta=NA, powerThreshold=NA, 
                       fdpNew=NA, fdpMeta=NA, fdpThreshold=NA)


for (prop_it in 1:3) {

# set the seed 
set.seed(aID * 10^3 + prop_it)
powerRes$seed[prop_it] <- aID * 10^3 + prop_it
col_name <- paste0("number", prop_it) 
sProp_name <- paste0("sProp", prop_it) 
powerRes$nCausal[prop_it] <- Hmat[,col_name][nrow(Hmat)]*2
powerRes$propCausal[prop_it] <- Hmat[,sProp_name ][nrow(Hmat)]*2


# load or save data
if (loadData) {
    setwd(outputDir)
    allZ <- fread(paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
  } else {
    
   allZ = c()
   for (case in 1:nrow(Hmat)) {
     tempZ = mvrnorm(n = Hmat[,col_name][case], mu = c(Hmat$mu1[case],Hmat$mu2[case]), Sigma = diag(K))
     allZ = rbind(allZ,tempZ)
   }
 

    # save it
  if (saveData) { 
    setwd(outputDir)
    write.table(allZ, paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    }
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

colnames(allZ) <- c("Z1","Z2")

###hold the results

# p-value matrix
allP <- 1- pchisq(as.matrix(allZ)^2, df=1)
colnames(allP) <- c("P1","P2")
testDat = as.data.frame(cbind(allZ,allP))

causalVec = c(rep(0,sum(Hmat[,col_name][c(1:(nrow(Hmat)-K^2))])),
              rep(1,Hmat[,col_name][(nrow(Hmat)-K^2) + 1]), 
              rep(0,Hmat[,col_name][nrow(Hmat)]*2),
              rep(1,Hmat[,col_name][nrow(Hmat)])
              )
totOut <- data.frame(Z1 = allZ[, 1], Z2 = allZ[, 2], P1 = allP[, 1], P2 = allP[, 2], origIdx = 1:nrow(allP), causal=causalVec) 



### teseting part

#First: threshold method 

same_dir <- subset(totOut, (totOut$Z1>0 & totOut$Z2>0) | (totOut$Z1<0 & totOut$Z2<0))
rejThres <- subset(same_dir, same_dir$P1 < 1e-8 & same_dir$P2 < 1e-8)

totOut <- totOut %>%
  mutate(threshold = ifelse(origIdx %in% rejThres$origIdx, 1, 0 ))

powerRes$fdpThreshold[prop_it] <- length(which(totOut$threshold == 1  & totOut$causal == 0)) /  length(which(totOut$threshold == 1))
powerRes$powerThreshold[prop_it] <- length(which(totOut$threshold == 1  & totOut$causal == 1)) / sum(causalVec)
powerRes$nRejThreshold[prop_it] <- nrow(rejThres)




#Second: new method
initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                       matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, sameDirAlt=TRUE, eps=10^(-5))

# save
# write.table(newRes$lfdrResults, paste0(fnameRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
# write.table(do.call(cbind, newRes$muInfo), paste0(fnameRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
# write.table(do.call(cbind, newRes$piInfo), paste0(fnameRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)

totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
  arrange(newLfdr) %>%
  mutate(newAvg = cummean(newLfdr)) %>%
  arrange(origIdx)
powerRes$fdpNew[prop_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$newAvg < 0.1))
powerRes$powerNew[prop_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
powerRes$nRejNew[prop_it] <- length(which(totOut$newAvg < 0.1))



#Third: Meta Analysis (Sample size based)
w1 = sqrt(85716) #McKay's
w2 = sqrt(332000) #UKB

totOut$tempZ1 = qnorm(totOut$P1/ 2) * sign(totOut$Z1)
totOut$tempZ2 = qnorm(totOut$P2/ 2) * sign(totOut$Z2)
totOut$Z_meta = (totOut$tempZ1*w1 + totOut$tempZ2*w2)/sqrt(w1^2 + w1^2)
totOut$P_meta =2*pnorm(abs(-totOut$Z_meta), lower.tail=FALSE)

rejMeta <- subset(totOut, totOut$P_meta < 1e-8)

totOut <- totOut %>%
  mutate(meta = ifelse(origIdx %in% rejMeta$origIdx, 1, 0 ))

powerRes$fdpMeta[prop_it] <- length(which(totOut$meta == 1 & totOut$causal == 0)) /  length(which(totOut$meta == 1 ))
powerRes$powerMeta[prop_it] <- length(which(totOut$meta == 1 & totOut$causal == 1)) / sum(causalVec)
powerRes$nRejMeta[prop_it] <- nrow(rejMeta)


 }
  
 
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')






