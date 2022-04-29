source("SEM_Gaussian.R")
source("EM_Gaussian.R")
source("EM_Gaussian_withoutNA.R")
source("Functions_misspecification.R")
library(VarSelLCM)
library(RMixtComp)
library(MASS)
library(mice)


n <- 100


K_true <- 3
pik <-  c(0.5,0.25,0.25)


source("SEM_Gaussian.R")
source("EM_Gaussian.R")
library(VarSelLCM)
library(RMixtComp)
library(MASS)


d <- 6 #9, 6, 3

tau <- 3.3




delta = matrix(0,nrow=K_true,ncol=d)
delta[1,1] = tau
delta[2,2] = tau
delta[3,3] = tau
delta[1,4] = tau
delta[2,5] = tau
delta[3,6] = tau

probmiss_z <- c(-0.55,0.25,1.7)

probmiss_y  <- matrix(0,nrow=K_true,ncol=d,byrow=TRUE)


intercept_y <- NULL


rmax_EM = 60
stop_EM = "classical"
tol_EM = 0.001
diag = TRUE

rmax = 60
rmax_MNARyzj = 30
rmax_MNARykzj = 5
init = NULL
stop = "loglikmax"
samplesize = NULL

Nbit_run <- 3

nb_it <- 50

#load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykzj_100.RData")

ARI_MCAR <- matrix(NA,nrow=nb_it,ncol=4)
ARI_MNARz <- matrix(NA,nrow=nb_it,ncol=4)

degen_MCAR <- matrix(0,nrow=nb_it,ncol=4)
degen_MNARz <- matrix(0,nrow=nb_it,ncol=4)

ICL_MCAR <- matrix(NA,nrow=nb_it,ncol=4)
ICL_MNARz <- matrix(NA,nrow=nb_it,ncol=4)

Nbit_run <- 3

#ARI_MNAR_true <- ARI_MNARy

for (it in 1:50){
  print("ITERATION GENERALE")
  print(it)
  set.seed(it)
  Z <- SimuZ(n=n,pik=pik)
  Partition_true <- apply(Z, 1, function(z) which(z==1))
  Y <- matrix(NA,nrow=n,ncol=d)
  for (j in 1:d){
    Y[,j] <- Z%*%delta[,j] + rnorm(n) 
  }
  
  C <- SimuC(pik,Y,Z,"MNARz",probmiss_z=probmiss_z,probmiss_y=probmiss_y,intercept_y=intercept_y)
  YNA <- Y
  YNA[C] <- NA
  
  indexallNA <- c()
  for (i in 1:n){
    if (sum(is.na(YNA[i,]))==d){
      indexallNA <- c(indexallNA,i)
      num <- sample(1:d,1)
      YNA[i,num] <- Y[i,num]
    }
  }
  
  print(sum(is.na(YNA))/(n*d))
  
  
  #### EM MCAR
  for (K in c(1,2,3,4)){
    print("MCAR")
    print(paste0("K is ",K))
    run_MCAR <- NULL
    it_MCAR <- 1
    while (it_MCAR<=Nbit_run){
      res_MCAR <- EM_Gaussian(YNA = YNA, K = K, mecha = "MCAR", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
      if (typeof(res_MCAR)=="character"  & degen_MCAR[it,K] < 10){
        degen_MCAR[it,K] <- degen_MCAR[it,K] + 1
        print(res_MCAR)
      }else{
        run_MCAR[[it_MCAR]] <- res_MCAR
        it_MCAR <- it_MCAR + 1
      }
      print(it_MCAR)
    }
    if (typeof(res_MCAR)=="character"){
      ARI_MCAR[it,K] <- NA
      ICL_MCAR[it,K] <- NA
    }else{
      crit_MCAR <- lapply(run_MCAR,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MCAR")
      whmax_MCAR <- which.max(sapply(crit_MCAR,function(x) x$ICL))
      ARI_MCAR[it,K] <- crit_MCAR[[whmax_MCAR]]$ARI
      ICL_MCAR[it,K] <- crit_MCAR[[whmax_MCAR]]$ICL
    }
  }
  
  
  
  #### EM MNARz
  for (K in c(1,2,3,4)){
    print("MNARz")
    print(paste0("K is ",K))
    run_MNARz <- NULL
    it_MNARz <- 1
    while (it_MNARz<=Nbit_run & degen_MNARz[it,K] < 10){
      res_MNARz <- EM_Gaussian(YNA = YNA, K = K, mecha = "MNARz", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
      if (typeof(res_MNARz)=="character"){
        degen_MNARz[it,K] <- degen_MNARz[it,K] + 1
      }else{
        run_MNARz[[it_MNARz]] <- res_MNARz
        it_MNARz <- it_MNARz + 1
      }
      print(it_MNARz)
    }
    if (typeof(res_MNARz)=="character"){
      ARI_MNARz[it,K] <- NA
      ICL_MNARz[it,K] <- NA
    }else{
      crit_MNARz <- lapply(run_MNARz,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MNARz")
      whmax_MNARz <- which.max(sapply(crit_MNARz,function(x) x$ICL))
      ARI_MNARz[it,K] <- crit_MNARz[[whmax_MNARz]]$ARI
      ICL_MNARz[it,K] <- crit_MNARz[[whmax_MNARz]]$ICL
    }
  }
  
  
}


save.image("xp4b_50perc.RData")

