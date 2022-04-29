source("SEM_Gaussian.R")
source("EM_Gaussian.R")
source("EM_Gaussian_withoutNA.R")
library(VarSelLCM)
library(RMixtComp)
library(MASS)
library(mice)


n <- 100


K <- 3
pik <-  c(0.5,0.25,0.25)


source("SEM_Gaussian.R")
source("EM_Gaussian.R")
library(VarSelLCM)
library(RMixtComp)
library(MASS)


d <- 6 #9, 6, 3


tau <- 2.31


delta = matrix(0,nrow=K,ncol=d)
delta[1,1] = tau
delta[2,2] = tau
delta[3,3] = tau
delta[1,4] = tau
delta[2,5] = tau
delta[3,6] = tau
#delta[1,7] = tau
#delta[2,8] = tau
#delta[3,9] = tau
probmiss_z <- matrix(0,nrow=K,ncol=d,byrow=TRUE)
probmiss_y  <- rep(c(1.45,0.2,-3),2)


intercept_y <- -1.38


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

ARI_MCAR <- numeric(nb_it)
ARI_MNARz <- numeric(nb_it)
ARI_MNARy <- numeric(nb_it)
ARI_MNARz_conc <- numeric(nb_it)
ARI_MNAR_true <- numeric(nb_it)
ARI_Mice <- numeric(nb_it)


degen_MCAR <- numeric(nb_it)
degen_MNARz <- numeric(nb_it)
degen_MNARy <- numeric(nb_it)
degen_MNAR_true <- numeric(nb_it)
degen_Mice <- numeric(nb_it)

Nbit_run <- 3

#ARI_MNAR_true <- ARI_MNARy

for (it in 1:nb_it){
  print("ITERATION GENERALE")
  print(it)
  set.seed(it)
  Z <- SimuZ(n=n,pik=pik)
  Partition_true <- apply(Z, 1, function(z) which(z==1))
  Y <- matrix(NA,nrow=n,ncol=d)
  for (j in 1:d){
    Y[,j] <- Z%*%delta[,j] + rnorm(n) 
  }
  
  C <- SimuC(pik,Y,Z,"MNARy",probmiss_z=probmiss_z,probmiss_y=probmiss_y,intercept_y = intercept_y)
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

  
  #### EM MCAR
  print("MCAR")
  run_MCAR <- NULL
  it_MCAR <- 1
  while (it_MCAR<=Nbit_run){
    res_MCAR <- EM_Gaussian(YNA = YNA, K = K, mecha = "MCAR", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
    if (typeof(res_MCAR)=="character"  & degen_MCAR[it] < 10){
      degen_MCAR[it] <- degen_MCAR[it] + 1
      print(res_MCAR)
    }else{
      run_MCAR[[it_MCAR]] <- res_MCAR
      it_MCAR <- it_MCAR + 1
    }
    print(it_MCAR)
  }
  if (typeof(res_MCAR)=="character"){
    ARI_MCAR[it] <- NA
  }else{
    crit_MCAR <- lapply(run_MCAR,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MCAR")
    whmax_MCAR <- which.max(sapply(crit_MCAR,function(x) x$ICL))
    ARI_MCAR[it] <- crit_MCAR[[whmax_MCAR]]$ARI
  }



  #### EM MNARz
  print("MNARz")

  run_MNARz <- NULL
  it_MNARz <- 1
  while (it_MNARz<=Nbit_run & degen_MNARz[it] < 10){
    res_MNARz <- EM_Gaussian(YNA = YNA, K = K, mecha = "MNARz", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
    if (typeof(res_MNARz)=="character"){
      degen_MNARz[it] <- degen_MNARz[it] + 1
    }else{
      run_MNARz[[it_MNARz]] <- res_MNARz
      it_MNARz <- it_MNARz + 1
    }
    print(it_MNARz)
  }
  if (typeof(res_MNARz)=="character"){
    ARI_MNARz[it] <- NA
  }else{
    crit_MNARz <- lapply(run_MNARz,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MNARz")
    whmax_MNARz <- which.max(sapply(crit_MNARz,function(x) x$ICL))
    ARI_MNARz[it] <- crit_MNARz[[whmax_MNARz]]$ARI
  }
  
  
  
  #### True mecha
  print("True mecha")

  run_MNAR_true <- NULL
  it_MNAR_true <- 1
  while (it_MNAR_true<=Nbit_run & degen_MNAR_true[it] < 10){
    res_MNAR_true <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARy", diag = diag, rmax = rmax, init = NULL, stop = stop, samplesize = samplesize)
      #EM_Gaussian(YNA = YNA, K = K, mecha = "MNARzj", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
    if (typeof(res_MNAR_true)=="character"){
      degen_MNAR_true[it] <- degen_MNAR_true[it] + 1
    }else{
      run_MNAR_true[[it_MNAR_true]] <- res_MNAR_true
      it_MNAR_true <- it_MNAR_true + 1
    }
    print(it_MNAR_true)
  }
  if (typeof(res_MNAR_true)=="character"){
    ARI_MNAR_true[it] <- NA
  }else{
    crit_MNAR_true <- lapply(run_MNAR_true,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MNARykzj")
    whmax_MNAR_true <- which.max(sapply(crit_MNAR_true,function(x) x$ICL))
    ARI_MNAR_true[it] <- crit_MNAR_true[[whmax_MNAR_true]]$ARI
  }

  
  #### MICE
  
  print("MICE")

  nbimp_multmice <- 5
  res_multmice = mice(YNA,m=nbimp_multmice,maxit=50)
  ARI_Mice_vec <- numeric(nbimp_multmice)
  for (imult in 1:nbimp_multmice){
    Y_SingleMice = mice::complete(res_multmice, imult)
    run_singlemice <- "error"
    it_singlemice <- 1
    seed_singlemice <- 0
    while(typeof(run_singlemice)=="character" & degen_Mice[it] < 15){
      run_singlemice <-  EM_Gaussian_withoutNA(Y = data.matrix(Y_SingleMice), K = K, diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
      if (typeof(run_singlemice)=="character"){
        degen_Mice[it] <- degen_Mice[it] + 1
        seed_singlemice <- seed_singlemice + 1
      }
      it_singlemice <- it_singlemice + 1
      print(it_singlemice)
    }
    if (typeof(run_singlemice)!="character"){
      crit_MultMice = Critere_Gaussian_withoutNA(run_singlemice,Y_SingleMice,Partition_true)
      ARI_Mice_vec[imult] <- crit_MultMice$ARI
    }else{
      ARI_Mice_vec[imult] <- NA
    }

  }
  ARI_Mice[it] <- mean(ARI_Mice_vec,na.rm=TRUE)


  
  #### EM MNARy
  # print("MNARy")
  # 
  # run_MNARy <- NULL
  # it_MNARy <- 1
  # while (it_MNARy<=Nbit_run & degen_MNARy[it] < 10){
  #   res_MNARy <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARyzj", diag = diag, rmax = rmax_MNARyzj, init = NULL, stop = stop, samplesize = samplesize)
  #   if (typeof(res_MNARy)=="character"){
  #     degen_MNARy[it] <- degen_MNARy[it] + 1
  #   }else{
  #     run_MNARy[[it_MNARy]] <- res_MNARy
  #     it_MNARy <- it_MNARy + 1
  #   }
  #   print(it_MNARy)
  # }
  # Tdiff_MNARy[it]= difftime(T2, T1)
  # if (typeof(res_MNARy)=="character"){
  #   ARI_MNARy[it] <- NA
  # }else{
  #   crit_MNARy <- lapply(run_MNARy,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MNARyzj")
  #   whmax_MNARy <- which.max(sapply(crit_MNARy,function(x) x$ICL))
  #   ARI_MNARy[it] <- crit_MNARy[[whmax_MNARy]]$ARI
  # }
  
  
  # #### MNARz with (Y|C)
  # print("MNARz concatenate")
  # 
  # YNA_mixtcomp = YNA
  # YNA_mixtcomp[C] = "?"
  # YNA_mixtcomp <- cbind(YNA_mixtcomp,matrix(as.character(C+1),ncol=d,nrow=n))
  # colnames(YNA_mixtcomp) <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12")
  # algo <- list(nbBurnInIter = 100,
  #             nbIter = 1000,
  #             nbGibbsBurnInIter = 100,
  #             nbGibbsIter = 100,
  #             nInitPerClass = 20,
  #             nSemTry = 20,
  #             confidenceLevel = 0.95)
  # model <- list(X1 = "Gaussian", X2 = "Gaussian",
  #               X3 = "Gaussian", X4 = "Gaussian", X5 = "Gaussian", 
  #               X6 = "Gaussian",X7 = "Multinomial", X8 = "Multinomial",
  #               X9 = "Multinomial", X10 = "Multinomial", X11 = "Multinomial", X12 = "Multinomial")
  # resLearn <- mixtCompLearn(YNA_mixtcomp, algo=algo, model=model, nClass = K, nRun = 100)
  # ARI(Partition_true,getPartition(resLearn))
  # Tdiff_MNARz_comp[it]= difftime(T2, T1)
  # ARI_MNARz_conc[it] <- ARI(Partition_true,getPartition(resLearn))
  
}


save.image("xp2_MNARy_6.RData")

