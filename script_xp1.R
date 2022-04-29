source("SEM_Gaussian.R")
source("EM_Gaussian.R")
library(VarSelLCM)
library(RMixtComp)
library(MASS)

res <- c()

n <- 10^5
d <- 6

K <- 3
pik <-  c(0.5,0.25,0.25)

Z <- SimuZ(n=n,pik=pik)
Partition_true <- apply(Z, 1, function(z) which(z==1))

delta = matrix(0,nrow=K,ncol=d)
tau <- 3.3



#559
delta[1,1] = tau
delta[2,2] = tau
delta[3,3] = tau
delta[1,4] = tau
delta[2,5] = tau
delta[3,6] = tau
#delta[1,7] = tau
#delta[2,8] = tau
#delta[3,9] = tau

Y <- matrix(NA,nrow=n,ncol=d)
mu <- list(rep(0,d),rep(0,d),rep(0,d))
for (j in 1:d){
  Y[,j] <- Z%*%delta[,j] + rnorm(n)
  mu[[1]][j] <- (Z%*%delta[,j])[Z[,1]==1][1]
  mu[[2]][j] <- (Z%*%delta[,j])[Z[,2]==1][1]
  mu[[3]][j] <- (Z%*%delta[,j])[Z[,3]==1][1]
}
sigma <- list(diag(1,d),diag(1,d),diag(1,d))

###For xp3b
# for (j in 1:d){
#   #Y[,j] <- Z%*%delta[,j] + rnorm(n)
#   mu[[1]][j] <- (Z%*%delta[,j])[Z[,1]==1][1]
#   mu[[2]][j] <- (Z%*%delta[,j])[Z[,2]==1][1]
#   mu[[3]][j] <- (Z%*%delta[,j])[Z[,3]==1][1]
# }
# rho=0.5
# sigma_k <- lower.tri(diag(1,d),diag=FALSE)*rho+upper.tri(diag(1,d),diag=FALSE)*rho+diag(1,d)
# sigma <- list(sigma_k,sigma_k,sigma_k)
# for (k in 1:K){
#   obs_k <- which(Z[,k]==1)
#   Y[obs_k, ] <- mvrnorm(n*pik[k] , mu[[k]] , sigma[[k]])
# }


mecha = "MNARz"

#gamma <- 0.005
probmiss_z <- c(-0.4,0.23,0.7) 

#c(-1,-0.3,0) 
probmiss_y <- matrix(0,nrow=K,ncol=d,byrow=TRUE)
#probmiss_y  <- matrix(c(c(-3,0.3,-3,-3,-2,1),c(0.5,-2,1,1,1,0.5),c(1,1,0.5,0.5,0.5,2)),nrow=K,ncol=d,byrow=TRUE)
#intercept_y <- -0.75

intercept_y <- NULL

C <- SimuC(pik,Y,Z,"MNARz",probmiss_z=probmiss_z,probmiss_y=probmiss_y,intercept_y=intercept_y)
#SimuC_miss(pik,Y,Z,probmiss_z=probmiss_z,opt="Laplace")
#SimuC(pik,Y,Z,"MNARyk",probmiss_z=probmiss_z,probmiss_y=probmiss_y,intercept_y=intercept_y)
#SimuC_miss(pik,Y,Z,probmiss_z=probmiss_z,opt="Laplace")

YNA <- Y

YNA[C] <- NA

sum(is.na(YNA))/(n*d)
sum(is.na(YNA[,1]))/(n)
sum(is.na(YNA[,2]))/(n)
sum(is.na(YNA[,3]))/(n)
sum(is.na(YNA[,4]))/(n)
sum(is.na(YNA[,5]))/(n)
sum(is.na(YNA[,6]))/(n)
#sum(is.na(YNA[,7]))/(n)
#sum(is.na(YNA[,8]))/(n)
#sum(is.na(YNA[,9]))/(n)
sum(is.na(YNA[Z[,1]==1]))/(sum(Z[,1]==1)*d)
sum(is.na(YNA[Z[,2]==1]))/(sum(Z[,2]==1)*d)
sum(is.na(YNA[Z[,3]==1]))/(sum(Z[,3]==1)*d)


#
indexallNA <- c()
for (i in 1:n){
  if (sum(is.na(YNA[i,]))==d){
    indexallNA <- c(indexallNA,i)
    num <- sample(1:d,1)
    YNA[i,num] <- Y[i,num]
  }
}

sum(is.na(YNA))/(n*d)

restheo <- list(mu=mu,sigma=sigma,pik=pik,alpha=matrix(probmiss_z,nrow=K,ncol=d),beta=probmiss_y)
#res <- c(res,Critere_Gaussian(restheo,YNA,mecha="MNARy",Partition_true)$ARI)
Critere_Gaussian(restheo,YNA,mecha="MNARz",Partition_true)$ARI


#########################
d <- 6

K <- 3
pik <-  c(0.5,0.25,0.25)


source("SEM_Gaussian.R")
source("EM_Gaussian.R")
library(VarSelLCM)
library(RMixtComp)
library(MASS)


n = 100

tau = 1.92
delta = matrix(0,nrow=K,ncol=d)
delta[1,1] = tau
delta[2,2] = tau
delta[3,3] = tau
delta[1,4] = tau
delta[2,5] = tau
delta[3,6] = tau
probmiss_z <- matrix(0,nrow=K,ncol=d,byrow=TRUE)
probmiss_y  <- matrix(c(c(-3,0.3,-3,-3,-2,1),c(0.5,-2,1,1,1,0.5),c(1,1,0.5,0.5,0.5,2)),nrow=K,ncol=d,byrow=TRUE)
intercept_y <- -0.75



rmax_EM = 60
stop_EM = "classical"
tol_EM = 0.001
diag = TRUE

rmax = 60
rmax_MNARyzj = 30
init = NULL
stop = "loglikmax"
samplesize = NULL

Nbit_run <- 5

nb_it <- 50

ARI_MCAR <- numeric(nb_it)
ARI_MNARz <- numeric(nb_it)
ARI_MNARy <- numeric(nb_it)
ARI_MNARz_conc <- numeric(nb_it)


degen_MCAR <- numeric(nb_it)
degen_MNARz <- numeric(nb_it)
degen_MNARy <- numeric(nb_it)

Tdiff_MCAR <- numeric(nb_it)
Tdiff_MNARz <- numeric(nb_it)
Tdiff_MNARy <- numeric(nb_it)
Tdiff_MNARz_comp <- numeric(nb_it)


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
  
  C <- SimuC(pik,Y,Z,"MNARyk",probmiss_z=probmiss_z,probmiss_y=probmiss_y,intercept_y = intercept_y)
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
  T1<-Sys.time() 
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
  T2<-Sys.time() 
  Tdiff_MCAR[it]= difftime(T2, T1)
  if (typeof(res_MCAR)=="character"){
    ARI_MCAR[it] <- NA
  }else{
    crit_MCAR <- lapply(run_MCAR,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MCAR")
    whmax_MCAR <- which.max(sapply(crit_MCAR,function(x) x$ICL))
    ARI_MCAR[it] <- crit_MCAR[[whmax_MCAR]]$ARI
  }
  
  
  
  #### EM MNARz ou EM MNARzj
  # print("MNARzj")
  # 
  # T1<-Sys.time()
  # run_MNARz <- NULL
  # it_MNARz <- 1
  # while (it_MNARz<=Nbit_run & degen_MNARz[it] < 10){
  #   res_MNARz <- EM_Gaussian(YNA = YNA, K = K, mecha = "MNARzj", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
  #   if (typeof(res_MNARz)=="character"){
  #     degen_MNARz[it] <- degen_MNARz[it] + 1
  #   }else{
  #     run_MNARz[[it_MNARz]] <- res_MNARz
  #     it_MNARz <- it_MNARz + 1
  #   }
  #   print(it_MNARz)
  # }
  # T2<-Sys.time()
  # Tdiff_MNARz[it]= difftime(T2, T1)
  # if (typeof(res_MNARz)=="character"){
  #   ARI_MNARz[it] <- NA
  # }else{
  #   crit_MNARz <- lapply(run_MNARz,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MNARzj")
  #   whmax_MNARz <- which.max(sapply(crit_MNARz,function(x) x$ICL))
  #   ARI_MNARz[it] <- crit_MNARz[[whmax_MNARz]]$ARI
  # }
  # 
  #### EM MNARy
  print("MNARy")
  
  T1<-Sys.time()
  run_MNARy <- NULL
  it_MNARy <- 1
  while (it_MNARy<=Nbit_run & degen_MNARy[it] < 10){
    res_MNARy <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARyk", diag = diag, rmax = rmax, init = NULL, stop = stop, samplesize = samplesize)
    if (typeof(res_MNARy)=="character"){
      degen_MNARy[it] <- degen_MNARy[it] + 1
    }else{
      run_MNARy[[it_MNARy]] <- res_MNARy
      it_MNARy <- it_MNARy + 1
    }
    print(it_MNARy)
  }
  # res_MNARy <- "error"
  # while(typeof(res_MNARy)=="character" & degen_MNARy[it] < 10){
  #   res_MNARy <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARyk", diag = diag, rmax = rmax, init = run_MCAR[[whmax_MCAR]], stop = stop, samplesize = samplesize)
  #   if (typeof(res_MNARy)=="character"){
  #     degen_MNARy[it] <- degen_MNARy[it] + 1
  #   }else{
  #     run_MNARy <- res_MNARy
  #   }
  #   print(degen_MNARy[it])
  # }
  T2<-Sys.time()
  Tdiff_MNARy[it]= difftime(T2, T1)
  if (typeof(res_MNARy)=="character"){
    ARI_MNARy[it] <- NA
  }else{
    crit_MNARy <- lapply(run_MNARy,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true,mecha="MNARyk")
    whmax_MNARy <- which.max(sapply(crit_MNARy,function(x) x$ICL))
    ARI_MNARy[it] <- crit_MNARy[[whmax_MNARy]]$ARI
  }
  # if (typeof(res_MNARy)=="character"){
  #   ARI_MNARy[it] <- NA
  # }else{
  #   crit_MNARy <- Critere_Gaussian(run_MNARy,YNA=YNA,Partition_true=Partition_true,mecha="MNARy")
  #   ARI_MNARy[it] <- crit_MNARy$ARI 
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
  # T1<-Sys.time() 
  # resLearn <- mixtCompLearn(YNA_mixtcomp, algo=algo, model=model, nClass = K, nRun = 100)
  # T2<-Sys.time() 
  # ARI(Partition_true,getPartition(resLearn))
  # Tdiff_MNARz_comp[it]= difftime(T2, T1)
  # ARI_MNARz_conc[it] <- ARI(Partition_true,getPartition(resLearn))
  
}


save.image("xp1_MNARyk_100.RData")

