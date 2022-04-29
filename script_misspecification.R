source("SEM_Gaussian.R")
source('EM_Gaussian.R')
library("MASS")
library('ExtDist')



n <- 10000
d <- 6
mu <- list(rep(0,d),rep(3,d),rep(6,d)) 
sigma <- list(diag(1,d)+matrix(0.8,ncol=d,nrow=d),diag(1,d)+matrix(0.8,ncol=d,nrow=d),diag(1,d)+matrix(0.8,ncol=d,nrow=d))

K <- 3
pik <-  c(0.5,0.2,0.3)

true_mecha = "MNARz"

diag = TRUE


#mu <- c(0,4,8)
#b <- c(1,1,1)

rmax = 60
init = NULL
stop = "loglikmax"
samplesize = NULL

rmax_EM = 60
stop_EM = "classical"
tol_EM = 0.001

Nbit_run <- 5


#probmiss_z <- c(-2,-1,-0.5)
#probmiss_y  =  c(-0.3,-10)
###PAS ASSEZ DE NA ?

probmiss_z <- c(-1,-0.3,1)


test_K <- c(1,2,3,4)
it_K <- 1:10



ICL_MCAR <- matrix(NA,ncol=length(it_K),nrow=length(test_K))
ICL_MNARz <- matrix(NA,ncol=length(it_K),nrow=length(test_K))
ICL_MNARy <- matrix(NA,ncol=length(it_K),nrow=length(test_K))
ICL_MNARyz <- matrix(NA,ncol=length(it_K),nrow=length(test_K))

ARI_MCAR <- matrix(NA,ncol=length(it_K),nrow=length(test_K))
ARI_MNARz <- matrix(NA,ncol=length(it_K),nrow=length(test_K))
ARI_MNARy <- matrix(NA,ncol=length(it_K),nrow=length(test_K))
ARI_MNARyz <- matrix(NA,ncol=length(it_K),nrow=length(test_K))

degen_MCAR <- matrix(0,ncol=length(it_K),nrow=length(test_K))
degen_MNARz <- matrix(0,ncol=length(it_K),nrow=length(test_K))
degen_MNARy <- matrix(0,ncol=length(it_K),nrow=length(test_K))
degen_MNARyz <- matrix(0,ncol=length(it_K),nrow=length(test_K))


C <- NULL
Z <- NULL
Y <- NULL
Partition_true <- NULL

for (it in it_K){
  
  
  T1 <- Sys.time() 
  
  print(paste0("iteration ",it))
  
  set.seed(it)
  
  Z[[it]] <- SimuZ(n=n,pik=pik)
  
  Partition_true[[it]] <- apply(Z[[it]], 1, function(z) which(z==1))
  
  Y[[it]] <- SimuY(Z=Z[[it]],d=d,pik=pik,mu=mu,sigma=sigma)
  
  #C_it <- Simulation_MNARyall(pik,Y[[it]],Z[[it]],probmiss_y=probmiss_y,mecha=true_mecha,probmiss_z=probmiss_z)
  C_it <- SimuC(pik,Y[[it]],Z[[it]],probmiss_z=probmiss_z,mecha=true_mecha)
  YNA <- Y[[it]]
  YNA[C_it] <- NA
  
  for (i in 1:n){
    if (sum(is.na(YNA[i,]))==d){
      num <- sample(1:d,1)
      YNA[i,num] <- Y[[it]][i,num]
    }
  }
  
  C[[it]] <- is.na(YNA)
  
  for (K in c(3)){
    
    print(K)
    
    
    Tdiff <- numeric(length(it_K))
    
    print("MCAR")
    
    run_MCAR <- NULL
    it_MCAR <- 1
    while (it_MCAR<=Nbit_run){
      res_MCAR <- EM_Gaussian(YNA = YNA, K = K, mecha = "MCAR", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
      if (typeof(res_MCAR)=="character"  & degen_MCAR[K,it] < 100){
        degen_MCAR[K,it] <- degen_MCAR[K,it] + 1
        print(res_MCAR)
      }else{
        run_MCAR[[it_MCAR]] <- res_MCAR
        it_MCAR <- it_MCAR + 1
      }
      print(it_MCAR)
    }
    if (typeof(res_MCAR)=="character"){
      ICL_MCAR[K,it] <- NA
      ARI_MCAR[K,it] <- NA
    }else{
      crit_MCAR <- lapply(run_MCAR,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true[[it]],mecha="MCAR")
      whmax_MCAR <- which.max(sapply(crit_MCAR,function(x) x$ICL))
      ICL_MCAR[K,it] <- max(sapply(crit_MCAR,function(x) x$ICL)) 
      ARI_MCAR[K,it] <- crit_MCAR[[whmax_MCAR]]$ARI
    }
    
    print("MNARz")
    
    run_MNARz <- NULL
    it_MNARz <- 1
    while (it_MNARz<=Nbit_run & degen_MNARz[K,it] < 30){
      res_MNARz <- EM_Gaussian(YNA = YNA, K = K, mecha = "MNARz", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
      if (typeof(res_MNARz)=="character"){
        degen_MNARz[K,it] <- degen_MNARz[K,it] + 1
      }else{
        run_MNARz[[it_MNARz]] <- res_MNARz
        it_MNARz <- it_MNARz + 1
      }
      print(it_MNARz)
    }
    if (typeof(res_MNARz)=="character"){
      ICL_MNARz[K,it] <- NA
      ARI_MNARz[K,it] <- NA
    }else{
      crit_MNARz <- lapply(run_MNARz,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true[[it]],mecha="MNARz")
      whmax_MNARz <- which.max(sapply(crit_MNARz,function(x) x$ICL))
      ICL_MNARz[K,it] <- max(sapply(crit_MNARz,function(x) x$ICL)) 
      ARI_MNARz[K,it] <- crit_MNARz[[whmax_MNARz]]$ARI
    }
    
    # res_MNARz <- "error"
    # while(typeof(res_MNARz)=="character" & degen_MNARz[K,it] < 10){
    #   res_MNARz <- EM_Gaussian(YNA = YNA, K = K, mecha = "MNARz", diag = diag, rmax = rmax_EM, init = run_MCAR[[whmax_MCAR]], stop = stop_EM, tol = tol_EM)
    #   if (typeof(res_MNARz)=="character"){
    #     degen_MNARz[K,it] <- degen_MNARz[K,it] + 1
    #   }else{
    #     run_MNARz <- res_MNARz
    #   }
    #   print(degen_MNARz[K,it])
    # }
    # if (typeof(res_MNARz)=="character"){
    #   while(typeof(res_MNARz)=="character" & degen_MNARz[K,it] < 20){
    #     res_MNARz <- EM_Gaussian(YNA = YNA, K = K, mecha = "MNARz", diag = diag, rmax = rmax_EM, init = NULL, stop = stop_EM, tol = tol_EM)
    #     if (typeof(res_MNARz)=="character"){
    #       degen_MNARz[K,it] <- degen_MNARz[K,it] + 1
    #     }else{
    #       run_MNARz <- res_MNARz
    #     }
    #     print(degen_MNARz[K,it])
    #   }
    # }
    # if (typeof(res_MNARz)=="character"){
    #   ICL_MNARz[K,it] <- NA
    # }else{
    #   crit_MNARz <- Critere_Gaussian(run_MNARz,YNA=YNA,Partition_true=Partition_true[[it]],mecha="MNARz")
    #   ICL_MNARz[K,it] <- crit_MNARz$ICL 
    #   ARI_MNARz[K,it] <- crit_MNARz$ARI 
    # }
    # 
    
    print("MNARy")
    
    run_MNARy <- NULL
    it_MNARy <- 1
    while (it_MNARy<=Nbit_run & degen_MNARy[K,it] < 30){
      res_MNARy <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARy", diag = diag, rmax = rmax, init = NULL, stop = stop, samplesize = samplesize)
      if (typeof(res_MNARy)=="character"){
        degen_MNARy[K,it] <- degen_MNARy[K,it] + 1
      }else{
        run_MNARy[[it_MNARy]] <- res_MNARy
        it_MNARy <- it_MNARy + 1
      }
      print(it_MNARy)
    }
    if (typeof(res_MNARy)=="character"){
      ICL_MNARy[K,it] <- NA
      ARI_MNARy[K,it] <- NA
    }else{
      crit_MNARy <- lapply(run_MNARy,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true[[it]],mecha="MNARy")
      whmax_MNARy <- which.max(sapply(crit_MNARy,function(x) x$ICL))
      ICL_MNARy[K,it] <- max(sapply(crit_MNARy,function(x) x$ICL)) 
      ARI_MNARy[K,it] <- crit_MNARy[[whmax_MNARy]]$ARI
    }
    
   
    print("MNARyz")
    
    
    res_MNARyz <- "error"
    while(typeof(res_MNARyz)=="character" & degen_MNARyz[K,it] < 10){
      res_MNARyz <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARyz", diag = diag, rmax = rmax, init = run_MNARy[[whmax_MNARy]], stop = stop, samplesize = samplesize)
      if (typeof(res_MNARyz)=="character"){
        degen_MNARyz[K,it] <- degen_MNARyz[K,it] + 1
      }else{
        run_MNARyz <- res_MNARyz
      }
      print(degen_MNARyz[K,it])
    }
    if (typeof(res_MNARyz)=="character"){
      while(typeof(res_MNARyz)=="character" & degen_MNARyz[K,it] < 20){
        res_MNARyz <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARyz", diag = diag, rmax = rmax, init = NULL, stop = stop, samplesize = samplesize)
        if (typeof(res_MNARyz)=="character"){
          degen_MNARyz[K,it] <- degen_MNARyz[K,it] + 1
        }else{
          run_MNARyz <- res_MNARyz
        }
        print(degen_MNARyz[K,it])
      }
    }
    if (typeof(res_MNARyz)=="character"){
      ICL_MNARyz[K,it] <- NA
    }else{
      crit_MNARyz <- Critere_Gaussian(run_MNARyz,YNA=YNA,Partition_true=Partition_true[[it]],mecha="MNARyz")
      ICL_MNARyz[K,it] <- crit_MNARyz$ICL 
      ARI_MNARyz[K,it] <- crit_MNARyz$ARI 
    }
    
    # run_MNARyz <- NULL
    # it_MNARyz <- 1
    # while (it_MNARyz<=Nbit_run & degen_MNARyz[K,it] < 30){
    #   res_MNARyz <- SEM_Gaussian(YNA = YNA, K = K, mecha = "MNARyz", diag = diag, rmax = rmax, init = run_MNARy[[whmax_MNARy]], stop = stop, samplesize = samplesize)
    #   if (typeof(res_MNARyz)=="character"){
    #     degen_MNARyz[K,it] <- degen_MNARyz[K,it] + 1
    #   }else{
    #     run_MNARyz[[it_MNARyz]] <- res_MNARyz
    #     it_MNARyz <- it_MNARyz + 1
    #   }
    #   print(it_MNARyz)
    # }
    # if (typeof(res_MNARyz)=="character"){
    #   ICL_MNARyz[K,it] <- NA
    #   ARI_MNARyz[K,it] <- NA
    # }else{
    #   crit_MNARyz <- lapply(run_MNARyz,Critere_Gaussian,YNA=YNA,Partition_true=Partition_true[[it]],mecha="MNARyz")
    #   whmax_MNARyz <- which.max(sapply(crit_MNARyz,function(x) x$ICL))
    #   ICL_MNARyz[K,it] <- max(sapply(crit_MNARyz,function(x) x$ICL)) 
    #   ARI_MNARyz[K,it] <- crit_MNARy[[whmax_MNARyz]]$ARI
    # }
    # 
    
    T2 <- Sys.time()
    
    Tdiff[it] <- difftime(T2, T1) 
    
    print(paste0("difftime ",Tdiff[it]))

  }
  
  save.image(paste0("n500_missdata",it,".Rdata"))
  
  
}


