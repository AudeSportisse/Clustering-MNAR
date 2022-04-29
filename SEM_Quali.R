source('Internal_Functions_Quali.R')
source('Internal_Functions_All.R')
source('initialization.R')

######
# Name: SEM_Quali
# Date: 2022/17/02
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARzj" (character)
## rmax: maximum number of iterations (numeric, default is 100)
## init: choice of the initialization  (NULL, "sample" or list)
## stop: choice of stop criterium "cut_average", "loglikmax" or "rmax"
## samplesize: size of the data sample used for the initialization, if init="sample" (numeric)
#####


SEM_Quali <- function(YNA, K, mecha, rmax = 100, init = NULL, stop = "rmax", samplesize = NULL) {
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)), ncol = ncol(YNA))
  
  Nb_levels <- sapply(1:d, function(j) length(unique(YNA[!is.na(YNA[,j]),j])))
  
  init <- init_SEM_Quali(YNA, K, mecha, init, samplesize)
  Z.new <- init$Z.init
  pi.new <- init$pi.init
  prob.new <- init$prob.init
  alpha.new <- init$alpha.init
  
  if (stop == "loglikmax"){
    pi.new.save <- pi.new
    prob.new.save <- prob.new
    alpha.new.save <- alpha.new
  }
  
  if (stop == "cut_average"){
    sumpar.pi <- rep(0,K)
    sumpar.prob <- lapply(1:K , function(z) lapply(1:d , function(j) rep(0,Nb_levels[j])))
    sumpar.alpha <- matrix(0,nrow=K,ncol=d)
  }
  
  if (mecha == "MNARz" | mecha == "MNARzj"){
    alpha.new <- Mechanism_SEM_GLM(YNA,YNA,Z.new,mecha)$alpha.new
  }
  
  
  Y.new <- SimulNA_SEM_Quali(YNA,Z.new,prob.new)
  if (typeof(Y.new)=="character"){
    return(Y.new)
  }
  
  res <- list(pik=pi.new, prob=prob.new, alpha=alpha.new)
  crit <- Critere_Quali(res,YNA,mecha)
  if (typeof(crit)=="character"){
    return(crit)
  }
  
  loglik_vec <- c(crit$loglik)
  
  if (stop == "loglikmax"){
    prec <- crit$loglik
  }
  
  r <- 0 
  
  while(r<rmax){
    
    #print(r)
    
    r <- r+1
    
    #################################
    #############SE STEP#############
    #################################
    
    
    Z.new <- SimulClass_SEM_Quali(YNA,Y.new,prob.new,alpha.new,pi.new)
    if (typeof(Z.new)=="character" ){
      return(Z.new)
    }
    
    Y.new <- SimulNA_SEM_Quali(YNA,Z.new,prob.new)
    if (typeof(Y.new)=="character" ){
      return(Y.new)
    }
    
    
    
    ################################
    #############M STEP#############
    ################################
    
    pi.new <- colSums(Z.new)/n
    
    prob.new <- lapply(1:K , function(k) lapply(1:d , function(j) sapply(1:Nb_levels[j] , function(l) sum(Y.new[Z.new[,k]==1,j]==l))/sum(Z.new[,k]==1)))
    
    
    if (mecha != "MCAR"){
      alpha.new <- Mechanism_SEM_GLM(YNA,Y.new,Z.new,mecha)$alpha.new
    }
    
    res <- list(pik=pi.new, prob=prob.new, alpha=alpha.new)
    crit <- Critere_Quali(res,YNA,mecha)
    if (typeof(crit)=="character"){
      return(crit)
    }
    loglik_vec <- c(loglik_vec,crit$loglik)
    
    
    if (stop == "loglikmax"){
      if(crit$loglik>prec){
        pi.new.save <- pi.new
        prob.new.save <- prob.new
        alpha.new.save <- alpha.new
        prec <- crit$loglik
      }
    }
    
    if (stop == "cut_average" & r>floor(0.8*rmax)){
      sumpar.pi <- sumpar.pi + pi.new
      sumpar.prob <- lapply(1:K , function(z) lapply(1:d , function(j) (sumpar.prob[[z]][[j]] + prob.new[[z]][[j]])))
      sumpar.alpha <- sumpar.alpha + alpha.new
    }
    
  }
  
  
  
  if (stop == "loglikmax"){
    pi.new <- pi.new.save
    prob.new <- prob.new.save
    alpha.new <- alpha.new.save
  }
  
  if (stop == "cut_average"){
    rcut <- floor(0.8*rmax)
    pi.new <- sumpar.pi/(rmax-rcut)
    prob.new <- lapply(1:K , function(z) lapply(1:d , function(j) sumpar.prob[[z]][[j]]/(rmax-rcut)))
    alpha.new <- sumpar.alpha/(rmax-rcut)
  }
  
  
  return(list(Ycomp=Y.new, Zcomp=Z.new, pik=pi.new, prob=prob.new, alpha=alpha.new, loglik_vec=loglik_vec))
  
  
}
  