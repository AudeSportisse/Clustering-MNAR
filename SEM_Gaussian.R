library(matrixcalc)
library(truncnorm)
library(tmvtnorm)
source('Internal_Functions_Gaussian.R')
source('Internal_Functions_All.R')
source('initialization.R')

######
# Name: SEM_Gaussian
# Date: 2021/12/21
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## rmax: maximum number of iterations (numeric, default is 100)
## init: choice of the initialization  (NULL, "sample" or list)
## stop: choice of stop criterium "cut_average", "loglikmax" or "rmax"
## samplesize: size of the data sample used for the initialization, if init="sample" (numeric)
#####

SEM_Gaussian <- function(YNA, K, mecha, diag, rmax = 100, init = NULL, stop = "rmax", samplesize = NULL) {
    
    n <- dim(YNA)[1]
    d <- dim(YNA)[2]
    C <- matrix(as.numeric(is.na(YNA)), ncol = ncol(YNA))
    
    
    ########################################
    #############INITIALIZATION#############
    ########################################
    
    init <- init_SEM_Gaussian(YNA, K, mecha, diag, init, samplesize)
    Z.new <- init$Z.init
    pi.new <- init$pi.init
    mu.new <- init$mu.init
    sigma.new <- init$sigma.init
    alpha.new <- init$alpha.init
    beta.new <- init$beta.init
    L.new <- init$L.init
    
    if (stop == "loglikmax"){
      pi.new.save <- pi.new
      mu.new.save <- mu.new
      sigma.new.save <- sigma.new
      alpha.new.save <- alpha.new
      beta.new.save <- beta.new 
    }
    
    if (stop == "cut_average"){
      sumpar.pi <- rep(0,K)
      sumpar.mu <- list()
      sumpar.sigma <- list()
      for (k in 1:K){
        sumpar.mu[[k]] <- rep(0,d)
        sumpar.sigma[[k]] <- matrix(0,nrow=d,ncol=d)
      }
      sumpar.alpha <- matrix(0,nrow=K,ncol=d)
      sumpar.beta <- matrix(0,nrow=K,ncol=d)
    }
    
    if(sum(unlist(lapply(1:K,function(z) is.positive.definite(sigma.new[[z]]))))<K){
      return("error: an estimator of a covariance matrix is not positive definite")
    }
    
    if (mecha == "MNARz" | mecha == "MNARzj"){
      param.mecha <- Mechanism_SEM_GLM(YNA,YNA,Z.new,mecha) #Mechanism_SEM_GLM(YNA,Y.new,Z.new,mecha) mais Y.new non défini
      alpha.new <- param.mecha$alpha.new
      beta.new <- param.mecha$beta.new
    }

    
    Y.new <- SimulNA_SEM_Gaussian(YNA,mecha,diag,Z.new,mu.new,sigma.new,alpha.new,beta.new,L.new)
    if (typeof(Y.new)=="character"){
      return(Y.new)
    }
    
    if (mecha == "MNARy" | mecha == "MNARyz" | mecha == "MNARyk" | mecha == "MNARykz" | mecha == "MNARyzj" | mecha == "MNARykzj"){
      param.mecha <- Mechanism_SEM_GLM(YNA,Y.new,Z.new,mecha)
      alpha.new <- param.mecha$alpha.new
      beta.new <- param.mecha$beta.new
    }
    
    res <- list(pik=pi.new, mu=mu.new, sigma=sigma.new, alpha=alpha.new, beta=beta.new)
    crit <- Critere_Gaussian(res,YNA,mecha)
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
      
      if (mecha == "MNARy" | mecha == "MNARyz" | mecha == "MNARyk" | mecha == "MNARykz" | mecha == "MNARyzj" | mecha == "MNARykzj"){
        L.new <- SimulL_SEM_Gaussian(YNA,Y.new,Z.new,alpha.new,beta.new)
      }
      
      Z.new <- SimulClass_SEM_Gaussian(YNA,Y.new,L.new,mu.new,sigma.new,alpha.new,beta.new,pi.new,mecha)
      if (typeof(Z.new)=="character" ){
        return(Z.new)
      }
      
      Y.new <- SimulNA_SEM_Gaussian(YNA,mecha,diag,Z.new,mu.new,sigma.new,alpha.new,beta.new,L.new)
      if (typeof(Y.new)=="character" ){
        return(Y.new)
      }

      
      
      ################################
      #############M STEP#############
      ################################
      
      pi.new <- colSums(Z.new)/n

      for (k in 1:K){
        Obs_classk <- which(Z.new[,k]==1)
        if (length(Obs_classk)<2){
          return("error: no enough observations for one class")
        }else{
          mu.new[[k]] <- colMeans(Y.new[Obs_classk,])

          if (diag==TRUE){
            sigma.new[[k]] <- diag(apply(Y.new[Obs_classk,],2,var))
          }else{
            sigma.new[[k]] <- cov(Y.new[Obs_classk,])
          }
          if (sum(is.na(sigma.new[[k]]))>0){return("error: an estimator of a covariance matrix contains NA")}
          if (!is.positive.definite(sigma.new[[k]])){
              return("error: an estimator of a covariance matrix is not positive definite")
          }
        }
      }

      if (mecha != "MCAR"){
       param.mecha <- Mechanism_SEM_GLM(YNA,Y.new,Z.new,mecha)
       alpha.new <- param.mecha$alpha.new
       beta.new <- param.mecha$beta.new
      }
      
      res <- list(pik=pi.new, mu=mu.new, sigma=sigma.new, alpha=alpha.new, beta=beta.new)
      crit <- Critere_Gaussian(res,YNA,mecha)
      if (typeof(crit)=="character"){
        return(crit)
      }
      loglik_vec <- c(loglik_vec,crit$loglik)
      
      
      if (stop == "loglikmax"){
        if(crit$loglik>prec){
          pi.new.save <- pi.new
          mu.new.save <- mu.new
          sigma.new.save <- sigma.new
          alpha.new.save <- alpha.new
          beta.new.save <- beta.new 
          prec <- crit$loglik
        }
      }
      
      if (stop == "cut_average" & r>floor(0.8*rmax)){
        sumpar.pi <- sumpar.pi + pi.new
        sumpar.mu <- lapply(1:K , function(z) sumpar.mu[[z]]+mu.new[[z]])
        sumpar.sigma <- lapply(1:K , function(z) sumpar.sigma[[z]]+sigma.new[[z]])
        sumpar.alpha <- sumpar.alpha + alpha.new
        sumpar.beta <- sumpar.beta + beta.new
      }
      
    }
    
    
    
    if (stop == "loglikmax"){
      pi.new <- pi.new.save
      mu.new <- mu.new.save
      sigma.new <- sigma.new.save
      alpha.new <- alpha.new.save
      beta.new <- beta.new.save
    }
    
    if (stop == "cut_average"){
      rcut <- floor(0.8*rmax)
      pi.new <- sumpar.pi/(rmax-rcut)
      mu.new <- lapply(1:K , function(z) sumpar.mu[[z]]/(rmax-rcut))
      sigma.new <- lapply(1:K , function(z) sumpar.sigma[[z]]/(rmax-rcut))
      alpha.new <- sumpar.alpha/(rmax-rcut)
      beta.new <- sumpar.beta/(rmax-rcut)
    }
    
    
    return(list(Ycomp=Y.new, Zcomp=Z.new, pik=pi.new, mu=mu.new, sigma=sigma.new, alpha=alpha.new, beta=beta.new, loglik_vec=loglik_vec))
    
    
  }



