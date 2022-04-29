library(matrixcalc)
source('Internal_Functions_Gaussian.R')
source('Internal_Functions_All.R')
source('initialization.R')


######
# Name: EM_Gaussian_withoutNA
# Date: 2021/12/27
# Description:
# Arguments:
## Y: complete matrix (data matrix of size n,d)
## K: number of classes to use (numeric)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## rmax: maximum number of iterations (numeric, default is 100)
## init: choice of the initialization  (NULL, "sample" or list)
## stop: choice of stop criterium "cut_average" or "classical"
## tol: tolerence level to stop the algorithm (numeric)
## samplesize: size of the data sample used for the initialization, if init="sample" (numeric)
#####



EM_Gaussian_withoutNA <- function(Y,K,diag,rmax,init,stop="classical",tol=0.0001,samplesize=NULL){
  
  n <- dim(Y)[1]
  d <- dim(Y)[2]

  ########################################
  #############INITIALIZATION#############
  ########################################
  
  init <- init_EM_Gaussian_withoutNA(Y, K, diag, init, samplesize)
  pi.new <- init$pi.init
  mu.new <- init$mu.init
  sigma.new <- init$sigma.init
  
  
  if (stop == "cut_average"){
    sumpar.pi <- rep(0,K)
    sumpar.mu <- list()
    sumpar.sigma <- list()
    for (k in 1:K){
      sumpar.mu[[k]] <- rep(0,d)
      sumpar.sigma[[k]] <- matrix(0,nrow=d,ncol=d)
    }
  }
  
  if(sum(unlist(lapply(1:K,function(z) is.positive.definite(sigma.new[[z]]))))<K){
    return("error: an estimator of a covariance matrix is not positive definite")
  }
  
  loglik_res <- Loglikelihood_Gaussian(Y,mu.new,sigma.new,pi.new)
  if (typeof(loglik_res)=="character"){
    return(loglik_res)
  }
  loglik_vec <- c(loglik_res$loglik)
  
  prec <- -Inf
  
  r <- 0
  
  while((loglik_res$loglik-prec)>tol | r<rmax){
    
    #print(r)
    
    r <- r+1
    
    #################################
    #############E STEP##############
    #################################
    

    tik <- loglik_res$tik
    
    #################################
    #############M STEP##############
    #################################
    
    pi.new <- colSums(tik) / n
    if (sum(colSums(tik)<=1)>=1){
      return("error: no enough observations for one class")
    }
    
    
    
    for (k in 1:K){
      mu.new.num <- rep(0,d)
      sigma.new.num <- matrix(0,nrow=d,ncol=d)
      for (i in 1:n){
        mu.new.num  <- mu.new.num + (Y[i,] * tik[i,k])
      }
      mu.new[[k]] <- mu.new.num / sum(tik[,k])
      for (i in 1:n){
        if (diag==TRUE){
          sigma.new.num <- sigma.new.num + (diag(diag((Y[i,] - mu.new[[k]]) %*% t(Y[i,] - mu.new[[k]]))) ) * tik[i,k]  
        }else{
          sigma.new.num <- sigma.new.num + ((Y[i,] - mu.new[[k]]) %*% t(Y[i,] - mu.new[[k]])) * tik[i,k]  
        }
      }
      sigma.new[[k]] <- sigma.new.num / sum(tik[,k]) 
    }
    
    
    prec <- loglik_res$loglik
    loglik_res <- Loglikelihood_Gaussian(Y,mu.new,sigma.new,pi.new)
    if (typeof(loglik_res)=="character"){
      return(loglik_res)
    }
    loglik_vec <- c(loglik_vec,loglik_res$loglik)
    
    if (stop == "cut_average" & r>floor(0.8*rmax)){
      sumpar.pi <- sumpar.pi + pi.new
      sumpar.mu <- lapply(1:K , function(z) sumpar.mu[[z]]+mu.new[[z]])
      sumpar.sigma <- lapply(1:K , function(z) sumpar.sigma[[z]]+sigma.new[[z]])
    }
    
    if(is.nan(loglik_res$loglik-prec)){
      return("error: observed loglikelihood is not computable")
    }
    
  }
  
  
  if (stop == "cut_average"){
    rcut <- floor(0.8*rmax)
    pi.new <- sumpar.pi/(r-rcut)
    mu.new <- lapply(1:K , function(z) sumpar.mu[[z]]/(r-rcut))
    sigma.new <- lapply(1:K , function(z) sumpar.sigma[[z]]/(r-rcut))
  }
  
  error <- loglik_res$loglik<prec
  
  
  return(list(pik=pi.new, mu=mu.new, sigma=sigma.new, loglik_vec=loglik_vec, prec=prec, error=error))
  
}

