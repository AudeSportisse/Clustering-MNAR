library(matrixcalc)
source('Internal_Functions_Gaussian.R')
source('Internal_Functions_All.R')
source('initialization.R')


######
# Name: EM_Gaussian
# Date: 2021/12/27
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR" or "MNARz" (character)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## rmax: maximum number of iterations (numeric, default is 100)
## init: choice of the initialization  (NULL, "sample" or list)
## stop: choice of stop criterium "cut_average" or "classical"
## tol: tolerence level to stop the algorithm (numeric)
## samplesize: size of the data sample used for the initialization, if init="sample" (numeric)
#####



EM_Gaussian <- function(YNA,K,mecha,diag,rmax,init,stop="classical",tol=0.0001,samplesize=NULL){
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)), ncol = ncol(YNA))
  
  ########################################
  #############INITIALIZATION#############
  ########################################
  
  init <- init_EM_Gaussian(YNA, K, mecha, diag, init, samplesize)
  pi.new <- init$pi.init
  mu.new <- init$mu.init
  sigma.new <- init$sigma.init
  alpha.new <- init$alpha.init
  
  
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
  
  loglik_res <- Loglikelihood_obs_Gaussian(YNA,mu.new,sigma.new,alpha.new,pi.new)
  if (typeof(loglik_res)=="character"){
    return(loglik_res)
  }
  loglik_vec <- c(loglik_res$loglik_obs)
  
  prec <- -Inf
  
  r <- 0
    
  while((loglik_res$loglik_obs-prec)>tol | r<rmax){
    
    #print(r)
    
    r <- r+1
    
    #################################
    #############E STEP##############
    #################################
    
    sigma.tilde <- list()
    for (i in 1:n){
      sigma.tilde[[i]] <- list()
      for (k in 1:K){
        sigma.tilde[[i]][[k]] <- matrix(0,nrow=d,ncol=d)
      }
    }
    for (k in 1:K){
      for (i in 1:n){
        indobs <- which(C[i,]==0)
        indmiss <- which(C[i,]==1)
        if (length(indmiss)>0){
          if(length(indobs>1)){
            #print(det(data.matrix(sigma.new[[k]][indobs, indobs])))
            if(abs(det(data.matrix(sigma.new[[k]][indobs, indobs])))<=10^(-15)){return("error: determinant of a sub-matrix of the covariance matrix is null")}
          }else{
            if(abs(sigma.new[[k]][indobs, indobs])<=10^-15){return("error: determinant of a sub-matrix of the covariance matrix is null")}
          }
          mu.tilde <- mu.new[[k]][indmiss] + sigma.new[[k]][indmiss,indobs]%*%solve(sigma.new[[k]][indobs,indobs])%*%(YNA[i,indobs]-mu.new[[k]][indobs])
          sigma.tilde.miss <- sigma.new[[k]][indmiss,indmiss]-sigma.new[[k]][indmiss,indobs]%*%solve(sigma.new[[k]][indobs,indobs])%*%sigma.new[[k]][indobs,indmiss]
          sigma.tilde[[i]][[k]][indmiss,indmiss] <- sigma.tilde.miss
        }
      }
    }
    
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
      mu.old = mu.new[[k]]
      for (i in 1:n){
        impute <- YNA[i,]
        impute[which(C[i,]==1)] <- mu.old[which(C[i,]==1)]
        mu.new.num  <- mu.new.num + (impute * tik[i,k])
      }
      mu.new[[k]] <- mu.new.num / sum(tik[,k])
      for (i in 1:n){
        impute <- YNA[i,]
        impute[which(C[i,]==1)] <- mu.old[which(C[i,]==1)]
        if (diag==TRUE){
          sigma.new.num <- sigma.new.num + (diag(diag((impute - mu.new[[k]]) %*% t(impute - mu.new[[k]]))) + sigma.tilde[[i]][[k]] ) * tik[i,k]  
        }else{
          sigma.new.num <- sigma.new.num + ((impute - mu.new[[k]]) %*% t(impute - mu.new[[k]]) + sigma.tilde[[i]][[k]] ) * tik[i,k]  
        }
      }
      sigma.new[[k]] <- sigma.new.num / sum(tik[,k]) 
      }#else{
        #for (j in 1:d){
         # impute <- YNA[,j]
         # impute[which(C[,j]==1)] <- mu.new[[k]][j]
         # mu.new[[k]][j] <-  sum(impute * tik[,k]) / sum(tik[,k])
         # sigma.new[[k]][j,j] <- sum((impute - mu.new[[k]][j])*(impute - mu.new[[k]][j])*tik[,k]) / sum(tik[,k]) 
        #}
      #}
      #if (!is.positive.definite(sigma.new[[k]])){
       # return("error: the estimated covariance matrix is not positive definite.")
      #}
      
    
    if (mecha != "MCAR"){
      alpha.new <- Mechanism_EM_GLM(YNA,tik,mecha)
    }
    
    prec <- loglik_res$loglik_obs
    loglik_res <- Loglikelihood_obs_Gaussian(YNA,mu.new,sigma.new,alpha.new,pi.new)
    if (typeof(loglik_res)=="character"){
      return(loglik_res)
    }
    loglik_vec <- c(loglik_vec,loglik_res$loglik_obs)
    
    if (stop == "cut_average" & r>floor(0.8*rmax)){
      sumpar.pi <- sumpar.pi + pi.new
      sumpar.mu <- lapply(1:K , function(z) sumpar.mu[[z]]+mu.new[[z]])
      sumpar.sigma <- lapply(1:K , function(z) sumpar.sigma[[z]]+sigma.new[[z]])
      sumpar.alpha <- sumpar.alpha + alpha.new
    }
    
    if(is.nan(loglik_res$loglik_obs-prec)){
      return("error: observed loglikelihood is not computable")
    }
    
  }
  
  
  if (stop == "cut_average"){
    rcut <- floor(0.8*rmax)
    pi.new <- sumpar.pi/(r-rcut)
    mu.new <- lapply(1:K , function(z) sumpar.mu[[z]]/(r-rcut))
    sigma.new <- lapply(1:K , function(z) sumpar.sigma[[z]]/(r-rcut))
    alpha.new <- sumpar.alpha/(r-rcut)
  }
  
  error <- loglik_res$loglik_obs<prec
  
  
  return(list(pik=pi.new, mu=mu.new, sigma=sigma.new, alpha=alpha.new,beta=matrix(0,nrow=K,ncol=d), loglik_vec=loglik_vec, prec=prec, error=error))
  
}

