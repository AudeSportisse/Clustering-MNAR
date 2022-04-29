
######
# Name: SimulNA_SEM_Gaussian
# Date: 2021/12/21
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## Z: current partition (data matrix of size n,K)
## mu: current value of the means (list of length K, with vectors of size d)
## sigma: current value of the covariance matrices (list of length K, with matrices of size d,d)
## alpha: current value of mechanism coefficients for the effect of the class membership (data matrix of size K,d)
## beta: current value of mechanism coefficients for the effect of the variables themselves (data matrix of size K,d)
## L: current value of the latent variable (for "MNARy" and "MNARyz") (data matrix of size n,d)
#####



SimulNA_SEM_Gaussian <- function(YNA,mecha,diag,Z,mu,sigma,alpha,beta,L){
  
  n <- dim(YNA)[1]
  K <- dim(Z)[2]
  
  Yres <- YNA 
  for (i in 1:n){
    var_na = which(is.na(Yres[i,]))
    if (length(var_na)>0){
      var_obs <- (1:dim(Yres)[2])[-var_na]
      k <- which(Z[i,]==1) 
      mu_tilde <- mu[[k]][var_na] + (sigma[[k]][var_na,var_obs]%*%solve(sigma[[k]][var_obs,var_obs]))%*%(Yres[i,var_obs]-mu[[k]][var_obs])
      
      if (length(var_obs)>1){
        sigma_tilde <- sigma[[k]][var_na,var_na] - (sigma[[k]][var_na,var_obs]%*%solve(sigma[[k]][var_obs,var_obs]))%*%sigma[[k]][var_obs,var_na] 
      }else{
        sigma_tilde <- sigma[[k]][var_na,var_na] - (sigma[[k]][var_na,var_obs])%*%t(sigma[[k]][var_na,var_obs])*(1/sigma[[k]][var_obs,var_obs])
      }
      
      if (mecha == "MNARy" | mecha == "MNARyz" | mecha == "MNARyk" | mecha == "MNARykz" | mecha == "MNARyzj" | mecha == "MNARykzj"){
        sigma_tilde_inv <- solve(sigma_tilde)
        if (length(var_na)==1){
          beta_diag <- beta[k,var_na]
        }else{
          beta_diag <- diag(beta[k,var_na]) 
        }
        Cov_Y <- solve(sigma_tilde_inv + beta_diag%*%beta_diag)
        Mean_Y <- Cov_Y %*% ( sigma_tilde_inv%*%mu_tilde + beta_diag%*%L[i,var_na] - beta_diag%*%alpha[k,var_na])
        if (length(var_na)>1){
          #if (sum(is.na(Cov_Y))>0){return("error: covariance matrix for sampling the missing values contains NA")}
          if(!is.positive.definite(Cov_Y)){return("error: covariance matrix for sampling the missing values is not definite positive")}
          Yres[i,var_na] <- mvrnorm(1,c(Mean_Y),Cov_Y)  
        }else{
          if(Cov_Y<=0){return("error: covariance matrix for sampling the missing values is not definite positive")}
          Yres[i,var_na] <- rnorm(1,Mean_Y,sqrt(Cov_Y)) 
        }
        
      }else{
        if (length(var_na)>1){
          if(!is.symmetric.matrix(sigma_tilde)){
            sigma_tilde.diag <- diag(sigma_tilde)
            sigma_tilde[lower.tri(sigma_tilde,diag=T)]=0
            sigma_tilde <- sigma_tilde+ t(sigma_tilde ) + diag(sigma_tilde.diag)
          }
          if(!is.positive.definite(sigma_tilde)){return("error: covariance matrix for sampling the missing values is not definite positive")}
          Yres[i,var_na] <- mvrnorm(1,mu_tilde,sigma_tilde)
        }else{
          if(sigma_tilde<=0){return("error: covariance matrix for sampling the missing values is not definite positive")}
          Yres[i,var_na] <- rnorm(1,mu_tilde,sqrt(sigma_tilde))
        } 
      }
    }
  }
  
  return(Yres)
}




######
# Name: SimulL_SEM_Gaussian
# Date: 2021/12/21
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## Y: current imputed matrix (data matrix of size n,d)
## Z: current partition (data matrix of size n,K)
## alpha: current value of mechanism coefficients for the effect of the class membership (data matrix of size K,d)
## beta: current value of mechanism coefficients for the effect of the variables themselves (data matrix of size K,d)
#####


SimulL_SEM_Gaussian <- function(YNA,Y,Z,alpha,beta){
  
  n = dim(YNA)[1]
  d = dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)),ncol=ncol(YNA))
  
  L = matrix(0, nrow=n, ncol = ncol(YNA))
  
  for (i in 1:n){
    k = which(Z[i,]==1)
    upper_trunc = numeric(d)
    lower_trunc = numeric(d)
    upper_trunc[C[i,]==0] <- 0
    upper_trunc[C[i,]==1] <- Inf
    lower_trunc[C[i,]==0] <- -Inf
    lower_trunc[C[i,]==1] <- 0
    L[i,] <- rtmvnorm(1, mean=alpha[k,]+beta[k,]*Y[i,],lower=lower_trunc,upper=upper_trunc,algorithm = "gibbs")
  }
  
  return(L)
}



######
# Name: SimulClass_SEM_Gaussian
# Date: 2021/12/21
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## Y: current imputed matrix (data matrix of size n,d)
## L: current value of the latent variable (for "MNARy" and "MNARyz") (data matrix of size n,d)
## mu: current value of the means (list of length K, with vectors of size d)
## sigma: current value of the covariance matrices (list of length K, with matrices of size d,d)
## alpha: current value of mechanism coefficients for the effect of the class membership (data matrix of size K,d)
## beta: current value of mechanism coefficients for the effect of the variables themselves (data matrix of size K,d)
## prop.pi: current proportion per class (vector of length K)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
#####


SimulClass_SEM_Gaussian <- function(YNA,Y,L,mu,sigma,alpha,beta,prop.pi,mecha){
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)),ncol=ncol(YNA))
  K <- length(prop.pi)
  
  Zres <- matrix(0,nrow=n,ncol=K)
  tau <- matrix(0,nrow=n,ncol=K)
  upper_trunc = matrix(NA,nrow=n,ncol=d)
  lower_trunc = matrix(NA,nrow=n,ncol=d)
  upper_trunc[C==0] <- 0
  upper_trunc[C==1] <- Inf
  lower_trunc[C==0] <- -Inf
  lower_trunc[C==1] <- 0

  for (k in 1:K){
    int_pn <- matrix(alpha[k,],nrow=n,ncol=d,byrow=TRUE)+Y%*%diag(beta[k,])
    pn <- pnorm(int_pn)
    logMecha <- rowSums(C*log(pn)+(1-C)*log(1-pn))
    if ((mecha == "MNARy" | mecha == "MNARyz" | mecha == "MNARyk" | mecha == "MNARykz" | mecha == "MNARyzj" | mecha == "MNARykzj")){
      logL <- sapply(1:n,function(i) dtmvnorm(L[i,],mean=int_pn[i,],lower=lower_trunc[i,],upper=upper_trunc[i,],log=TRUE)) 
    }else{
      logL <- 0
    }
    tau[,k] = logMecha+dmvnorm(Y,mu[[k]],sigma[[k]],log=TRUE)+log(prop.pi[k]) + logL
  }
  
  if (length(which(is.infinite(tau)))>0 | length(which(is.na(tau)))>0){
    print(length(which(is.infinite(tau))))
    print(length(which(is.na(tau))))
    if ((length(which(is.infinite(tau)))+length(which(is.na(tau))))>=n*K){
      return("Degenerescence loglikelihood, computation of tau.")
    }
    val_remp = min(tau[-which(is.infinite(tau)|is.na(tau))])-1000
    tau[which(is.infinite(tau))]=val_remp
    tau[which(is.na(tau))]=val_remp
  }
  
  IndMax <- apply(tau,1,which.max)
  LogFmax <- sapply(1:length(IndMax) , function(i) tau[i,IndMax[i]])
  LogSum <- LogFmax + log(rowSums(exp(tau-LogFmax)))
  ProbaClust <- exp(tau-LogSum)
  Zres <- t(sapply(1:nrow(ProbaClust) , function(i) {
    return(rmultinom(1,1,ProbaClust[i,]))
  }))
  if (K==1){
    Zres=t(Zres)
  }
  return(Zres)
}



######
# Name: Critere_Gaussian
# Date: 2021/12/21
# Description:
# Arguments:
## ResAlgo: object of the type returning by the function SEM_Gaussian (list)
## YNA: matrix containing missing values (data matrix of size n,d)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
## Partition_true: true partition (vector of size n)
#####


Critere_Gaussian <- function(ResAlgo,YNA,mecha,Partition_true=NULL){
  
  C <- is.na(YNA)
  
  Pattern <- matrix(C[!duplicated(C),],nrow=sum(!duplicated(C)))
  
  Y <- YNA 
  n <- nrow(Y)
  d <- ncol(Y)
  mu <- ResAlgo$mu
  sigma <- ResAlgo$sigma
  pik <- ResAlgo$pik
  K <- length(pik)
  alpha <- ResAlgo$alpha
  beta <- ResAlgo$beta
  
  
  valcrit <- 0
  Zall <- matrix(NA,ncol=K,nrow=n)
  proba_class <- matrix(NA,ncol=K,nrow=n)
  loglik=0
  
  for (i in 1:nrow(Pattern)){
    wherePat <- which(sapply(1:nrow(C),function(x) prod(C[x,]==Pattern[i,])==1))
    Y_Pat <- rbind(Y[wherePat,])
    C_Pat <- rbind(C[wherePat,])
    
    var_obs <- which(Pattern[i,]==0)
    Y_obs <- as.matrix(Y_Pat[,var_obs])
    mu_obs <- lapply(mu, function(x) x[var_obs])
    sigma_obs <- lapply(sigma , function(x) matrix(x[var_obs,var_obs],ncol=length(var_obs),nrow=length(var_obs)))
    d_varobs <- length(var_obs)
    C_obs <- matrix(C_Pat[,var_obs],ncol=length(var_obs))
    
    if (nrow(Y_Pat)==1){
      Y_obs <- matrix(Y_Pat[,var_obs],nrow=nrow(Y_Pat),ncol=length(var_obs))
      C_obs <- matrix(C_Pat[,var_obs],nrow=nrow(Y_Pat),ncol=length(var_obs))
    }
    
    
    LogF <- matrix(NA,nrow=length(wherePat),ncol=K)
    for (ik in 1:K){
      term1_log = log(pik[ik]) + dmvnorm(Y_obs,mean=mu_obs[[ik]],sigma=sigma_obs[[ik]],log=TRUE)
      term2_log = 0
      for (j in 1:d){
        if (Pattern[i,j]==FALSE){
          term_mecha = pnorm(alpha[ik,j]+beta[ik,j]*Y_Pat[,j])
          if(sum(sapply(term_mecha, function(z) z>= 1))==length(term_mecha)){
            return("Degenerescence (mechanism): loglikelihood is not computable")
          }
          res_mecha = log(1-term_mecha)
          term2_log = term2_log + (1-C_Pat[,j]) * res_mecha
        }else{
          if (mecha == "MCAR" | mecha == "MNARz" | mecha == "MNARzj"){
            res_miss <- log(pnorm(alpha[ik,j]))
          }else{
            res_miss <- c()
            for (ipat in 1:nrow(Y_Pat)){
              mu_tilde <- mu[[ik]][j] + sigma[[ik]][j,var_obs]%*%solve(sigma[[ik]][var_obs,var_obs])%*%(Y_obs[ipat,]-mu[[ik]][var_obs])
              sigma_tilde <- sigma[[ik]][j,j] - sigma[[ik]][j,var_obs]%*%solve(sigma[[ik]][var_obs,var_obs])%*%sigma[[ik]][j,var_obs]
              integrand <- function(u) dnorm(u, mu_tilde, sqrt(sigma_tilde)) * pnorm(alpha[ik,j]+beta[ik,j]*u)
              res <- integrate(integrand, lower = -Inf, upper = Inf, stop.on.error = FALSE)
              if (typeof(res)=="character" | res$value==0){
                return("Degenerescence (mechanism): loglikelihood (numerator) is not computable")
              }else{
                res_miss = c(res_miss,log(res$value)) 
              }
            }
          }
          term2_log = term2_log + res_miss
        }
      }
      LogF[,ik] = term1_log + term2_log
    }
    
    tauik_help <- apply(LogF, 1, max)
    tauik <- exp(sweep(LogF, 1, tauik_help, "-"))
    
    logcal <- sum(tauik_help) + sum(log(rowSums(tauik)))
    tauik <- tauik / rowSums(tauik)
    proba_class[wherePat,] <- tauik
    
    maxT <- apply(tauik,1,which.max)
    Zres <- matrix(0,ncol=K,nrow=length(maxT))
    for (m in 1:length(maxT)){
      Zres[m,maxT[m]] <- 1
    }
    
    loglik <- loglik + logcal
    valcrit <- valcrit + sum(sapply(1:nrow(Zres) , function(l) LogF[l,Zres[l,]==1])) 
    Zall[wherePat,] <- Zres
  }
  
  df <- K-1 + K*d + K*d ###K*d*(d+1)/2 A CHANGER SI DIAG!=TRUE
  if (mecha == "MCAR"){
    df <- df + 1
  }else if(mecha == "MNARz"){
    df <- df + K
  }else if(mecha == "MNARy"){
    df <- df + d + 1
  }else if(mecha == "MNARyz"){
    df <- df + K + d
  }else if(mecha == "MNARyk"){
    df <- df + K*d+1 
  }else if(mecha == "MNARzj"){
    df <-  df + K*d
  }else if(mecha == "MNARykz"){
    df <- df + K*(d+1)
  }else if(mecha == "MNARyzj"){
    df <- df + (K+1)*d
  }else if(mecha == "MNARykzj"){
    df <- df + 2*K*d
  }
  
  ICL <- valcrit - df/2 * log(n)
  Partition_estim <- apply(Zall, 1, function(z) which(z==1))
  
  if (is.null(Partition_true)){
    ARI_fin <- NULL
  }else{
    ARI_fin <- ARI(Partition_estim,Partition_true) 
  }
  return(list(ICL=ICL,Z=Zall,ARI=ARI_fin,loglik=loglik,proba_class=proba_class))
}


######
# Name: Loglikelihood_obs_Gaussian
# Date: 2021/12/21
# Description: used for the EM algorithm ("MCAR" or "MNARz")
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## mu: current value of the means (list of length K, with vectors of size d)
## sigma: current value of the covariance matrices (list of length K, with matrices of size d,d)
## alpha: current value of mechanism coefficients for the effect of the class membership (data matrix of size K,d)
## prop.pi: current proportion per class (vector of length K)
#####


Loglikelihood_obs_Gaussian <- function(YNA,mu,sigma,alpha,prop.pi){
  
  n = dim(YNA)[1]
  d = dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)),ncol=ncol(YNA))
  K <- length(prop.pi)
  
  logprobcond <- matrix(0, n, K)
  
  Pattern <- matrix(C[!duplicated(C),],nrow=sum(!duplicated(C)))
  
  for (i in 1:nrow(Pattern)){
    wherePat <- which(sapply(1:nrow(C),function(x) prod(C[x,]==Pattern[i,])==1))
    Y_Pat <- rbind(YNA[wherePat,])
    C_Pat <- rbind(C[wherePat,])
    
    var_obs <- which(Pattern[i,]==0)
    Y_obs <- as.matrix(Y_Pat[,var_obs])
    mu_obs <- lapply(mu, function(x) x[var_obs])
    sigma_obs <- lapply(sigma , function(x) matrix(x[var_obs,var_obs],ncol=length(var_obs),nrow=length(var_obs)))
    d_varobs <- length(var_obs)
    C_obs <- matrix(C_Pat[,var_obs],ncol=length(var_obs))
    
    if (nrow(Y_Pat)==1){
      Y_obs <- matrix(Y_Pat[,var_obs],nrow=nrow(Y_Pat),ncol=length(var_obs))
      C_obs <- matrix(C_Pat[,var_obs],nrow=nrow(Y_Pat),ncol=length(var_obs))
    }
    
    
    for (ik in 1:K){
      term1_log = log(prop.pi[ik]) + dmvnorm(Y_obs,mean=mu_obs[[ik]],sigma=sigma_obs[[ik]],log=TRUE)
      term2_log = 0
      for (j in 1:d){
        if (Pattern[i,j]==FALSE){
          term_mecha = pnorm(alpha[ik,j])
          if(sum(sapply(term_mecha, function(z) z>= 1))==length(term_mecha)){
            return("Degenerescence (mechanism): loglikelihood is not computable")
          }
          res_mecha = log(1-term_mecha)
          term2_log = term2_log + (1-C_Pat[,j]) * res_mecha
        }else{
          term2_log = term2_log + log(pnorm(alpha[ik,j]))
        }
      }
      logprobcond[wherePat,ik] = term1_log + term2_log
    }
  }
    
  
  normval <- apply(logprobcond, 1, max)
  logprobcond <- exp(sweep(logprobcond, 1, normval, "-"))
  loglik_obs <- sum(normval) + sum(log(rowSums(logprobcond))) 
  tik <- logprobcond / rowSums(logprobcond)
  
  return(list(loglik_obs=loglik_obs,tik=tik))
  
}


######
# Name: Critere_Gaussian_withoutNA
# Date: 2021/12/21
# Description:
# Arguments:
## ResAlgo: object of the type returning by the function SEM_Gaussian (list)
## Y: complete matrix (data matrix of size n,d)
## Partition_true: true partition (vector of size n)
#####


Critere_Gaussian_withoutNA <- function(ResAlgo,Y,Partition_true=NULL){
  
  n <- nrow(Y)
  d <- ncol(Y)
  mu <- ResAlgo$mu
  sigma <- ResAlgo$sigma
  pik <- ResAlgo$pik
  K <- length(pik)

    
  LogF <- matrix(NA,nrow=n,ncol=K)
  for (ik in 1:K){
    LogF[,ik] = log(pik[ik]) + dmvnorm(Y,mean=mu[[ik]],sigma=sigma[[ik]],log=TRUE)
  }  
    
  tauik_help <- apply(LogF, 1, max)
  tauik <- exp(sweep(LogF, 1, tauik_help, "-"))
    
  logcal <- sum(tauik_help) + sum(log(rowSums(tauik)))
  tauik <- tauik / rowSums(tauik)
  proba_class <- tauik
    
  maxT <- apply(tauik,1,which.max)
  Zres <- matrix(0,ncol=K,nrow=length(maxT))
  for (m in 1:length(maxT)){
    Zres[m,maxT[m]] <- 1
  }
    
  loglik <- logcal
  valcrit <- sum(sapply(1:nrow(Zres) , function(l) LogF[l,Zres[l,]==1])) 
  Zall <- Zres
  
  
  df <- K-1 + K*d + K*d*(d+1)/2
  
  ICL <- valcrit - df/2 * log(n)
  Partition_estim <- apply(Zall, 1, function(z) which(z==1))
  
  if (is.null(Partition_true)){
    ARI_fin <- NULL
  }else{
    ARI_fin <- ARI(Partition_estim,Partition_true) 
  }
  return(list(ICL=ICL,Z=Zall,ARI=ARI_fin,loglik=loglik,proba_class=proba_class))
}


######
# Name: Loglikelihood_Gaussian
# Date: 2021/12/21
# Description: used for the EM algorithm ("MCAR" or "MNARz")
# Arguments:
## Y: complete matrix (data matrix of size n,d)
## mu: current value of the means (list of length K, with vectors of size d)
## sigma: current value of the covariance matrices (list of length K, with matrices of size d,d)
## prop.pi: current proportion per class (vector of length K)
#####


Loglikelihood_Gaussian <- function(Y,mu,sigma,prop.pi){
  
  n = dim(Y)[1]
  d = dim(Y)[2]
  K <- length(prop.pi)
  
  logprobcond <- matrix(0, n, K)
  for (ik in 1:K){
    logprobcond[,ik] = log(prop.pi[ik]) + dmvnorm(Y,mean=mu[[ik]],sigma=sigma[[ik]],log=TRUE)
  }
  

  normval <- apply(logprobcond, 1, max)
  logprobcond <- exp(sweep(logprobcond, 1, normval, "-"))
  loglik <- sum(normval) + sum(log(rowSums(logprobcond))) 
  tik <- logprobcond / rowSums(logprobcond)
  
  return(list(loglik=loglik,tik=tik))
  
}

