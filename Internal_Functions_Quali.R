######
# Name: SimulNA_SEM_Quali
# Date: 2022/02/17
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## Z: current partition (data matrix of size n,K)
## prob: current value of the probabilities (list of length K, with list of length d of vectors of size the number of levels)
#####

SimulNA_SEM_Quali <- function(YNA,Z,prob){ 
  
  n <- dim(YNA)[1]
  K <- dim(Z)[2]
  
  Yres <- YNA
  for (k in 1:K){
    obs_k <- which(Z[,k]==1)
    Y_obs_k <- Yres[obs_k,]
    for (j in 1:ncol(Y_obs_k)){
      obs_na <- which(is.na(Y_obs_k[,j]))
      if (length(obs_na)>0){
        Y_obs_k[obs_na,j] <- apply(rmultinom(length(obs_na),1,prob[[k]][[j]]),2,which.max)
        }
    }
    Yres[obs_k,] <- Y_obs_k
  }
  
  return(Yres)
}






######
# Name: SimulClass_SEM_Quali
# Date: 2022/02/17
# Description:
# Arguments:
## YNA: matrix containing missing values (data matrix of size n,d)
## Y: current imputed matrix (data matrix of size n,d)
## prob: current value of the probabilities (list of length K, with list of length d of vectors of size the number of levels)
## alpha: current value of mechanism coefficients for the effect of the class membership (data matrix of size K,d)
## prop.pi: current proportion per class (vector of length K)
#####


SimulClass_SEM_Quali <- function(YNA,Y,prob,alpha,prop.pi){
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)),ncol=ncol(YNA))
  K <- length(prop.pi)
  
  Zres <- matrix(0,nrow=n,ncol=K)
  tau <- matrix(0,nrow=n,ncol=K)

  
  for (k in 1:K){
    int_pn <- matrix(alpha[k,],nrow=n,ncol=d,byrow=TRUE)
    pn <- pnorm(int_pn)
    logMecha <- rowSums(C*log(pn)+(1-C)*log(1-pn))
    logY <- 0
    for (j in 1:d){
      logY <- logY + log(prob[[k]][[j]][Y[,j]])
    }
    tau[,k] = logMecha + log(prop.pi[k]) + logY
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
# Name: Critere_Quali
# Date: 2022/02/17
# Description:
# Arguments:
## ResAlgo: object of the type returning by the function SEM_Gaussian (list)
## YNA: matrix containing missing values (data matrix of size n,d)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
## Partition_true: true partition (vector of size n)
#####


Critere_Quali <- function(ResAlgo,YNA,mecha,Partition_true=NULL){
  
  C <- is.na(YNA)
  
  Pattern <- matrix(C[!duplicated(C),],nrow=sum(!duplicated(C)))
  
  Y <- YNA 
  n <- nrow(Y)
  d <- ncol(Y)
  prob <- ResAlgo$prob
  pik <- ResAlgo$pik
  K <- length(pik)
  alpha <- ResAlgo$alpha
  
  Nb_levels <- sapply(1:d, function(j) length(unique(Y[!is.na(Y[,j]),j])))
  
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
    d_varobs <- length(var_obs)
    C_obs <- matrix(C_Pat[,var_obs],ncol=length(var_obs))
    
    if (nrow(Y_Pat)==1){
      Y_obs <- matrix(Y_Pat[,var_obs],nrow=nrow(Y_Pat),ncol=length(var_obs))
      C_obs <- matrix(C_Pat[,var_obs],nrow=nrow(Y_Pat),ncol=length(var_obs))
    }
    
    
    LogF <- matrix(NA,nrow=length(wherePat),ncol=K)
    for (ik in 1:K){
      term1_log = log(pik[ik]) +  rowSums(matrix(sapply(1:length(var_obs) , function(j) log(prob[[ik]][[var_obs[j]]][Y_obs[,j]])),ncol=length(var_obs)))
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
          res_miss <- log(pnorm(alpha[ik,j]))
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
  
  df <-  K-1 + K*sum(Nb_levels-1)
  if (mecha == "MCAR"){
    df <- df + 1
  }else if(mecha == "MNARz"){
    df <- df + K
  }else if(mecha == "MNARzj"){
    df <-  df + K*d
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







