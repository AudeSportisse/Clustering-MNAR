SimulClass_SEM_Mixed <- function(YNA,Y,mu,sigma,prob,alpha,prop.pi,var_cont,var_cat){
  
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
    logYcat <- 0
    for (j in 1:length(var_cat)){
      logYcat <- logYcat + log(prob[[k]][[j]][Y[,var_cat[j]]])
    }
    tau[,k] = logMecha+dmvnorm(Y[,var_cont],mu[[k]],sigma[[k]],log=TRUE)+log(prop.pi[k]) + logYcat
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

  


SimulNA_SEM_Mixed <- function(YNA,mecha,Z,mu,sigma,prob,alpha,var_cont,var_cat){
  
  n <- dim(YNA)[1]
  K <- dim(Z)[2]
  
  resCont <- SimulNA_SEM_Gaussian(YNA[,var_cont],mecha,diag=TRUE,Z=Z,mu=mu,sigma=sigma,alpha=alpha[,var_cont],beta=matrix(0,nrow=K,ncol=var_cont),L=NULL)
  resCat <- SimulNA_SEM_Quali(YNA[,var_cat],Z,prob)
  
  Yres <- YNA
  if (typeof(resCont)=="character"){
    return("error: a matrix not symmetric is used as covariance matrix.")
  }else{
    Yres[,var_cat] <- resCat
    Yres[,var_cont] <- resCont 
  }
  
  return(Yres)
}






Critere_Mixed <- function(ResAlgo,YNA,mecha,type_data,Partition_true=NULL){
  
  C <- is.na(YNA)
  
  Pattern <- matrix(C[!duplicated(C),],nrow=sum(!duplicated(C)))
  
  
  var_cat <- which(type_data=="cat")
  var_cont <- which(type_data=="cont")
  d_cat <- length(var_cat)
  d_cont <- length(var_cont)
  index_var <- rowSums(sapply(unique(type_data),function(x) {cumsum(type_data==x)*(type_data==x)}))
  
  Nb_levels <- sapply(1:d_cat, function(j) length(unique(YNA[!is.na(YNA[,j]),j])))
  
  Y <- YNA 
  n <- nrow(Y)
  d <- ncol(Y)
  mu <- ResAlgo$mu
  sigma <- ResAlgo$sigma
  pik <- ResAlgo$pik
  K <- length(pik)
  prob <- ResAlgo$prob
  alpha <- ResAlgo$alpha
  
  
  valcrit <- 0
  Zall <- matrix(NA,ncol=K,nrow=n)
  proba_class <- matrix(NA,ncol=K,nrow=n)
  loglik=0
  
  for (i in 1:nrow(Pattern)){
    wherePat <- which(sapply(1:nrow(C),function(x) prod(C[x,]==Pattern[i,])==1))
    Y_Pat <- rbind(Y[wherePat,])
    C_Pat <- rbind(C[wherePat,])
    
    
    
    var_obs <- which(Pattern[i,]==0)
    
    var_obs_cat <- intersect(var_obs,var_cat)
    if (length(var_obs_cat)>0){
      Y_obs_cat <- matrix(Y_Pat[,var_obs_cat],ncol=length(var_obs_cat)) 
    }
    var_obs_cat_order <- index_var[var_obs_cat]
    
    var_obs_cont <- intersect(var_obs,var_cont)
    if (length(var_obs_cont)>0){
      Y_obs_cont <- matrix(Y_Pat[,var_obs_cont],ncol=length(var_obs_cont)) 
    }
    var_obs_cont_order <- index_var[var_obs_cont]
    
    mu_obs <- lapply(mu, function(x) x[var_obs_cont_order])
    sigma_obs <- lapply(sigma , function(x) matrix(x[var_obs_cont_order,var_obs_cont_order],ncol=length(var_obs_cont),nrow=length(var_obs_cont)))
    
    
    
    
    LogF <- matrix(NA,nrow=length(wherePat),ncol=K)
    for (ik in 1:K){
      if (length(var_obs_cat)>0){
        term1_log_multi =  rowSums(matrix(sapply(1:length(var_obs_cat) , function(j) log(prob[[ik]][[var_obs_cat_order[j]]][Y_obs_cat[,j]])),ncol=length(var_obs_cat))) 
      }else{
        term1_log_multi = 0
      }
      if (length(var_obs_cont)>0){
        term1_log_quanti = dmvnorm(Y_obs_cont,mean=t(mu_obs[[ik]]),sigma=sigma_obs[[ik]],log=TRUE) 
      }else{
        term1_log_quanti = 0
      }
      
      term1_log = log(pik[ik]) + term1_log_quanti + term1_log_multi
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
  
  df <- K-1 + K*sum(Nb_levels-1) + K*d_cont + K*d_cont*(d_cont+1)/2
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
