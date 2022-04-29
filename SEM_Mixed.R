source("SEM_Quali.R")
source("SEM_Gaussian.R")
source("Internal_Functions_Mixed.R")


SEM_Mixed <- function(YNA, K, type_data, mecha, rmax = 100, init = NULL, stop = "rmax", samplesize = NULL) {
  
  ## When the type of data is the same for all the covariables
  if (length(unique(type_data))==1){
    if (unique(type_data)=="cont") {
      return(SEM_Gaussian(YNA, K, mecha, diag = TRUE, rmax = rmax, init = init, stop = stop, samplesize = samplesize))
    }else{
      return(SEM_Quali(YNA, K, mecha, rmax, init, stop, samplesize))
    }
  }
  
  ## Check if type_data is of length the number of variables
  if (length(type_data)!=ncol(YNA)) stop("Number of variables different from the number of variable types")
  
  
  var_cat <- which(type_data=="cat")
  var_cont <- which(type_data=="cont")
  d_cat <- length(var_cat)
  d_cont <- length(var_cont)
  index_var_type <- rowSums(sapply(unique(type_data),function(x) {cumsum(type_data==x)*(type_data==x)}))
  
  for (j in var_cat){
    YNA[,j] <- as.numeric(YNA[,j])
  }
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)),ncol=ncol(YNA))
  
  ## Number of possible values for the categorical variables
  Nb_levels <- sapply(var_cat, function(j) length(unique(YNA[!is.na(YNA[,j]),j])))
  
  ########################################
  #############INITIALIZATION#############
  ########################################

  init <- init_SEM_Mixed(YNA, K, mecha, init, samplesize, var_cont, var_cat)
  Z.new <- init$Z.init
  pi.new <- init$pi.init
  prob.new <- init$prob.init
  mu.new <- init$mu.init
  sigma.new <- init$sigma.init
  alpha.new <- init$alpha.init
  
  if (stop == "loglikmax"){
    pi.new.save <- pi.new
    prob.new.save <- prob.new
    mu.new.save <- mu.new
    sigma.new.save <- sigma.new
    alpha.new.save <- alpha.new
  }
  
  if (stop == "cut_average"){
    sumpar.pi <- rep(0,K)
    sumpar.prob <- lapply(1:K , function(z) lapply(1:d_cat , function(j) rep(0,Nb_levels[j])))
    sumpar.mu <- list()
    sumpar.sigma <- list()
    for (k in 1:K){
      sumpar.mu[[k]] <- rep(0,d_cont)
      sumpar.sigma[[k]] <- matrix(0,nrow=d,ncol=d_cont)
    }
    sumpar.alpha <- matrix(0,nrow=K,ncol=d)
  }
  
  if(sum(unlist(lapply(1:K,function(z) is.positive.definite(sigma.new[[z]]))))<K){
    return("error: an estimator of a covariance matrix is not positive definite")
  }
  
  if (mecha == "MNARz" | mecha == "MNARzj"){
    alpha.new <- Mechanism_SEM_GLM(YNA,YNA,Z.new,mecha)$alpha.new 
  }
  
  Y.new <- SimulNA_SEM_Mixed(YNA,mecha,Z.new,mu.new,sigma.new,prob.new,alpha.new,var_cont,var_cat)
  if (typeof(Y.new)=="character"){
    return(Y.new)
  }
  
  res <- list(pik=pi.new, mu=mu.new, sigma=sigma.new, prob=prob.new, alpha=alpha.new)
  crit <- Critere_Mixed(res,YNA,mecha,type_data)
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
    
    ###SE-step: one-Gibbs
    Z.new <- SimulClass_SEM_Mixed(YNA,Y.new,mu.new,sigma.new,prob.new,alpha.new,pi.new,var_cont,var_cat)
    if (typeof(Z.new)=="character" ){
      return(Z.new)
    }
    
    Y.new <- SimulNA_SEM_Mixed(YNA,mecha,Z.new,mu.new,sigma.new,prob.new,alpha.new,var_cont,var_cat)
    if (typeof(Y.new)=="character" ){
      return(Y.new)
    }
    
    ### M-step
    pi.new <- colSums(Z.new)/n 
    for (k in 1:K){
      Obs_classk <- which(Z.new[,k]==1)
      if (length(Obs_classk)<2){
        return("error: no enough observations for class k.")
      }else{
        mu.new[[k]] <- colMeans(Y.new[Obs_classk,var_cont])
        sigma.new[[k]] <- diag(apply(Y.new[Obs_classk,var_cont],2,var))
        if (sum(is.na(sigma.new[[k]]))>0){return("error: an estimator of a covariance matrix contains NA")}
        if (!is.positive.definite(sigma.new[[k]])){
          return("error: an estimator of a covariance matrix is not positive definite")
        }
      }
    }
    prob.new <- lapply(1:K , function(k) lapply(1:d_cat , function(j) sapply(1:Nb_levels[j] , function(l) sum(Y.new[Z.new[,k]==1,var_cat[j]]==l))/sum(Z.new[,k]==1)))
    
    if (mecha != "MCAR"){
      alpha.new <- Mechanism_SEM_GLM(YNA,Y.new,Z.new,mecha)$alpha.new
    }
    
    res <- list(pik=pi.new, mu=mu.new, sigma=sigma.new, prob=prob.new, alpha=alpha.new)
    crit <- Critere_Mixed(res,YNA,mecha,type_data)
    if (typeof(crit)=="character"){
      return(crit)
    }
    loglik_vec <- c(loglik_vec,crit$loglik)
    
    
    if (stop == "loglikmax"){
      if(crit$loglik>prec){
        pi.new.save <- pi.new
        mu.new.save <- mu.new
        sigma.new.save <- sigma.new
        prob.new.save <- prob.new
        alpha.new.save <- alpha.new
        prec <- crit$loglik
      }
    }
    
    if (stop == "cut_average" & r>floor(0.8*rmax)){
      sumpar.pi <- sumpar.pi + pi.new
      sumpar.mu <- lapply(1:K , function(z) sumpar.mu[[z]]+mu.new[[z]])
      sumpar.sigma <- lapply(1:K , function(z) sumpar.sigma[[z]]+sigma.new[[z]])
      sumpar.prob <- lapply(1:K , function(z) lapply(1:d_cat , function(j) (sumpar.prob[[z]][[j]] + prob.new[[z]][[j]])))
      sumpar.alpha <- sumpar.alpha + alpha.new
      sumpar.beta <- sumpar.beta + beta.new
    }
    
  }
  
  if (stop == "loglikmax"){
    pi.new <- pi.new.save
    mu.new <- mu.new.save
    sigma.new <- sigma.new.save
    prob.new <- prob.new.save
    alpha.new <- alpha.new.save
  }
  
  if (stop == "cut_average"){
    rcut <- floor(0.8*rmax)
    pi.new <- sumpar.pi/(rmax-rcut)
    mu.new <- lapply(1:K , function(z) sumpar.mu[[z]]/(rmax-rcut))
    sigma.new <- lapply(1:K , function(z) sumpar.sigma[[z]]/(rmax-rcut))
    prob.new <- lapply(1:K , function(z) lapply(1:d_cat , function(j) sumpar.prob[[z]][[j]]/(rmax-rcut)))
    alpha.new <- sumpar.alpha/(rmax-rcut)
  }
  
  
  return(list(Ycomp=Y.new, Zcomp=Z.new, pik=pi.new, mu=mu.new, sigma=sigma.new, prob=prob.new, alpha=alpha.new,  loglik_vec=loglik_vec))
  
}



