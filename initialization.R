######
# Name: init_SEM_Gaussian
# Date: 2021/12/21
# Description: 
# Arguments: 
## YNA: matrix containing missing values (data matrix)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## init: choice of the initialization (NULL, "sample" or list)
## samplesize: size of the data sample used for the initialization (numeric)
#####


init_SEM_Gaussian <- function(YNA, K, mecha, diag, init, samplesize){

  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  
  
  ################################################################
  #############DEFAULT INITIALIZATION USING FULL DATA#############
  ################################################################
  
  if (is.null(init)){
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    
    if (K==1){
      Z.init=t(Z.init)
    }
    
    Obs_classk <- lapply(1:K, function(z){return(which(Z.init[,z]==1))})
    pi.init = colSums(Z.init)/n 
    
    mu.init <- lapply(1:K, function(z){return(colMeans(YNA[Obs_classk[[z]],],na.rm=T))})
    
    if (diag == TRUE){
      sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[Obs_classk[[z]],],2,var,na.rm=T)))})
    }else{
      sigma.init <- lapply(1:K, function(z){return(cov(YNA[Obs_classk[[z]],],use="complete.obs"))}) 
    }
    
    if (mecha == "MCAR"){
      alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
    }else{
      alpha.init <- matrix(0,nrow=K,ncol=d) 
    }
    
    
  ############################################################
  #############INITIALIZATION USING A DATA SAMPLE#############
  ############################################################
    
  }else if(typeof(init)=="character"){
    
    if (init!="sample"){
      stop("Type of initialization unknwown.") 
    }

    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    if (K==1){
      Z.init=t(Z.init)
    }
    
    if (is.null(samplesize)){
      samplesize = floor(n/5)
    }
    
    if (K==1){
      sampleUsed <- sample(1:n,samplesize)
      pi.init = sum(Z.init[sampleUsed])/samplesize
      mu.init <- lapply(1:K, function(z){return(colMeans(YNA[sampleUsed,],na.rm=T))})
      
      if (diag == TRUE){
        sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[sampleUsed,],2,var,na.rm=T)))})
      }else{
        sigma.init <- lapply(1:K, function(z){return(cov(YNA[sampleUsed,],use="complete.obs"))})
      }
      
      if (mecha == "MCAR"){
        alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }else{
      
      sampleUsed <- sample(1:n,samplesize)
      Z.init_sampleUsed <- Z.init[sampleUsed,]
      Obs_classk <- lapply(1:K, function(z){return(which(Z.init_sampleUsed[,z]==1))})
      pi.init = colSums(Z.init_sampleUsed)/samplesize
      
      mu.init <- lapply(1:K, function(z){return(colMeans(YNA[Obs_classk[[z]],],na.rm=T))})
      
      if (diag == TRUE){
        sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[Obs_classk[[z]],],2,var,na.rm=T)))})
      }else{
        sigma.init <- lapply(1:K, function(z){return(cov(YNA[Obs_classk[[z]],],use="complete.obs"))}) 
      }
      
      if (mecha == "MCAR"){
        alpha.init <- matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }
    
  ##########################################################
  #############INITIALIZATION GIVEN BY THE USER#############
  ##########################################################
    
  }else{
    
    if (!is.null(init$Zcomp)){
      Z.init <- init$Zcomp
    }else{
      Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
        res <- rep(0, K)
        res[i] <- 1
        return(res)
      }))
      if (K==1){
        Z.init=t(Z.init)
      }
    }
    pi.init <- init$pik
    mu.init <- init$mu
    sigma.init <- init$sigma
    alpha.init <- init$alpha
    
  }


  if (mecha == "MNARy" | mecha == "MNARyz" | mecha == "MNARyk" | mecha == "MNARykz" | mecha == "MNARyzj" | mecha == "MNARykzj"){
    
    beta.init <- matrix(0.1,nrow=K,ncol=d) 
    
    L.init <- matrix(0,nrow=n,ncol=d)
    for (i in 1:n){
      k <- which(Z.init[i,]==1)
      upper_trunc = numeric(d)
      lower_trunc = numeric(d)
      for (j in 1:d){
        if (is.na(YNA[i,j])==1){
          upper_trunc[j] = Inf
          lower_trunc[j] = 0
        }else{
          upper_trunc[j] = 0
          lower_trunc[j] = -Inf
        }
      }
      L.init[i,] <- rtmvnorm(1,alpha.init[k,]+beta.init[k,]*mu.init[[k]],lower=lower_trunc,upper=upper_trunc,algorithm="gibbs") 
    }
  
  }else{
    
    beta.init <- matrix(0,nrow=K,ncol=d)
    
    L.init <- NULL
  
  }
  
  return(list(Z.init=Z.init,pi.init=pi.init,mu.init=mu.init,sigma.init=sigma.init,alpha.init=alpha.init,beta.init=beta.init,L.init=L.init))
  
}



######
# Name: init_EM_Gaussian
# Date: 2021/12/21
# Description: 
# Arguments: 
## YNA: matrix containing missing values (data matrix)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR", "MNARz", "MNARy" or "MNARyz" (character)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## init: choice of the initialization (NULL, "sample" or list)
## samplesize: size of the data sample used for the initialization (numeric)
#####


init_EM_Gaussian <- function(YNA, K, mecha, diag, init, samplesize){
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  
  ################################################################
  #############DEFAULT INITIALIZATION USING FULL DATA#############
  ################################################################
  
  if (is.null(init)){
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    
    if (K==1){
      Z.init=t(Z.init)
    }
    
    Obs_classk <- lapply(1:K, function(z){return(which(Z.init[,z]==1))})
    pi.init = colSums(Z.init)/n 
    
    mu.init <- lapply(1:K, function(z){return(colMeans(YNA[Obs_classk[[z]],],na.rm=T))})
    
    if (diag == TRUE){
      sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[Obs_classk[[z]],],2,var,na.rm=T)))})
    }else{
      sigma.init <- lapply(1:K, function(z){return(cov(YNA[Obs_classk[[z]],],use="complete.obs"))}) 
    }
    
    if (mecha == "MCAR"){
      alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
    }else{
      alpha.init <- matrix(0,nrow=K,ncol=d) 
    }
    
    
    ############################################################
    #############INITIALIZATION USING A DATA SAMPLE#############
    ############################################################
    
  }else if(typeof(init)=="character"){
    
    if (init!="sample"){
      stop("Type of initialization unknwown.") 
    }
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    if (K==1){
      Z.init=t(Z.init)
    }
    
    if (is.null(samplesize)){
      samplesize = floor(n/5)
    }
    
    if (K==1){
      sampleUsed <- sample(1:n,samplesize)
      pi.init = sum(Z.init[sampleUsed])/samplesize
      mu.init <- lapply(1:K, function(z){return(colMeans(YNA[sampleUsed,],na.rm=T))})
      
      if (diag == TRUE){
        sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[sampleUsed,],2,var,na.rm=T)))})
      }else{
        sigma.init <- lapply(1:K, function(z){return(cov(YNA[sampleUsed,],use="complete.obs"))})
      }
      
      if (mecha == "MCAR"){
        alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }else{
      
      sampleUsed <- sample(1:n,samplesize)
      Z.init_sampleUsed <- Z.init[sampleUsed,]
      Obs_classk <- lapply(1:K, function(z){return(which(Z.init_sampleUsed[,z]==1))})
      pi.init = colSums(Z.init_sampleUsed)/samplesize
      
      mu.init <- lapply(1:K, function(z){return(colMeans(YNA[Obs_classk[[z]],],na.rm=T))})
      
      if (diag == TRUE){
        sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[Obs_classk[[z]],],2,var,na.rm=T)))})
      }else{
        sigma.init <- lapply(1:K, function(z){return(cov(YNA[Obs_classk[[z]],],use="complete.obs"))}) 
      }
      
      if (mecha == "MCAR"){
        alpha.init <- matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }
    
    ##########################################################
    #############INITIALIZATION GIVEN BY THE USER#############
    ##########################################################
    
  }else{
    
    pi.init <- init$pik
    mu.init <- init$mu
    sigma.init <- init$sigma
    alpha.init <- init$alpha
    
  }
  
  
  return(list(pi.init=pi.init,mu.init=mu.init,sigma.init=sigma.init,alpha.init=alpha.init))
  
  
}




######
# Name: init_EM_Gaussian_withoutNA
# Date: 2021/12/21
# Description: 
# Arguments: 
## Y: complete matrix (data matrix)
## K: number of classes to use (numeric)
## diag: if TRUE, the variables are independent given the class membership (boolean)
## init: choice of the initialization (NULL, "sample" or list)
## samplesize: size of the data sample used for the initialization (numeric)
#####


init_EM_Gaussian_withoutNA <- function(Y, K, diag, init, samplesize){
  
  n <- dim(Y)[1]
  d <- dim(Y)[2]
  
  ################################################################
  #############DEFAULT INITIALIZATION USING FULL DATA#############
  ################################################################
  
  if (is.null(init)){
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    
    if (K==1){
      Z.init=t(Z.init)
    }
    
    Obs_classk <- lapply(1:K, function(z){return(which(Z.init[,z]==1))})
    pi.init = colSums(Z.init)/n 
    
    mu.init <- lapply(1:K, function(z){return(colMeans(Y[Obs_classk[[z]],]))})
    
    if (diag == TRUE){
      sigma.init <- lapply(1:K, function(z){return(diag(apply(Y[Obs_classk[[z]],],2,var)))})
    }else{
      sigma.init <- lapply(1:K, function(z){return(cov(Y[Obs_classk[[z]],]))}) 
    }
    
    
    
    ############################################################
    #############INITIALIZATION USING A DATA SAMPLE#############
    ############################################################
    
  }else if(typeof(init)=="character"){
    
    if (init!="sample"){
      stop("Type of initialization unknwown.") 
    }
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    if (K==1){
      Z.init=t(Z.init)
    }
    
    if (is.null(samplesize)){
      samplesize = floor(n/5)
    }
    
    if (K==1){
      sampleUsed <- sample(1:n,samplesize)
      pi.init = sum(Z.init[sampleUsed])/samplesize
      mu.init <- lapply(1:K, function(z){return(colMeans(Y[sampleUsed,]))})
      
      if (diag == TRUE){
        sigma.init <- lapply(1:K, function(z){return(diag(apply(Y[sampleUsed,],2,var)))})
      }else{
        sigma.init <- lapply(1:K, function(z){return(cov(Y[sampleUsed,]))})
      }
      

    }else{
      
      sampleUsed <- sample(1:n,samplesize)
      Z.init_sampleUsed <- Z.init[sampleUsed,]
      Obs_classk <- lapply(1:K, function(z){return(which(Z.init_sampleUsed[,z]==1))})
      pi.init = colSums(Z.init_sampleUsed)/samplesize
      
      mu.init <- lapply(1:K, function(z){return(colMeans(Y[Obs_classk[[z]],]))})
      
      if (diag == TRUE){
        sigma.init <- lapply(1:K, function(z){return(diag(apply(Y[Obs_classk[[z]],],2,var)))})
      }else{
        sigma.init <- lapply(1:K, function(z){return(cov(Y[Obs_classk[[z]],]))}) 
      }
      

    }
    
    ##########################################################
    #############INITIALIZATION GIVEN BY THE USER#############
    ##########################################################
    
  }else{
    
    pi.init <- init$pik
    mu.init <- init$mu
    sigma.init <- init$sigma
    alpha.init <- init$alpha
    
  }
  
  
  return(list(pi.init=pi.init,mu.init=mu.init,sigma.init=sigma.init))
  
  
}



######
# Name: init_SEM_Quali
# Date: 2022/02/17
# Description: 
# Arguments: 
## YNA: matrix containing missing values (data matrix)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR", "MNARz" or "MNARzj" (character)
## init: choice of the initialization (NULL, "sample" or list)
## samplesize: size of the data sample used for the initialization (numeric)
#####


init_SEM_Quali <- function(YNA, K, mecha, diag, init, samplesize){
  
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)), ncol = ncol(YNA))
  
  Nb_levels <- sapply(1:d, function(j) length(unique(YNA[!is.na(YNA[,j]),j])))
  
  ################################################################
  #############DEFAULT INITIALIZATION USING FULL DATA#############
  ################################################################
  
  if (is.null(init)){
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    
    if (K==1){
      Z.init=t(Z.init)
    }
    
    Obs_classk <- lapply(1:K, function(z){return(which(Z.init[,z]==1))})
    pi.init = colSums(Z.init)/n 

    prob.init <- lapply(1:K , function(z) lapply(1:d , function(j) {
      t <- sapply(1:Nb_levels[j] , function(d) sum(YNA[Obs_classk[[z]],j]==d,na.rm=TRUE))/sum(1-C[Obs_classk[[z]],j])
      return((t+0.05)/sum(t+0.05))
    }))
    
    
    if (mecha == "MCAR"){
      alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
    }else{
      alpha.init <- matrix(0,nrow=K,ncol=d) 
    }
    
    
    ############################################################
    #############INITIALIZATION USING A DATA SAMPLE#############
    ############################################################
    
  }else if(typeof(init)=="character"){
    
    if (init!="sample"){
      stop("Type of initialization unknwown.") 
    }
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    if (K==1){
      Z.init=t(Z.init)
    }
    
    if (is.null(samplesize)){
      samplesize = floor(n/5)
    }
    
    if (K==1){
      sampleUsed <- sample(1:n,samplesize)
      pi.init = sum(Z.init[sampleUsed])/samplesize
      prob.init <- lapply(1:K , function(z) lapply(1:d , function(j) {
        t <- sapply(1:Nb_levels[j] , function(d) sum(YNA[sampleUsed,j]==d,na.rm=TRUE))/sum(1-C[sampleUsed,j])
        return((t+0.05)/sum(t+0.05))
      }))
      
      if (mecha == "MCAR"){
        alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }else{
      
      sampleUsed <- sample(1:n,samplesize)
      Z.init_sampleUsed <- Z.init[sampleUsed,]
      Obs_classk <- lapply(1:K, function(z){return(which(Z.init_sampleUsed[,z]==1))})
      pi.init = colSums(Z.init_sampleUsed)/samplesize
      
      
      prob.init <- lapply(1:K , function(z) lapply(1:d , function(j) {
        t <- sapply(1:Nb_levels[j] , function(d) sum(YNA[Obs_classk[[z]],j]==d,na.rm=TRUE))/sum(1-C[Obs_classk[[z]],j])
        return((t+0.05)/sum(t+0.05))
      }))
      
      if (mecha == "MCAR"){
        alpha.init <- matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }
    
    ##########################################################
    #############INITIALIZATION GIVEN BY THE USER#############
    ##########################################################
    
  }else{
    
    if (!is.null(init$Zcomp)){
      Z.init <- init$Zcomp
    }else{
      Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
        res <- rep(0, K)
        res[i] <- 1
        return(res)
      }))
      if (K==1){
        Z.init=t(Z.init)
      }
    }
    pi.init <- init$pik
    prob.init <- init$prob
    alpha.init <- init$alpha
    
  }

  
  return(list(Z.init=Z.init,pi.init=pi.init,prob.init=prob.init,alpha.init=alpha.init))
  
}


######
# Name: init_SEM_Mixed
# Date: 2021/02/17
# Description: 
# Arguments: 
## YNA: matrix containing missing values (data matrix)
## K: number of classes to use (numeric)
## mecha: mechanism to use, "MCAR", "MNARz" or "MNARzj" (character)
## init: choice of the initialization (NULL, "sample" or list)
## samplesize: size of the data sample used for the initialization (numeric)
## var_cont: indexes of the quantitative variables
## var_cat: indexes of the categorical variables
#####


init_SEM_Mixed <- function(YNA, K, mecha, init, samplesize, var_cont, var_cat){
  
  
  n <- dim(YNA)[1]
  d <- dim(YNA)[2]
  C <- matrix(as.numeric(is.na(YNA)), ncol = ncol(YNA))
  
  Nb_levels <- sapply(var_cat, function(j) length(unique(YNA[!is.na(YNA[,j]),j])))
  
  
  ################################################################
  #############DEFAULT INITIALIZATION USING FULL DATA#############
  ################################################################
  
  if (is.null(init)){
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    
    if (K==1){
      Z.init=t(Z.init)
    }
    
    Obs_classk <- lapply(1:K, function(z){return(which(Z.init[,z]==1))})
    pi.init = colSums(Z.init)/n 
    
    prob.init <- lapply(1:K , function(z) lapply(1:d_cat , function(j) {
      t <- sapply(1:Nb_levels[j] , function(d) sum(YNA[Obs_classk[[z]],var_cat[j]]==d,na.rm=TRUE))/sum(1-C[Obs_classk[[z]],var_cat[j]])
      return((t+0.05)/sum(t+0.05))
    }))
    
    mu.init <- lapply(1:K, function(z){return(colMeans(YNA[Obs_classk[[z]],var_cont],na.rm=T))})
    
    sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[Obs_classk[[z]],var_cont],2,var,na.rm=T)))})
    
    
    if (mecha == "MCAR"){
      alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
    }else{
      alpha.init <- matrix(0,nrow=K,ncol=d) 
    }
    
    
    ############################################################
    #############INITIALIZATION USING A DATA SAMPLE#############
    ############################################################
    
  }else if(typeof(init)=="character"){
    
    if (init!="sample"){
      stop("Type of initialization unknwown.") 
    }
    
    Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
      res <- rep(0, K)
      res[i] <- 1
      return(res)
    }))
    if (K==1){
      Z.init=t(Z.init)
    }
    
    if (is.null(samplesize)){
      samplesize = floor(n/5)
    }
    
    if (K==1){
      sampleUsed <- sample(1:n,samplesize)
      pi.init = sum(Z.init[sampleUsed])/samplesize
      prob.init <- lapply(1:K , function(z) lapply(1:d_cat , function(j) {
        t <- sapply(1:Nb_levels[j] , function(d) sum(YNA[sampleUsed,var_cat[j]]==d,na.rm=TRUE))/sum(1-C[sampleUsed,var_cat[j]])
        return((t+0.05)/sum(t+0.05))
      }))
      
      mu.init <- lapply(1:K, function(z){return(colMeans(YNA[sampleUsed,var_cont],na.rm=T))})
      
      sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[sampleUsed,var_cont],2,var,na.rm=T)))})
     
      if (mecha == "MCAR"){
        alpha.init <-  matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }else{
      
      sampleUsed <- sample(1:n,samplesize)
      Z.init_sampleUsed <- Z.init[sampleUsed,]
      Obs_classk <- lapply(1:K, function(z){return(which(Z.init_sampleUsed[,z]==1))})
      pi.init = colSums(Z.init_sampleUsed)/samplesize
      prob.init <- lapply(1:K , function(z) lapply(1:d_cat , function(j) {
        t <- sapply(1:Nb_levels[j] , function(d) sum(YNA[Obs_classk[[z]],var_cat[j]]==d,na.rm=TRUE))/sum(1-C[Obs_classk[[z]],var_cat[j]])
        return((t+0.05)/sum(t+0.05))
      }))
      
      mu.init <- lapply(1:K, function(z){return(colMeans(YNA[Obs_classk[[z]],var_cont],na.rm=T))})
      
      sigma.init <- lapply(1:K, function(z){return(diag(apply(YNA[Obs_classk[[z]],var_cont],2,var,na.rm=T)))})
      
      
      if (mecha == "MCAR"){
        alpha.init <- matrix(sum(is.na(YNA))/(n*d),nrow=K,ncol=d)
      }else{
        alpha.init <- matrix(0,nrow=K,ncol=d) 
      }
      
    }
    
    ##########################################################
    #############INITIALIZATION GIVEN BY THE USER#############
    ##########################################################
    
  }else{
    
    if (!is.null(init$Zcomp)){
      Z.init <- init$Zcomp
    }else{
      Z.init <- t(sapply(sample(1:K, n, replace = T) , function(i) {
        res <- rep(0, K)
        res[i] <- 1
        return(res)
      }))
      if (K==1){
        Z.init=t(Z.init)
      }
    }
    pi.init <- init$pik
    prob.init <- init$prob
    mu.init <- init$mu
    sigma.init <- init$sigma
    alpha.init <- init$alpha
    
  }
  
  

  return(list(Z.init=Z.init,pi.init=pi.init,mu.init=mu.init,sigma.init=sigma.init,prob.init=prob.init,alpha.init=alpha.init))
  
}




