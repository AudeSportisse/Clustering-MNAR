library('ExtDist')

SimuC_miss <- function(pik,Y,Z,probmiss_z=NULL,opt="probit"){
   #######Arguments
  ##pik: vector of size K, proportion of the K classes
  ##Y: complete data matrix of size (n,d)
  ##Z: matrix of size (n,K) for the partition
  ##probiss_z: parameters for the missing-data mechanism, related to the class membership
  ##opt: link for the missing-data mechanism, probit, logit or Laplace
  #######Values
  ##Return an imcomplete data matrix by introducing missing values with the mechanism mecha
  
  if(is.null(probmiss_z)){
    return("error: no parameter for the mechanism is given.")
  }
  if (is.matrix(Y)==TRUE){
    n <- dim(Y)[1]
    d <- dim(Y)[2] 
  }else{
    n <- length(Y)
    d <- 1
    Y <- matrix(Y)
  }
  K <- length(pik)
  C = matrix(TRUE,ncol=d,nrow=n)
  
  if(opt == "probit"){
    for (k in 1:K){
      obs_k <- which(Z[,k]==1)
      C[obs_k,] <-   matrix(runif(length(obs_k)*d,0,1) < pnorm(probmiss_z[k]),ncol=d)
    }
  }else if(opt == "logit"){
    for (k in 1:K){
      obs_k <- which(Z[,k]==1)
      C[obs_k,] <-   matrix(runif(length(obs_k)*d,0,1) < 1/(1+exp(-probmiss_z[k])),ncol=d)
    }
  }else if(opt == "Laplace"){
    for (k in 1:K){
      obs_k <- which(Z[,k]==1)
      C[obs_k,] <-   matrix(runif(length(obs_k)*d,0,1) < dLaplace(probmiss_z[k]),ncol=d)
    }
  }
  #else if(opt == "uniform"){
    #for (k in 1:K){
      #obs_k <- which(Z[,k]==1)
      #ind1 = 0<=probmiss_z[k] & probmiss_z[k]<1
      #ind2 = probmiss_z[k]<0 | probmiss_z[k]>=1
      #C[obs_k,] <-   matrix(runif(length(obs_k)*d,0,1) < probmiss_z[k]*ind1 + 1*ind2,ncol=d)
    #}
  #}
  
  return(C)  
}

