Impute0 <- function(X){
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[2]){
      if (is.na(X[i,j]==TRUE)){
        X[i,j]=0 
      }
    }
  }
  return(X)
}

MSE <- function(X, Xtrue) {
  return(sqrt(sum((X-Xtrue) ^ 2)))
}


imputeMean <- function(df, args_list = NULL){
  df <- data.frame(df)
  df_imp <- df
  vars_real <- colnames(df)[sapply(df, is.numeric)]
  vars_factor <- colnames(df)[!sapply(df, is.numeric)]
  if (length(vars_real) >0){
    if (length(vars_real)>1){
      mean_real <- sapply(df[,vars_real], mean, na.rm=T)
    } else {
      mean_real <- c(mean(df[,vars_real], na.rm=T))
    }
    df_imp[ , vars_real] <- sapply(1:length(vars_real), 
                                   function(x) ifelse(is.na(df_imp[,vars_real[x]]), mean_real[x], df_imp[,vars_real[x]]))
  }
  
  if (length(vars_factor) > 0 ){
    fact.levels <- sapply(vars_factor, FUN=function(x) levels(df[,x]))
    mode_factor <- as.vector(apply(data.frame(df[, vars_factor]), MARGIN=2, 
                                   FUN = function(x) as.data.frame(table(x))[which.max(table(x)),1]))
    df_imp[ , vars_factor] <- sapply(1:length(vars_factor), 
                                     function(x) ifelse(is.na(df_imp[,vars_factor[x]]), mode_factor[x], as.character(df_imp[,vars_factor[x]])))
    for (i in 1:length(vars_factor)){
      df_imp[,vars_factor[i]] <- as.factor(df_imp[,vars_factor[i]])
      levels(df_imp[,vars_factor[i]]) <- levels(df[,vars_factor[i]])
    }
  }
  return(df_imp)
}