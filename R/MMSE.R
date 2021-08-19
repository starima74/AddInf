MMSE <-
function(dd, theta.f, Add.Inf,
                nboots = 500, 
                eig.cutoff = 0.9, ...) {
  p <- nrow(Add.Inf)
  n <- nrow(dd)
  thetahat <- theta.f(dd, ...)
  thetahat.length <- length(thetahat)
  betahat <- c()
  for(i in 1:p) {
    betahat.add <- Add.Inf$Functions[[i]](dd)
    betahat <- c(betahat,betahat.add)
  }
  betahat.length <- length(betahat)
  bootres <- matrix(NA, nrow = nboots, 
                    ncol = thetahat.length + betahat.length)
  for(i in 1:nboots) {
    ind <- sample(1:n, replace = T)
    boot.hat <- theta.f(dd[ind,], ...)
    for(j in 1:p) 
      boot.hat <- c(boot.hat,Add.Inf$Functions[[j]](dd[ind,]))
    bootres[i,] <- boot.hat
  }
  K = var(bootres)/nrow(dd)
  ind1 <- 1:thetahat.length
  ind2 <- (thetahat.length+1):(betahat.length+thetahat.length)
  K11 <- K[ind1,ind1] * n
  K12 <- K[ind1,ind2] * n
  K21 <- K[ind2,ind1] * n
  if(is.null(dim(bootres[,ind2])))  Biases2 <- (mean(bootres[,ind2]) - unlist(Add.Inf$Means))^2
  if(!is.null(dim(bootres[,ind2]))) Biases2 <- (colMeans(bootres[,ind2]) - unlist(Add.Inf$Means))^2
  Biases2 <- Biases2 %o% Biases2
  Biases2[Add.Inf$Biases == 0, ] <- 0
  Biases2[, Add.Inf$Biases == 0] <- 0
  K22 <- K[ind2,ind2] * n + Biases2
  K.add <- as.matrix(Matrix::bdiag(Add.Inf$Vars))
  betatilde <- unlist(Add.Inf$Means)
  SVD_var <- svd(K22 + K.add)
  ind.svd <- rep(1,length(SVD_var$d))
  if(length(SVD_var$d)>1)
    for(i in 2:length(SVD_var$d))
      ind.svd[i] <- sum(SVD_var$d[1:(i-1)])/sum(SVD_var$d) < eig.cutoff
  diag_elem <- rep(0,length(SVD_var$d))
  diag_elem[SVD_var$d>0 & ind.svd == 1] <- 
    1/SVD_var$d[SVD_var$d>0 & ind.svd == 1]
  if(length(SVD_var$d>0 & ind.svd == 1) == 1)  
    Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
  if(length(SVD_var$d>0 & ind.svd == 1) > 1)  
    Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
  thetahat.MVAR = thetahat + K12 %*% Kinv %*% (betahat - betatilde)
  thetahat.MVAR.Var = K11 - K12 %*% Kinv %*% K21
  list(Theta.Est = thetahat.MVAR, 
       Theta.Est.Var = thetahat.MVAR.Var,
       Theta.Hat = thetahat, Theta.Hat.Var = K11)
}
