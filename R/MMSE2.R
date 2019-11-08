MMSE2 <-
  function(K, c, Add.Inf, thetahat, betahat,
           eig.cutoff = 1) {
    thetahat.length <- length(thetahat)
    ind1 <- 1:thetahat.length
    betahat.length <- nrow(K) - thetahat.length 
    ind2 <- (thetahat.length+1):(betahat.length+thetahat.length)
    thetahat.cov <- K[ind1,ind1]
    theta_delta_hat.cov <- K[ind1,ind2]
    delta_theta_hat.cov <- K[ind2,ind1]
    c2 <- c %o% c
    c2[Add.Inf$Biases == 0, ] <- 0
    c2[, Add.Inf$Biases == 0] <- 0
    
    eta.cov <- K[ind2,ind2] + c2
    
    K.add <- as.matrix(Matrix::bdiag(Add.Inf$Vars))
    betatilde <- unlist(Add.Inf$Means)
    SVD_var <- svd(eta.cov + K.add)
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
    thetahat.MMSE = thetahat - theta_delta_hat.cov %*% Kinv %*% (betahat - betatilde)
    thetahat.MMSE.MSE = thetahat.cov - theta_delta_hat.cov %*% Kinv %*% delta_theta_hat.cov
    list(Theta.Est = thetahat.MMSE, 
         Theta.Est.MSE = thetahat.MMSE.MSE,
         Theta.Hat = thetahat, Theta.Hat.Var = thetahat.cov)
}


