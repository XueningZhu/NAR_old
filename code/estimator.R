
################################### Estimator Function for NAR model
################################### ##########################################

### OLS estimation for theta: eq (2.8)
betaOLS <- function(Ymat, W, Z) {
    Ymat1 = W %*% Ymat  ### obtain WY
    Time = ncol(Ymat) - 1  ### Time doesn't count Y0
    ### the intercept, WY_{t-1}, Y_{t-1}, and nodal covariates
    X = cbind(rep(1, nrow(Ymat) * Time), as.vector(Ymat1[, -ncol(Ymat)]), 
        as.vector(Ymat[, -ncol(Ymat)]), do.call("rbind", rep(list(Z), Time)))
    invXX = solve(crossprod(X))  ### {t(X)X}^{-1}
    Yvec = as.vector(Ymat[, -1])  ### the response vector
    thetaEst = invXX %*% colSums(X * Yvec)  ### estimation equation (2.8)
    sigmaHat2 = mean((Yvec - X %*% thetaEst)^2)  ### estimation for hat sigma^2
    covHat = invXX * sigmaHat2  ### covariance for hat theta
    return(list(theta = thetaEst, covHat = covHat, sigmaHat = sqrt(sigmaHat2)))  ### return the result
} 
