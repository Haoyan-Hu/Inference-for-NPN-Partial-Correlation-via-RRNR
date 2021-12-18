library(MASS)
library(glmnet)
#---AR covariance structure---
ARcov = function(p, rho){
    Cov = matrix(0, p, p)
    for (i in 1 : p){
        for (j in 1 : p){
            Cov[i, j] = rho^(abs(i - j))
        }
    }
    return(Cov)
}
#---Block diagonal covariance structure 2---
BD2 = function(p, k, tau2 = 1, rho){
    # k is the number of blocks; rho is a vector, giving the lower and upper limits of the coefficients
    C = matrix(0, p, p)
    d = p / k
    for (m in 1 : k){
        rhotemp = runif(1, rho[1], rho[2])
        for (i in ((m - 1) * d + 1) : (m * d)) {
            for (j in ((m - 1) * d + 1) : (m * d)){
                if (i == j) C[i, j] = tau2
                else C[i, j] = rhotemp
            }	
        }
    }
    return(C)
}
# Normal Transformation #
normal_transform<-function(x){
    n=length(x)
    empirical_cdf<-ecdf(x)  
    x_ecdf_hat<-empirical_cdf(x)
    for (i in 1:n){
        if (x_ecdf_hat[i]<1/n^2){x_ecdf_hat[i]=1/n^2}
        if (x_ecdf_hat[i]>1-1/n^2){x_ecdf_hat[i]=1-1/n^2}
    }
    return(qnorm(x_ecdf_hat, mean = 0, sd = 1))
}
# Power h(X)=X^3 
power3_trans<-function(x){x^3}
# Exponential h(X)=exp(X), same with exp
exp_trans<-function(x){exp(x)}
# Identity h(X)=X
identity_trans<-function(x){x}

rrnr_sim = function(n, p, covariance, trans_type, trans_func, R){
    tau = seq(0.1, 3.1, 0.01); alpha = c(0.05, 0.1, 0.15, 0.2)
    lentau = length(tau); lenalpha = length(alpha)
    
    # Using precision matrix
    precisionMatrix = solve(covariance)
    IndMatrix = matrix(1, p, p) - diag(rep(1, p))
    STrue = 1 * (abs(precisionMatrix) > 10^(-3)); NoNSTrue = 1 * (STrue == 0)
    STrueNum = sum(STrue) - p

    result = list()
    for (rep in 1:R){
        X = mvrnorm(n, rep(0, p), covariance)
        
        X_original<-X
        if (trans_type == "all"){ 
            X_original <- trans_func(X) # X_original : Z
        }else if (trans_type == "alternate"){
            for (j in seq(1,p,2)){ X_original[,j] <- trans_func(X[,j]) } 
        }
        
        X_est_transform<-X_original
        for (j in 1:p){ X_est_transform[,j]<-normal_transform(X_original[,j]) } 
        
        #_______________RRNR______________#
        {
            X = X_est_transform
            Eresidual = matrix(0, n, p)
            CoefMatrix = matrix(0, p, p - 1)
            
            # Standardization
            meanX = colMeans(X)
            X = t(t(X) - meanX)
            # XS = matrix(0, n, p)
            # for (i in 1 : p){
            #     XS[, i] = X[, i] / sd(X[, i])
            # }
            
            for (i in 1 : p){
                # glmnet will scale x variable by default, not y
                node_reg = glmnet(X[, -i], y = X[, i], lambda = 2*sqrt(log(p)/n), alpha =1)
                # X is centered, thus no intercept
                # residual is calculated at the original scale
                Eresidual[, i] = as.vector(X[,i] - X[, -i] %*% node_reg$beta )
                # the returned beta is at the original scale
                CoefMatrix[i, ] = as.vector(node_reg$beta)
            }
            
            # V matrix in paper, covariance of residuals
            CovRes = t(Eresidual) %*% Eresidual / n
            Est = matrix(1, p, p)
            
            for (i in 1 : (p - 1)){
                for (j in (i + 1) : p){
                    temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
                    Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
                    Est[j, i] = Est[i, j]
                }
            }
            
            EstThresh = Est * ( abs(Est) >= (2 * sqrt(log(p) / n) * IndMatrix) )
            
            CovX = t(X) %*% X / n - matrix(colMeans(X), p, 1) %*% matrix(colMeans(X), 1, p)
            
            res = c(); reject = c()
            # Need result, FDPresult, resultC0, FDPresultC0
            for (i in 1 : lentau){
                Threshold = tau[i] * sqrt(log(p) / n) * (1 - EstThresh^2)
                # Not C0 method
                SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0)
                TruePositive = sum(SRec * STrue) - p; FasleNegative = sum(NoNSRec * STrue); FalsePositive = sum(SRec * NoNSTrue); TrueNegative = sum(NoNSRec * NoNSTrue)
                temp = c(FalsePositive, FasleNegative, TruePositive, TrueNegative)
                #reject : # of reject - p
                reject = c(reject, max(1, (sum(SRec) - p)))
                res = c(res, temp)
            }
            # result
            result_est_transform = c(STrueNum, res)

            # FDPresult
            FDP = 2 * p * (p - 1) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / reject
            FDPres = c()
            for (i in 1 : lenalpha){
                if (sum(FDP <= alpha[i]) > 0) tau0 = min(c(2, tau[FDP <= alpha[i]]))
                else tau0 = 2
                Threshold = tau0 * sqrt(log(p) / n) * (1 - EstThresh^2)
                SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0)
                TruePositive = sum(SRec * STrue) - p; FalsePositive = sum(SRec * NoNSTrue)
                power = TruePositive / STrueNum
                FDPalpha = FalsePositive / max(1, (sum(SRec) - p))
                tempFDP = c(tau0, FDPalpha, power)
                FDPres = c(FDPres, tempFDP)
            }
            FDPresult_est_transform = FDPres
        }
        result[[rep]] = list(result_est_transform = result_est_transform, FDPresult_est_transform = FDPresult_est_transform)
    }
    return(result)
}

# Example: simulation with R = 1000 repetitions
# True transformation: exponential transformation on all random variables.
n = 60; p = 100; 
cov_ind = 'AR'
if (cov_ind == 'AR'){ AR_coef = 0.5; Sigma = ARcov(p, AR_coef)} else{ AR_coef = ''; Sigma = BD2(p, p/4, 1, c(0.3, 0.9)) }
true_transformation_type = 'all'
true_transformation_func_name = 'exp'
true_transformation_func = get(paste0(true_transformation_func_name,'_trans'))

result = rrnr_sim(n, p, Sigma, true_transformation_type, true_transformation_func, R = 1000)
