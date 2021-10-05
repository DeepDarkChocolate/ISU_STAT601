getwd()
soil = read.csv("soilresptemp.txt", sep = " ")

Rsoil = soil$Rsoil
Tsoil = soil$Tsoil

Rsoil_rmv = Rsoil[Rsoil > 0]
Tsoil_rmv = Tsoil[Rsoil > 0]

getest = function(x, y, beta_ini, ..., eps = 10^(-10)){
  beta = beta_ini
  iter = 0
  
  while(TRUE){
    iter = iter + 1
    eta = x %*% beta
    mu = exp(eta)
    D = c(mu) * x
    
    #mu = ifelse(abs(mu) > eps, mu, eps)
    #mu = ifelse(abs(1 - mu) > eps, mu, 1 - eps)
    
    V_inv = diag(1 / c(mu))
    z = y - mu
    
    beta_2 = beta + solve(t(D) %*% V_inv %*% D, t(D) %*% V_inv %*% z)
    #xi_2 = solve(t(x) %*% W %*% x, t(x) %*% W %*% (x %*% xi + z))
    if(any(is.na(beta_2))){
      stop("Fatal error:: NA's generated")
    }
    
    if(norm(beta - beta_2) < 10^(-10)){
      beta = beta_2
      break
    }
    else{
      beta =  beta_2
    }
  }
  eta = x %*% beta
  mu = exp(eta)
  D = c(mu) * x
  V_inv = diag(1 / c(mu))
  phi = sum((y - mu)^2 / mu) / (length(mu) - 2)
  
  return(list(beta = beta, phi = phi, cov = solve(t(D) %*% V_inv %*% D) * phi, iteration = iter))
}

quasilik = function(beta, x, y){
  eta = x %*% beta
  mu = exp(eta)
  return(sum(y * log(mu / y) - mu + y))
}

quasilik_neg = function(beta, x, y){
  return(-quasilik(beta, x, y))
}

loglik = function(theta, x, y){
  gamma = theta[1:2]
  beta = theta[3]
  n = length(y)
  eta = x %*% gamma
  mu = exp(eta)
  alpha = mu * beta
  return(sum(sapply(1:n, function(i) dgamma(y[i], shape = alpha[i], rate = beta,  log = TRUE))))
}

loglik_neg = function(theta, x, y){
  return(-loglik(theta, x, y))
} 

grad = function(theta, x, y){
  gamma = theta[1:2]
  beta = theta[3]
  eta = x %*% gamma
  mu = exp(eta)
  alpha = mu * beta
  
  l_alpha = log(beta) - digamma(alpha) + log(y)
  l_beta = sum(mu - y)
  l_gamma = beta * t(x) %*% (l_alpha * mu)
  
  return(c(l_gamma, l_beta))
}


getest2 = function(x, y, theta_ini, ..., eps = 10^(-5)){
  gamma_ini = theta_ini[1:2]
  beta_ini = theta_ini[3]
  iter = 0
  theta = c(gamma_ini, beta_ini)
  
  while(TRUE){
    gamma = theta[1:2]
    beta = theta[3]
    iter = iter + 1
    eta = x %*% gamma
    mu = exp(eta)
    alpha = mu * beta
    
    l_alpha = log(beta) - digamma(alpha) + log(y)
    l_beta = sum(mu - y)
    l_gamma = beta * t(x) %*% (l_alpha * mu)
    
    l_theta = c(l_gamma, l_beta)
    
    l_alpha2 = - trigamma(alpha)
    l_beta2 = -sum(alpha^2 / beta)
    
    
    w11 = l_alpha2 * mu^2 * beta^2 + l_alpha * mu * beta
    
    l_gamma2 = t(x) %*% diag(c(w11)) %*% x
    l_gammabeta = t(x) %*% mu
    
    l_theta2 = rbind(cbind(l_gamma2, l_gammabeta), c(l_gammabeta, l_beta2))
    
    #print(iter)
    #print(l_theta)
    #print(l_theta2)
    
    theta_2 = c(gamma, beta) - solve(l_theta2, l_theta)
    
    logliktmp = loglik(theta, x, y)
    logliktmp2 = loglik(theta_2, x, y)
    
    if(any(is.na(theta_2))){
      stop("Fatal error:: NA's generated")
    }
    
    if(norm(theta - theta_2, "2") < eps & norm(logliktmp - logliktmp2, "2") < eps){
      theta = theta_2
      break
    }
    else{
      theta =  theta_2
    }
  }
  gamma = theta_2[1:2]
  beta = theta_2[3]
  logliktmp = loglik(theta, x, y)
  
  return(list(gamma = gamma, beta = beta, iteration = iter, gradient = l_theta, 
              loglik = logliktmp, info = -l_theta2))
}

#Estimates
x = cbind(1, Tsoil_rmv)
y = Rsoil_rmv

model_ini = glm(Rsoil_rmv ~ Tsoil_rmv, family = Gamma(link = "log"))
s_model_ini <- summary(model_ini)
beta_ini = model_ini$coefficients
#confint(model_ini, trace = TRUE)
moe = qnorm(1 - 0.025) *sqrt(diag(vcov(model_ini)))
cbind(beta_ini - moe, beta_ini, beta_ini + moe)

res = getest(x, y, beta_ini)
beta = res$beta

res2 = optim(par = beta_ini, fn = quasilik_neg, x = x, y = y)
beta2 = res2$par

moe = qnorm(1 - 0.025) * sqrt(diag(res$cov))
cbind(res$beta - moe, res$beta, res$beta + moe)


plot(Rsoil_rmv ~ Tsoil_rmv)

xvec = seq(from = min(Tsoil_rmv), to = max(Tsoil_rmv), len = 100)
xmat = cbind(1, xvec)
lines(exp(xmat %*% beta_ini)~ xvec , type = "l", col = "red")
lines(exp(xmat %*% beta)~ xvec , type = "l", col = "blue")
beta3 = c(0.06985412, 0.06838407)
lines(exp(xmat %*% beta3)~ xvec , type = "l", col = "green")
legend("topleft", legend=c("glm1", "quasi-likelihood", "glm2"),
       col=c("red", "blue", "green"), lty=1, cex=0.8)

mu = exp(x %*% beta)
plot((y - mu) / sqrt(mu)~y, main = "Residual Plot for quasi-likelihood",
     ylab = "Pearson residuals")

mu_ini = exp(x %*% beta_ini)
plot((y - mu_ini) / mu_ini~y, main = "Residual Plot for GLM",
     ylab = "Pearson residuals")
