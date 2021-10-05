getwd()
heartdat = read.csv("heartdat0.txt", sep = " ")


  
  
heartdat = read.csv("heartdat10_fixed.txt", sep = " ")
if(!is.matrix(heartdat)) heartdat = as.matrix(heartdat)
heartdat_uncen = heartdat[heartdat[,4] == 1,]
heartdat_cen = heartdat[heartdat[,4] == 0,]

loglik_cen = function(pars, heartdat){
  theta = pars[1:2]
  sigma = pars[3]
  if(length(heartdat) == 0) return(0)
  x = heartdat[,1:2]
  y = heartdat[,3]
  h = exp(x %*% theta)
  n = nrow(heartdat)
  z = (y / h)^(1 / sigma)
  return(sum(- z))
}

loglik_uncen = function(pars, heartdat){
  theta = pars[1:2]
  sigma = pars[3]
  if(length(heartdat) == 0) return(0)
  x = heartdat[,1:2]
  y = heartdat[,3]
  h = exp(x %*% theta)
  n = nrow(heartdat)
  zL = ((y - 10^(-3))/ h)^(1 / sigma)
  zU = ((y + 10^(-3))/ h)^(1 / sigma)
  #z = (y / h)^(1 / sigma)
  #return(sum(-log(sigma) + (1 / sigma - 1) * log(z) - z / sigma))
  return(sum(log(exp(-zL) - exp(-zU))))
}

loglik = function(pars, heartdat){
  heartdat_uncen = heartdat[heartdat[,4] == 1,]
  heartdat_cen = heartdat[heartdat[,4] == 0,]
  return(loglik_uncen(pars, heartdat_uncen) + 
           loglik_cen(pars, heartdat_cen))
} 

loglik_neg = function(pars, heartdat){
  return(-loglik(pars, heartdat))
}

theta_ini = unname(coef(lm(log(heartdat[,"y"]) ~ 0 + heartdat[,c("trt", "ldl")])))
sigma_vec = seq(from = 0.01, to = 10, length = 100)
loglik_vec = sapply(sigma_vec, function(k) loglik(c(theta_ini, k), heartdat))
#plot(sigma_vec, loglik_vec)
sigma_ini = sigma_vec[which.max(loglik_vec)]
pars_ini = c(theta_ini, sigma_ini)

res = optim(par = pars_ini, fn = loglik_neg, heartdat = heartdat, hessian = TRUE)
(pars = res$par)
loglik(res$par, heartdat)
#res$hessian #Information matrix
var_est = diag(solve(res$hessian)) #Estimated variance
ME = sqrt(var_est) * qnorm(1 - 0.025)
res_table = cbind(pars - ME, pars, pars + ME)
row.names(res_table) <- c("theta1", "theta2", "sigma")
colnames(res_table) <- c("95%LL", "Estimates", "95%UL")
res_table

