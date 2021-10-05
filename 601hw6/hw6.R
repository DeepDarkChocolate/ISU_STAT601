y = read.csv("gammadat_assign6.txt", sep = " ")
y = y$y
n = length(y)

loglik = function(theta, y){
  alpha = theta[1]
  beta = theta[2]
  if(alpha <= 0 | beta <= 0) return(-Inf)
  n = length(y)
  #if(dgamma(y[i], shape = alpha, rate = beta,  log = TRUE)) print()
  return(sum(sapply(1:n, function(i) dgamma(y[i], shape = alpha, rate = beta,  log = TRUE))))
}

loglik_slice = function(alpha, beta, y){
  if(alpha <= 0 | beta <= 0) return(-Inf)
  n = length(y)
  #if(dgamma(y[i], shape = alpha, rate = beta,  log = TRUE)) print()
  return(sum(sapply(1:n, function(i) dgamma(y[i], shape = alpha, rate = beta,  log = TRUE))))
}

loglik_neg = function(theta, y){
  return(-loglik(theta, y))
} 

loglik_slice_neg = function(alpha, beta, y){
  return(-loglik_slice(alpha, beta, y))
} 


alphahat = mean(y)^2 / var(y)
betahat = mean(y) / var(y)

theta_ini = c(alphahat, betahat)

res = optim(par = theta_ini, fn = loglik_neg, y = y)

res$par
alphahat = res$par[1]
betahat = res$par[2]

mean(rgamma(1000, alphahat, rate = betahat))
var(rgamma(1000, alphahat, rate = betahat))

B = 2 * alphahat
lambda = 1  
gamma = lambda * betahat

#alphavec = seq(from = 0.1, to = B, length = 3000)
#sapply(alphavec, function(x) betahat^(n*x) / gamma(x)^n * prod(y)^(x - 1))
#plot(alphavec , sapply(alphavec, function(x) betahat^(n*x) / gamma(x)^n * prod(y)^(x - 1)))
#plot(alphavec , sapply(alphavec, function(x) n*x*log(betahat) - n * log(gamma(x)) + (x - 1) * sum(log(y))))

loglik_alpha = function(x, beta_t, ..., y2 = y, n2 = n){
  return(n2*x*log(beta_t) - n2 * log(gamma(x)) + (x - 1) * sum(log(y2)))
}

lik_alpha = function(x, beta_t, alpha_max, ..., y2 = y, n2 = n){
  if(x <= 0 | x >= B) return(0)
  const = integrate(lik_alpha_tmp, beta_t = beta_t, alpha_max = alpha_max, y2 = y, n2 = n, lower = 0, upper = B)$value
  return(exp(loglik_alpha(x, beta_t, y2 = y, n2 = n) - loglik_alpha(alpha_max, beta_t, y2 = y, n2 = n)) / const)
}

lik_alpha_tmp = function(x, beta_t, alpha_max, ..., y2 = y, n2 = n){
  return(exp(loglik_alpha(x, beta_t, y2 = y, n2 = n) - loglik_alpha(alpha_max, beta_t, y2 = y, n2 = n)))
}

optimize(interval = c(0, B), loglik_slice_neg, y = y, beta = betahat)

Gibbs <- function(alpha_ini, beta_ini, MCsize1, MCsize2){
  alpha_t = alpha_ini
  beta_t = beta_ini
  alpha_sample = c()
  beta_sample = c()
  #alphavec = seq(from = 0.1, to = B, length = 300)
  for(i in 1:MCsize1){
    #max_idx = which.max(sapply(alphavec, function(x) loglik_alpha(x, beta_t, y2 = y, n2 = n)))
    optim_res = optimize(interval = c(0, B), loglik_slice_neg, y = y, beta = beta_t)
    mu_alpha = optim_res$minimum

    #mu_alpha = alphavec[max_idx]
    sigma2_alpha = 1 / (lik_alpha(mu_alpha, beta_t, mu_alpha, y2 = y, n2 = n))^2 / 2 / pi
    for(j in 1:MCsize2){
      alpha_star = rnorm(1, mu_alpha, sqrt(sigma2_alpha))
      R = lik_alpha(alpha_star, beta_t, mu_alpha, y2 = y, n2 = n) / lik_alpha(alpha_t, beta_t, mu_alpha, y2 = y, n2 = n) *
        dnorm(alpha_t, mu_alpha, sqrt(sigma2_alpha)) / dnorm(alpha_star, mu_alpha, sqrt(sigma2_alpha))
      R = min(c(R, 1))
      
      alpha_t <- ifelse(runif(1) < R, alpha_star, alpha_t)
      #print(alpha_t)
      if(is.na(alpha_t)) stop("alpha_t is NA")
    }
    beta_t <- rgamma(1, shape = n * alpha_t + gamma, rate = sum(y) + lambda)
    alpha_sample <- c(alpha_sample, alpha_t)
    beta_sample <- c(beta_sample, beta_t)
  }
  return(list(alpha_sample, beta_sample))
}
res <- Gibbs(alphahat, betahat, 50000, 1)
res <- Gibbs(alphahat, betahat, 2, 5)
res <- Gibbs(alphahat, betahat, 10, 10)

res <- Gibbs(alphahat, betahat, 10000, 1)
#res <- Gibbs(alphahat, betahat, 1000, 10)
alpha_sampled <- res[[1]][-(1:100)]
beta_sampled <- res[[2]][-(1:100)]

alphahat
mean(alpha_sampled)
betahat
mean(beta_sampled)
mean(res[[1]])
mean(res[[2]])

plot(res[[1]], type = "l")
plot(res[[2]], type = "l")

acf(res[[1]], lag.max = 1000)
acf(res[[2]], lag.max = 1000)

loglik(c(mean(res[[1]]), mean(res[[2]])), y)
loglik(c(alphahat, betahat), y)

res <- Gibbs(alphahat, betahat, 500, 10)

Gibbs(5.815472, 0.5321557, 500, 10)

alpha_star = rnorm(1, mu_alpha, sqrt(sigma2_alpha))
alpha_t = mu_alpha

R = lik_alpha(alpha_star, beta_t, mu_alpha) / lik_alpha(alpha_t, beta_t, mu_alpha) *
  dnorm(alpha_t, mu_alpha, sqrt(sigma2_alpha)) / dnorm(alpha_star, mu_alpha, sqrt(sigma2_alpha))
R = min(c(R, 1))

alpha_t <- ifelse(runif(1) < R, alpha_star, alpha_t)

rgamma(1, shape = n * alpha_t + gamma, rate = sum(y) + lambda)



integrate(function(x) lik_alpha(x, beta_t, alphavec[max_idx]), lower = 0, upper = 20)$value

plot(alphavec , sapply(alphavec, function(x) lik_alpha(x, beta_t, mu_alpha)))

plot(alphavec , sapply(alphavec, function(x) dnorm(x, mu_alpha, sqrt(sigma2_alpha))))


plot(alphavec , sapply(alphavec, function(x) lik_alpha(x, beta_t)))

plot(alphavec , sapply(alphavec, function(x) exp(n*x*log(beta_t) - n * log(gamma(x)) + (x - 1) * sum(log(y)) - 332.2805165)))
max_idx = which.max(sapply(alphavec, function(x) exp(n*x*log(beta_t) - n * log(gamma(x)) + (x - 1) * sum(log(y)) - 332.2805165)))



