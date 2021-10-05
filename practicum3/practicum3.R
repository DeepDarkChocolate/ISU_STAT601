data = read.csv("greenbeandatfor601.txt", sep = " ")

price = data$x
mvm = data$y
store = data$store

length(table(store))
y_list <- tapply(mvm, store, c)
x_list <- tapply(price, store, c)

library(extRemes)

# loglik_yj = function(beta, gamma, yj, xj){
#   pj = exp(beta[1] + beta[2] * xj) / (1 + exp(beta[1] + beta[2] * xj))
#   lambdaj = exp(gamma[1] + gamma[2] * log(xj))
#   
#   return(sum(log(ifelse(yj == 0, 
#                 (1 - pj) + pj * exp(-lambdaj),
#                 pj * dpois(yj, lambda = lambdaj)))))
# }

loglik_yj = function(beta0j, beta1j, gamma0j, gamma1j, yj, xj){
  pj = exp(beta0j + beta1j * xj) / (1 + exp(beta0j + beta1j * xj))
  lambdaj = exp(gamma0j + gamma1j * log(xj))
  return(sum(log(ifelse(yj == 0, 
                        (1 - pj) + pj * exp(-lambdaj),
                        pj * dpois(yj, lambda = lambdaj)))))
}

loglik_yj_thetaj = function(thetaj, yj, xj){
  beta0j = thetaj[1]
  beta1j = thetaj[2]
  gamma0j = thetaj[3]
  gamma1j = thetaj[4]
  pj = exp(beta0j + beta1j * xj) / (1 + exp(beta0j + beta1j * xj))
  lambdaj = exp(gamma0j + gamma1j * log(xj))
  return(sum(log(ifelse(yj == 0, 
                        (1 - pj) + pj * exp(-lambdaj),
                        pj * dpois(yj, lambda = lambdaj)))))
}


lik_yj_norm = function(beta0j, beta1j, gamma0j, gamma1j, yj, xj, ...,
                       beta0j_ref = 0, beta1j_ref = 0, gamma0j_ref = 0, gamma1j_ref = 0){
  exp(loglik_yj(beta0j, beta1j, gamma0j, gamma1j, yj, xj) -
        loglik_yj(beta0j_ref, beta1j_ref, gamma0j_ref, gamma1j_ref, yj, xj))
}


# lik_yj = function(beta0j, beta1j, gamma0j, gamma1j, yj, xj){
#   pj = exp(beta0j + beta1j * xj) / (1 + exp(beta0j + beta1j * xj))
#   lambdaj = exp(gamma0j + gamma1j * log(xj))
#   
#   return(exp(sum(log(ifelse(yj == 0, 
#                         (1 - pj) + pj * exp(-lambdaj),
#                         pj * dpois(yj, lambda = lambdaj)))) + 14216))
# }

lik_beta0j = function(beta0j, xi0, theta0, beta1j, gamma0j, gamma1j, yj, xj, beta0j_mle){
  devd(beta0j, loc = xi0, scale = theta0) * 
    lik_yj_norm(beta0j, beta1j, gamma0j, gamma1j, yj, xj,
                beta0j_ref = beta0j_mle, beta1j_ref = beta1j, gamma0j_ref = gamma0j, gamma1j_ref = gamma1j)
}

lik_beta1j = function(beta1j, xi1, theta1, beta0j, gamma0j, gamma1j, yj, xj, beta1j_mle){
  devd(beta1j, loc = xi1, scale = theta1) * 
    lik_yj_norm(beta0j, beta1j, gamma0j, gamma1j, yj, xj,
                beta0j_ref = beta0j, beta1j_ref = beta1j_mle, gamma0j_ref = gamma0j, gamma1j_ref = gamma1j)
}

lik_gamma0j = function(gamma0j, mu0, sigma2_0, beta0j, beta1j, gamma1j, yj, xj, gamma0j_mle){
  dnorm(gamma0j, mean = mu0, sd = sqrt(sigma2_0)) *
    lik_yj_norm(beta0j, beta1j, gamma0j, gamma1j, yj, xj,
                beta0j_ref = beta0j, beta1j_ref = beta1j, gamma0j_ref = gamma0j_mle, gamma1j_ref = gamma1j)
}

lik_gamma1j = function(gamma1j, mu1, sigma2_1, beta0j, beta1j, gamma0j, yj, xj, gamma1j_mle){
  dnorm(gamma1j, mean = mu1, sd = sqrt(sigma2_1)) *
    lik_yj_norm(beta0j, beta1j, gamma0j, gamma1j, yj, xj,
                beta0j_ref = beta0j, beta1j_ref = beta1j, gamma0j_ref = gamma0j, gamma1j_ref = gamma1j_mle)
}

#function(xi0, G0, VG0, beta0_sample, )

loglik_yj(7, -9, 2, 0, y_list[[1]], x_list[[1]])
loglik_yj(-2, -7, 2, 0, y_list[[1]], x_list[[1]])
loglik_yj(7, -9, 0, -4, y_list[[1]], x_list[[1]])
loglik_yj(7, -9, -1.38, -5.22, y_list[[1]], x_list[[1]])

loglik_yj(0.4, -3, 0.3, -4.6, y_list[[34]], x_list[[34]])
loglik_yj(0, 0, -3, 4.6, y_list[[34]], x_list[[34]])
loglik_yj(-5.445738, -5.623919,  7.681465,  4.049629, y_list[[34]], x_list[[34]])


glm(ifelse(y_list[[1]] == 0, 0, 1) ~ log(x_list[[1]]), family = binomial)

glm(ifelse(y_list[[1]] == 0, 0, 1) ~ log(x_list[[1]]), family = binomial)
glm(y_list[[1]][y_list[[1]] != 0] ~ log(x_list[[1]])[y_list[[1]] != 0], family = poisson)

optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[1]], xj = x_list[[1]], 
      par = c(7, -9, 0, -4))

optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[1]], xj = x_list[[1]], 
      par = c(0, 0, 0, 0))

glm(ifelse(y_list[[15]] == 0, 0, 1) ~ log(x_list[[15]]), family = binomial)$coefficients
glm(y_list[[15]][y_list[[15]] != 0] ~ log(x_list[[15]])[y_list[[15]] != 0], family = poisson)$coefficients

optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[15]], xj = x_list[[15]], 
      par = c(-4.764304, -11.653234, -1.263197, -6.097698))

glm(ifelse(y_list[[34]] == 0, 0, 1) ~ log(x_list[[34]]), family = binomial)$coefficients
glm(y_list[[34]][y_list[[34]] != 0] ~ log(x_list[[34]])[y_list[[34]] != 0], family = poisson)$coefficients
glm(y_list[[34]] ~ log(x_list[[34]]), family = poisson)$coefficients
optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[34]], xj = x_list[[34]], 
      par = c(0, 0, 5, -1.6666))

index = 83
glm(ifelse(y_list[[index]] == 0, 0, 1) ~ log(x_list[[index]]), family = binomial)$coefficients
glm(y_list[[index]][y_list[[index]] != 0] ~ log(x_list[[index]])[y_list[[index]] != 0], family = poisson)$coefficients
glm(y_list[[index]] ~ log(x_list[[index]]), family = poisson)$coefficients
optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[index]], xj = x_list[[index]], 
      par = c(0, 0, 0, 0))
res = gridSearch(function(...) -loglik_yj_thetaj(...), levels = list(
  seq(from = -100, to = 100, len = 10), seq(from = -100, to = 100, len = 10),
  seq(from = -5, to = 5, len = 10), seq(from = -5, to = 5, len = 10)
), yj = y_list[[index]], xj = x_list[[index]])

theta_mle_mat[,index]
res$minlevels
res$minfun

glm(y_list[[1]][y_list[[1]] == 0] ~ log(x_list[[1]][y_list[[1]] == 0]), family = poisson)

y_list[[15]][y_list[[15]] == 0]

loglik_yj(-5, -11, -1.263, -6.098, y_list[[15]], x_list[[15]])

lik_gamma0j(0.001, 0, 12, 7, -9, -4, y_list[[1]], x_list[[1]], 0)
lik_gamma1j(-3.8, 0, 12, 7, -9, 0, y_list[[1]], x_list[[1]])
lik_gamma1j(-4.2, 0, 12, 6.6165299, -9.3677621, -0.3231684, y_list[[1]], x_list[[1]])

lik_beta0j(0, 0, 2, 0, 0, 0, y_list[[1]], x_list[[1]])
lik_beta0j(-1, 0, 2, 0, 0, 0, y_list[[1]], x_list[[1]])
lik_beta0j(2, 0, 2, 0, 0, 0, y_list[[1]], x_list[[1]])
lik_beta0j(9, 0, 2, 0, 0, 0, y_list[[1]], x_list[[1]])
lik_beta0j(1, 0, 2, 0, 0, 0, y_list[[1]], x_list[[1]])

optimize(function(x) return(-x^2+1), interval = c(-10^10, 10^10), maximum = TRUE)$objective

runif(1, -1, 10^10)
min(4,2,3)

rejection <- function(f, ..., lower = -10, upper = 10){
  ff <- function(x) f(x, ...)
  a_opt = optimize(interval = c(lower, upper), ff, maximum = TRUE)
  if(a_opt$objective == 0){
    print(sprintf("a_opt$maximum = %g", a_opt$maximum))
    print(sprintf("lower = %g", lower))
    print(sprintf("upper = %g", upper))
    print(sprintf("ff(optimized) = %g", ff(optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum)))
    print(sprintf("a_opt$objective = %g", a_opt$objective))
    stop("a == 0")
  }
  cnt = 0
  while(TRUE){
    cnt = cnt + 1
    Y = runif(1, lower, upper)
    U = runif(1)
    if(is.nan(ff(Y) / a_opt$objective)){
      #warning("ff(Y) / a_opt$objective is NaN")
      print(sprintf("U = %g", U))
      print(sprintf("Y = %g", Y))
      print(sprintf("ff(Y) = %g", ff(Y)))
      print(sprintf("a_opt$objective = %g", a_opt$objective))
      stop("ff(Y) / a_opt$objective is NaN")
    }
    if(U <= ff(Y) / a_opt$objective){
      res = Y
      break
    }
    if(cnt > 1000){
      print(sprintf("U = %g", U))
      print(sprintf("Y = %g", Y))
      print(sprintf("ff(Y) = %g", ff(Y)))
      print(sprintf("a_opt$objective = %g", a_opt$objective))
      stop("infinite loop in rejection sampling algorithm")
    }
  }
  return(res)
}  
  
ratiounif <- function(f, ..., lower = -10, upper = 10){
  
  ff <- function(x) f(x, ...)
  ff2 <- function(x) f(x, ...) * x^2

  a_opt = optimize(interval = c(lower, upper), ff, maximum = TRUE)
  a = sqrt(a_opt$objective)
  # if(a_opt$maximum == lower | a_opt$maximum == upper){
  #   warning("a_opt$maximum == lower, upper")
  #   # lower2 = lower
  #   # upper2 = upper
  #   # while(a != 0){
  #   #   lower2 = lower2 / 2
  #   #   upper2 = upper2 / 2
  #   #   a = sqrt(optimize(interval = c(lower2, upper2), ff, maximum = TRUE)$objective)
  #   # }
  #   #stop("a == 0")
  # }
  if(lower < 0){
    b1_opt = optimize(interval = c(lower, min(upper, 0)), ff2, maximum = TRUE)
    b1 = -sqrt(b1_opt$objective)
    # if(b1_opt$maximum == lower | b1_opt$maximum == min(upper, 0)){
    #   warning("b1_opt$maximum == lower,  min(upper, 0))")
    # }
  }else{
    b1 = 0
  }
  if(upper > 0){
    b2_opt = optimize(interval = c(max(lower, 0), upper), ff2, maximum = TRUE)
    b2 = sqrt(b2_opt$objective)
    # if(b2_opt$maximum == max(lower, 0) | b2_opt$maximum == upper){
    #   warning("b2_opt$maximum == max(lower, 0),  upper)")
    # }
  }else{
    b2 = 0
  }
  
  if(is.infinite(a) | is.nan(a)){
    print(sprintf("a_opt$maximum = %g", a_opt$maximum))
    print(sprintf("lower = %g", lower))
    print(sprintf("upper = %g", upper))
    print(sprintf("ff(optimized) = %g", ff(optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum)))
    print(sprintf("a = %g", a))
    print(sprintf("optimize(interval = c(-M, M), ff, maximum = TRUE)$maximum = %g", optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum))
    print(sprintf("b1 = %g", b1))
    print(sprintf("b2 = %g", b2))
    stop("a is infinite or nan")
  }
  
  if(a == 0){
    print(sprintf("a_opt$maximum = %g", a_opt$maximum))
    print(sprintf("lower = %g", lower))
    print(sprintf("upper = %g", upper))
    print(sprintf("ff(optimized) = %g", ff(optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum)))
    print(sprintf("a = %g", a))
    print(sprintf("optimize(interval = c(-M, M), ff, maximum = TRUE)$maximum = %g", optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum))
    print(sprintf("b1 = %g", b1))
    print(sprintf("b2 = %g", b2))
    stop("a == 0")
  }
  
  # print(sprintf("ff(-1) = %g", ff(-1)))
  # print(sprintf("ff(0) = %g", ff(0)))
  # print(sprintf("ff(1) = %g", ff(1)))
  # print(sprintf("ff(optimized) = %g", ff(optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum)))
  # print(sprintf("a = %g", a))
  # print(sprintf("optimize(interval = c(-M, M), ff, maximum = TRUE)$maximum = %g", optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum))
  # print(sqrt(optimize(interval = c(lower, upper), ff, maximum = TRUE)$objective))
  # print(sprintf("b1 = %g", b1))
  # print(sprintf("b2 = %g", b2))
  
  cnt = 0
  while(TRUE){
    cnt = cnt + 1
    u <- runif(1, 0, a)
    v <- runif(1, b1, b2)
    
    if(is.nan(ff(v / u))){
      # warning("ff(v / u) is NaN")
      # print(sprintf("u = %g", u))
      # print(sprintf("v = %g", v))
      # print(sprintf("v / u = %g", v / u))
      next
    }else{
      if(u <= sqrt(ff(v / u))) break 
    }
    
    if(cnt > 1000){
      print(sprintf("a_opt$maximum = %g", a_opt$maximum))
      print(sprintf("lower = %g", lower))
      print(sprintf("upper = %g", upper))
      print(sprintf("ff(optimized) = %g", ff(optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum)))
      print(sprintf("a = %g", a))
      print(sprintf("optimize(interval = c(-M, M), ff, maximum = TRUE)$maximum = %g", optimize(interval = c(lower, upper), ff, maximum = TRUE)$maximum))
      print(sprintf("b1 = %g", b1))
      print(sprintf("b2 = %g", b2))
      
      print(sprintf("u = %g", u))
      print(sprintf("v = %g", v))
      print(sprintf("v / u = %g", v / u))
      stop("infinite loop in ratio of uniform algorithm")
    }
  }
  #print(cnt)
  return(v / u)
}

ratiounif(lik_beta0j, xi0 = 0, theta0 = 2, beta1j = -9, gamma0j = 0, gamma1j = -4,
          yj = y_list[[1]], xj = x_list[[1]], beta0j_mle = 7, lower = 0, upper = 20)

ratiounif(lik_beta1j, xi1 = 0, theta1 = 2, beta0j = 7, gamma0j = 0, gamma1j = -4,
          yj = y_list[[1]], xj = x_list[[1]], beta1j_mle = -9, lower = -20, upper = 0)

ratiounif(lik_gamma0j, mu0 = 0, sigma2_0 = 12, beta0j = 7, beta1j = -9, gamma1j = -4,
          yj = y_list[[1]], xj = x_list[[1]], gamma0j_mle = 0,lower = -1, upper = 1)

ratiounif(lik_gamma1j, mu1 = 0, sigma2_1 = 12, beta0j = 7, beta1j = -9, gamma0j = 0,
          yj = y_list[[1]], xj = x_list[[1]], gamma1j_mle = -4, lower = -5, upper = -3)

######################

theta_mle_mat[,15]

ratiounif(lik_beta0j, xi0 = 0, theta0 = 2, beta1j = -15, gamma0j = -1, gamma1j = -5,
          yj = y_list[[1]], xj = x_list[[1]], beta0j_mle = 10, lower = 0, upper = 20)

ratiounif(lik_beta1j, xi1 = 0, theta1 = 2, beta0j = 7, gamma0j = 0, gamma1j = -4,
          yj = y_list[[1]], xj = x_list[[1]], beta1j_mle = -9, lower = -20, upper = 0)

ratiounif(lik_gamma0j, mu0 = 0, sigma2_0 = 12, beta0j = 1.774711, beta1j = 2.903434, gamma1j = -2.675319,
          yj = y_list[[13]], xj = x_list[[13]], gamma0j_mle = 1.732733,lower = 1.732733+0.5, upper = 1.732733 - 0.5)

lower_tmp = 1.732733
while(lik_gamma0j(lower_tmp, 0, 12, 1.774711, 2.903434, -2.675319, y_list[[13]], x_list[[13]], 1.732733) > 0) lower_tmp = lower_tmp - 0.1
upper_tmp = 1.732733
while(lik_gamma0j(upper_tmp, 0, 12, 1.774711, 2.903434, -2.675319, y_list[[13]], x_list[[13]], 1.732733) > 0) upper_tmp = upper_tmp + 0.1


ratiounif(lik_gamma1j, mu1 = 0, sigma2_1 = 2, beta0j = 14.803634, beta1j = -22.745317, gamma0j = -0.7878141,
          yj = y_list[[1]], xj = x_list[[1]], gamma1j_mle = -4.6019937, lower = -5, upper = -3)



lik_yj(0, 0, 0, 0, y_list[[1]], x_list[[1]])
loglik_yj(0, 0, 0, 0, y_list[[1]], x_list[[1]])
loglik_yj(1, 0, 0, 0, y_list[[1]], x_list[[1]])

lik_yj_norm(2, 0, 0, 0, y_list[[1]], x_list[[1]])

library(NMOF)
res <- gridSearch(function(...) -loglik_yj_thetaj(...), levels = list(7, -9,
                                          seq(from = -10, to = 10, len = 20),
                                          seq(from = -10, to = 10, len = 20)),
           yj = y_list[[1]], xj = x_list[[1]])

res$minfun
res$minlevels
theta_mle_mat = NULL
for(j in 1:length(y_list)){
  beta_initial = glm(ifelse(y_list[[j]] == 0, 0, 1) ~ log(x_list[[j]]), family = binomial)$coefficients
  gamma_initial = glm(y_list[[j]][y_list[[j]] != 0] ~ log(x_list[[j]])[y_list[[j]] != 0], family = poisson)$coefficients
  if(is.infinite(loglik_yj(beta_initial[1], beta_initial[2], 
                           gamma_initial[1], gamma_initial[2], 
                           y_list[[j]], x_list[[j]]))){
    print(j)
    for(i2 in seq(from = -10, to = 10, len = 20)){
      for (j2 in seq(from = -10, to = 10, len = 20)){
        if(is.finite(loglik_yj(beta_initial[1], beta_initial[2], i2, j2, y_list[[j]], x_list[[j]]))){
          gamma_initial[1] = i2
          gamma_initial[2] = j2
        }
      }
    }
  }
  thetaj_mle = optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[j]], xj = x_list[[j]], 
                     par = c(beta_initial[1], beta_initial[2],gamma_initial[1], gamma_initial[2]))$par
  theta_mle_mat = cbind(theta_mle_mat, thetaj_mle)
}

theta_mle_mat = NULL
for(j in 1:length(y_list)){
  beta_initial = glm(ifelse(y_list[[j]] == 0, 0, 1) ~ log(x_list[[j]]), family = binomial)$coefficients
  gamma_initial = glm(y_list[[j]][y_list[[j]] != 0] ~ log(x_list[[j]])[y_list[[j]] != 0], family = poisson)$coefficients
  if(is.infinite(loglik_yj(beta_initial[1], beta_initial[2], 
                           gamma_initial[1], gamma_initial[2], 
                           y_list[[j]], x_list[[j]]))){
    print(j)
    for(i2 in seq(from = -10, to = 10, len = 20)){
      for (j2 in seq(from = -10, to = 10, len = 20)){
        if(is.finite(loglik_yj(0, 0, i2, j2, y_list[[j]], x_list[[j]]))){
          gamma_initial[1] = i2
          gamma_initial[2] = j2
          break
        }
      }
    }
  }
  thetaj_mle = optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[j]], xj = x_list[[j]], 
                     par = c(0, 0,gamma_initial[1], gamma_initial[2]))$par
  theta_mle_mat = cbind(theta_mle_mat, thetaj_mle)
}
row.names(theta_mle_mat) <- NULL
#theta_mle_mat[,104]
theta_mle_mat[1:2,83] <- c(0, 0)



Gibbs <- function(y_list, x_list, theta_mle_mat, MCsize){
  G0 = G1 = B0 = B1 = 0
  VG0 = VG1 = VB0 = VB1 = 25
  sqrtVG0 = sqrtVG1 = sqrtVB0 = sqrtVB1 = sqrt(VG0)
  A0 = A1 = 8
  C0 = C1 = 25
  S = length(y_list)
  
  xi_t = c(G0, G1)
  theta_t = c(A0 / 2, A1 / 2)
  mu_t = c(B0, B1)
  sigma2_t = c(C0 / 2, C1 / 2)
  
  beta0_sample = vector(mode = "numeric", length = S)
  beta1_sample = vector(mode = "numeric", length = S)
  gamma0_sample = vector(mode = "numeric", length = S)
  gamma1_sample = vector(mode = "numeric", length = S)
  para_sample = vector(mode = "numeric", length = 8)
  
  beta0_res = NULL
  beta1_res = NULL
  gamma0_res = NULL
  gamma1_res = NULL
  para_res = NULL
  
  # theta_mle_mat = NULL
  # for(j in 1:S){
  #   beta_initial = glm(ifelse(y_list[[j]] == 0, 0, 1) ~ log(x_list[[j]]), family = binomial)$coefficients
  #   gamma_initial = glm(y_list[[j]][y_list[[j]] != 0] ~ log(x_list[[j]])[y_list[[j]] != 0], family = poisson)$coefficients
  #   if(is.infinite(loglik_yj(beta_initial[1], beta_initial[2], 
  #                            gamma_initial[1], gamma_initial[2], 
  #                            y_list[[j]], x_list[[j]]))){
  #     print(j)
  #     for(i2 in seq(from = -10, to = 10, len = 20)){
  #       for (j2 in seq(from = -10, to = 10, len = 20)){
  #         if(is.finite(loglik_yj(beta_initial[1], beta_initial[2], i2, j2, y_list[[j]], x_list[[j]]))){
  #           gamma_initial[1] = i2
  #           gamma_initial[2] = j2
  #         }
  #       }
  #     }
  #   }
  #   
  #   
  #   thetaj_mle = optim(function(...) -loglik_yj_thetaj(...), yj = y_list[[j]], xj = x_list[[j]], 
  #                      par = c(beta_initial[1], beta_initial[2],gamma_initial[1], gamma_initial[2]))$par
  #   theta_mle_mat = cbind(theta_mle_mat, thetaj_mle)
  # }
  
  for(MCnum in 1:MCsize){
    for(j in 1:S){
      print(j)
      thetaj_mle = theta_mle_mat[,j]
      if(MCnum == 1){
        beta1j = thetaj_mle[2]
        gamma0j = thetaj_mle[3]
        gamma1j = thetaj_mle[4]
      }else{
        beta1j = beta1_sample[j]
        gamma0j = gamma0_sample[j]
        gamma1j = gamma1_sample[j]
      }
      print("beta1j"); print(beta1j)
      print("gamma0j"); print(gamma0j)
      print("gamma1j"); print(gamma1j)
      print("beta0_sample[j]")
      if(is.nan(lik_beta0j(thetaj_mle[1], xi_t[1], theta_t[1], beta1j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[1]))){
        beta0_sample[j] = beta0j = thetaj_mle[1]
      }else{
        lower = thetaj_mle[1]
        gap = 0.1
        while(TRUE){
          gap = gap * 2
          tmpval = lik_beta0j(lower - gap, xi_t[1], theta_t[1], beta1j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[1])
          if(is.nan(tmpval)) break
          else if(tmpval <= 1.0e-10){
            lower = lower - gap
            break
          } 
        }
        
        upper = thetaj_mle[1]
        gap = 0.1
        while(TRUE){
          gap = gap * 2
          tmpval = lik_beta0j(upper + gap, xi_t[1], theta_t[1], beta1j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[1])
          if(is.nan(tmpval)) break
          else if(tmpval <= 1.0e-10){
            upper = upper + gap
            break
          } 
        }

        #while(lik_beta0j(upper, xi_t[1], theta_t[1], beta1j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[1]) > 1.0e-10) upper = upper + gap / 4
        beta0_sample[j] = beta0j = ratiounif(lik_beta0j, xi0 = xi_t[1], theta0 = theta_t[1], 
                                             beta1j = beta1j, gamma0j = gamma0j, gamma1j = gamma1j,
                                             yj = y_list[[j]], xj = x_list[[j]], beta0j_mle = thetaj_mle[1], 
                                             #lower = min(thetaj_mle[1] - 10, 0), upper = max(thetaj_mle[1] + 10, 0))
                                             lower = lower, upper = upper)
      }
      print("beta1_sample[j]")
      print("beta0j"); print(beta0j)
      
      if(is.nan(lik_beta1j(thetaj_mle[2], xi_t[2], theta_t[2], beta0j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[2]))){
        beta1_sample[j] = beta1j = thetaj_mle[2]
      }else{
        lower = thetaj_mle[2]
        gap = 0.1
        while(TRUE){
          gap = gap * 2
          tmpval = lik_beta1j(lower - gap, xi_t[2], theta_t[2], beta0j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[2])
          if(is.nan(tmpval)) break
          else if(tmpval <= 1.0e-10){
            lower = lower - gap
            break
          } 
        }

        upper = thetaj_mle[2]
        gap = 0.1
        while(TRUE){
          gap = gap * 2
          tmpval = lik_beta1j(upper + gap, xi_t[2], theta_t[2], beta0j, gamma0j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[2])
          if(is.nan(tmpval)) break
          else if(tmpval <= 1.0e-10){
            upper = upper + gap
            break
          } 
        }
        
        beta1_sample[j] = beta1j = ratiounif(lik_beta1j, xi1 = xi_t[2], theta1 = theta_t[2], 
                                             beta0j = beta0j, gamma0j = gamma0j, gamma1j = gamma1j,
                                             yj = y_list[[j]], xj = x_list[[j]], beta1j_mle = thetaj_mle[2], 
                                             #lower = min(thetaj_mle[2] - 10, 0), upper = max(thetaj_mle[2] + 10, 0))
                                             lower = lower, upper = upper)
      }
      print("gamma0_sample[j]")
      if(is.nan(lik_gamma0j(thetaj_mle[3], mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))){
        gamma0_sample[j] = gamma0j = thetaj_mle[3]
      }else{
        lower = thetaj_mle[3]
        gap = 0.1
        if(is.nan(lik_gamma0j(lower - gap, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))){
          lower = thetaj_mle[3]
        }else{
          while(lik_gamma0j(lower - gap, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]) <= 1.0e-10) gap = gap / 2
          while(lik_gamma0j(lower, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]) > 1.0e-10) lower = lower - gap / 4
        }
        upper = thetaj_mle[3]
        gap = 0.1
        if(is.nan(lik_gamma0j(upper + gap, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))){
          upper = thetaj_mle[3]
        }else{
          #print(lik_gamma0j(upper, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.25, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.5, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.75, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.1, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          while(lik_gamma0j(upper + gap, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]) <= 1.0e-10) gap = gap / 2
          while(lik_gamma0j(upper, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]) > 1.0e-10) upper = upper + gap / 4
          
        }
        # print("lower"); print(lower)
        # print("thetaj_mle[3]"); print(thetaj_mle[3])
        # print("upper"); print(upper)
        gamma0_sample[j] = 
          gamma0j = 
          rejection(lik_gamma0j, mu0 = mu_t[1], sigma2_0 = sigma2_t[1],
                    beta0j = beta0j, beta1j = beta1j, gamma1j = gamma1j,
                    yj = y_list[[j]], xj = x_list[[j]], gamma0j_mle = thetaj_mle[3],
                    lower = lower, upper = upper)
      }
      # print("max(log(x_list[[j]]) * gamma1j)")
      # print(max(log(x_list[[j]]) * gamma1j))
      # print("min(log(x_list[[j]]) * gamma1j)")
      # print(min(log(x_list[[j]]) * gamma1j))
      # seqvec <- seq(from = thetaj_mle[3] - .1, thetaj_mle[3] + .1, len = 200)
      # tmpvec <- sapply(seqvec,
      #                  function(x) loglik_yj(beta0j, beta1j, x, gamma1j, y_list[[j]], x_list[[j]]))
      # theta3_ref = seqvec[which.max(tmpvec)]
      # theta3_ref = thetaj_mle[3] + 10^(-5)
      # lower = theta3_ref
      # gap = 0.1
      # print("beta0j"); print(beta0j)
      # print("beta1j"); print(beta1j)
      # print("gamma0j"); print(gamma0j)
      # print("gamma1j"); print(gamma1j)
      # print("thetaj_mle"); print(thetaj_mle)
      # 
      # if(is.infinite(max(tmpvec))) warning("is.infinite(max(tmpvec))")
      # print("which.max(tmpvec)"); print(which.max(tmpvec))
      # print("seqvec[which.max(tmpvec)]"); print(seqvec[which.max(tmpvec)])
      # print("loglik_yj(beta0j, beta1j, theta3_ref, gamma1j, y_list[[j]], x_list[[j]])")
      # print(loglik_yj(beta0j, beta1j, theta3_ref, gamma1j, y_list[[j]], x_list[[j]]))
      # print(tmpvec)

      # while(lik_gamma0j(lower - gap, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], theta3_ref) <= 1.0e-10) gap = gap / 2
      # print("gap"); print(gap)
      # while(lik_gamma0j(lower, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], theta3_ref) > 1.0e-10) lower = lower - gap / 4
      # upper = theta3_ref
      # gap = 0.1
      # while(lik_gamma0j(upper + gap, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], theta3_ref) <= 1.0e-10) gap = gap / 2
      # while(lik_gamma0j(upper, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], theta3_ref) > 1.0e-10) upper = upper + gap / 4
      # print("lower"); print(lower)
      # print("upper"); print(upper)
      # print("thetaj_mle[3]"); print(thetaj_mle[3])
      # gamma0_sample[j] = gamma0j = ratiounif(lik_gamma0j, mu0 = mu_t[1], sigma2_0 = sigma2_t[1], 
      #                    beta0j = beta0j, beta1j = beta1j, gamma1j = gamma1j,
      #                    yj = y_list[[j]], xj = x_list[[j]], gamma0j_mle = thetaj_mle[3], 
      #                    lower = lower, upper = upper)
      
      # print("gamma0j"); print(gamma0j)
      print("gamma1_sample[j]")
      
      if(is.nan(lik_gamma1j(thetaj_mle[4], mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]))){
        gamma1_sample[j] = gamma1j = thetaj_mle[4]
      }else{
        lower = thetaj_mle[4]
        gap = 0.1
        if(is.nan(lik_gamma1j(lower - gap, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]))){
          lower = thetaj_mle[4]
        }else{
          while(lik_gamma1j(lower - gap, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) <= 1.0e-10) gap = gap / 2
          while(lik_gamma1j(lower, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) > 1.0e-10) lower = lower - gap / 4
        }
        upper = thetaj_mle[4]
        gap = 0.1
        if(is.nan(lik_gamma1j(upper + gap, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]))){
          upper = thetaj_mle[4]
        }else{
          #print(lik_gamma0j(upper, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.25, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.5, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.75, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          # print(lik_gamma0j(upper + 0.1, mu_t[1], sigma2_t[1], beta0j, beta1j, gamma1j, y_list[[j]], x_list[[j]], thetaj_mle[3]))
          while(lik_gamma1j(upper + gap, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) <= 1.0e-10) gap = gap / 2
          while(lik_gamma1j(upper, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) > 1.0e-10) upper = upper + gap / 4
          
        }
        # print("lower"); print(lower)
        # print("thetaj_mle[3]"); print(thetaj_mle[3])
        # print("upper"); print(upper)
        gamma1_sample[j] = 
          gamma1j = 
          rejection(lik_gamma1j, mu1 = mu_t[2], sigma2_1 = sigma2_t[2],
                    beta0j = beta0j, beta1j = beta1j, gamma0j = gamma0j,
                    yj = y_list[[j]], xj = x_list[[j]], gamma1j_mle = thetaj_mle[4],
                    lower = lower, upper = upper)
      }
      
      
      
      # lower = thetaj_mle[4]
      # gap = 0.1
      # while(lik_gamma1j(lower - gap, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) <= 1.0e-10) gap = gap / 2
      # while(lik_gamma1j(lower, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) > 1.0e-10) lower = lower - gap / 4
      # upper = thetaj_mle[4]
      # gap = 0.1
      # while(lik_gamma1j(upper + gap, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) <= 1.0e-10) gap = gap / 2
      # while(lik_gamma1j(upper, mu_t[2], sigma2_t[2], beta0j, beta1j, gamma0j, y_list[[j]], x_list[[j]], thetaj_mle[4]) > 1.0e-10) upper = upper + gap / 4
      # #print("lower"); print(lower)
      # #print("upper"); print(upper)
      # gamma1_sample[j] = gamma1j = rejection(lik_gamma1j, mu1 = mu_t[2], sigma2_1 = sigma2_t[2], 
      #                     beta0j = beta0j, beta1j = beta1j, gamma0j = gamma0j,
      #                     yj = y_list[[j]], xj = x_list[[j]], gamma1j_mle = thetaj_mle[4], 
      #                     lower = lower, upper = upper)
      #print("thetaj_mle"); print(thetaj_mle)
    }
    beta0_res = cbind(beta0_res, beta0_sample)
    beta1_res = cbind(beta1_res, beta1_sample)
    gamma0_res = cbind(gamma0_res, gamma0_sample)
    gamma1_res = cbind(gamma1_res, gamma1_sample)
  }
  
  return(list(beta0_res, beta1_res, gamma0_res, gamma1_res))
}

set.seed(1)
Gibbs_res = Gibbs(y_list, x_list, theta_mle_mat, 3)


theta_mle_mat[1:2,72] <- c(0, 0)
theta_mle_mat[1:2,83] <- c(0, 0)
theta_mle_mat[1:2,31] <- c(0, 0)
