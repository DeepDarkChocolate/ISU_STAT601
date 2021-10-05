data = read.csv("greenbeandat.txt", sep = " ")

price = data$price
mvm = data$mvm
store = data$store

hist(mvm)
summary(mvm)
sum(mvm == 0) / length(mvm)

idx = c(1202, 1400, 1866, 1026)
price[store == 1400][1:100]
plot(mvm[store == 1202])
plot(mvm[store == 1400])
plot(mvm[store == 1866])
plot(mvm[store == 1026])

plot(mvm[store == 1026] ~ price[store == 1026])
boxplot(mvm[store == 1026] ~ price[store == 1026])
summary(mvm[store == 1202])

acf(mvm[store == 1400])
acf(price[store == 1400])
acf(mvm[store == 1202])
acf(mvm[store == 1866])
acf(mvm[store == 1026])

pacf(mvm[store == 1400])
pacf(price[store == 1400])
pacf(mvm[store == 1202])
pacf(mvm[store == 1866])
pacf(mvm[store == 1026])

library(forecast)

#for(i in idx) print(auto.arima(price[store == i]))
auto.arima(mvm[store == 1202])
auto.arima(mvm[store == 1866])
auto.arima(mvm[store == 1026])

plot(price[store == 1400][2:1000] - price[store == 1400][1:999])

library(tseries)
arma(price[store == 1400], order = c(8, 0))
arima(price[store == 1400], order = c(0, 1, 0))
arima(price[store == 1866], order = c(0, 1, 0))
arima(price[store == 1026], order = c(0, 1, 0))
arima(price[store == 1026], order = c(0, 1, 0), xreg = mvm[store == 1026])
arima(price[store == 1400], order = c(8,0, 0))

negloglik = function(theta, x, y){
  p = theta[1]
  beta0 = theta[2]
  beta1 = theta[3]
  lambda = exp(beta0 + x * beta1)
  if(p < 0 | p > 1) return(-Inf)
  return(-sum(ifelse(y == 0, log((1 - p) + p * exp(-lambda)), log(p) + dpois(y, lambda, log = TRUE))))
}

mean(mvm)

glm(mvm ~ price, family = poisson())
res = optim(c(0.20, 6, -8),negloglik, x = price[store == 1400], y = mvm[store == 1400])

optim(c(0.20, 6, -8),negloglik, x = price[store == 1400], y = mvm[store == 1400])

index = 1202
res = optim(c(0.20, 6, -8),negloglik, x = price[store == index], y = mvm[store == index])
plot(mvm[store == index] ~ price[store == index], main = index)
p = res$par[1]
beta0 = res$par[2]
beta1 = res$par[3]
lambda = exp(beta0 + price[store == index] * beta1)
Ey = p *(lambda - exp(-lambda))
points(Ey ~ price[store == index], col = "red")
legend("topright", c("oberved", "conditional expectation"), col = c("black", "red"), pch = 1)
res$par

plot(Ey - mvm[store == index])
acf(Ey - mvm[store == index])

index = 1400
res = optim(c(0.20, 6, -8),negloglik, x = price[store == index], y = mvm[store == index])
plot(mvm[store == index] ~ price[store == index], main = index)
p = res$par[1]
beta0 = res$par[2]
beta1 = res$par[3]
lambda = exp(beta0 + price[store == index] * beta1)
Ey = p *(lambda - exp(-lambda))
points(Ey ~ price[store == index], col = "red")
legend("topright", c("oberved", "conditional expectation"), col = c("black", "red"), pch = 1)
res$par

index = 1866
res = optim(c(0.20, 6, -8),negloglik, x = price[store == index], y = mvm[store == index])
plot(mvm[store == index] ~ price[store == index], main = index)
p = res$par[1]
beta0 = res$par[2]
beta1 = res$par[3]
lambda = exp(beta0 + price[store == index] * beta1)
Ey = p *(lambda - exp(-lambda))
points(Ey ~ price[store == index], col = "red")
legend("topright", c("oberved", "conditional expectation"), col = c("black", "red"), pch = 1)
res$par

index = 1026
res = optim(c(0.20, 6, -8),negloglik, x = price[store == index], y = mvm[store == index])
plot(mvm[store == index] ~ price[store == index], main = index)
p = res$par[1]
beta0 = res$par[2]
beta1 = res$par[3]
lambda = exp(beta0 + price[store == index] * beta1)
Ey = p *(lambda - exp(-lambda))
points(Ey ~ price[store == index], col = "red")
legend("topright", c("oberved", "conditional expectation"), col = c("black", "red"), pch = 1)
res$par
