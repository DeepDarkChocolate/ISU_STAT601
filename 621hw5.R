getwd()


data = read.csv2("Table6-2.csv", sep = ",")
names(data)[1] <- c("Stratum")

rep()

tab1 <- table(data[1:2])
as.matrix(tab1)
tab1[,1]
w = c(9000 / 20, 9000 / 10, 9000/ 20)
p = 1 / w
w1 = sum(w * tab1[,1]) / sum(tab1[,1])
w2 = sum(w * tab1[,2]) / sum(tab1[,2])
w3 = sum(w * tab1[,3]) / sum(tab1[,3])
w_tilde = c(w1, w2, w3)
pi_tilde = 1 / w_tilde
c(data[2])
data.frame(data)
data[2] <- factor(c(data[2])$Type)

X = with(data,model.matrix(~Type-1))
y = as.numeric(c(data[3])$Y)

Dinv = diag(1 / rep(p, times = c(20, 10, 20)))
Winv = diag(c(1 / (X %*% w_tilde)))
solve((t(X) %*% Dinv %*% Winv %*% X), (t(X) %*% Dinv %*% Winv %*% y))
Winv %*% y
data
w %*% tab1  %*% pi_tilde

muhat = sum(w * c(sum((c(X %*% pi_tilde) * y)[1:20]), sum((c(X %*% pi_tilde) * y)[21:30]),
  sum((c(X %*% pi_tilde) * y)[31:50]))) / 50

e = y - muhat
N_h = rep(9000, 3)
n_h = c(20, 10, 20)
sum(N_h^2 / n_h * (1 - n_h / N_h) * c(var((c(X %*% pi_tilde) * e)[1:20]), var((c(X %*% pi_tilde) * e)[21:30]),
  var((c(X %*% pi_tilde) * e)[31:50]))) / 50^2
