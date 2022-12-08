## -----------------------------------------------------------------------------
library(StatComp22067)

## -----------------------------------------------------------------------------
#散点图
n = 50
x = runif(n, min=-2, max=2)
y = x^2 + rnorm(n)
plot(x, y, type="p")

## -----------------------------------------------------------------------------
#用表格展示R中mtcars数据集
library(knitr)
kable(head(mtcars),align="c")

## -----------------------------------------------------------------------------
#取出数据框的列名即为字符串
colnames(mtcars)

## -----------------------------------------------------------------------------
set.seed(123)
a = 2
b = 2
n = 1000
u = runif(n)
x = b/(1-u)^(1/a)
hist(x, probability=TRUE, main="f(x)")
y = seq(0, 80, 0.5)
lines(y, a*b^a/y^(a+1))

## -----------------------------------------------------------------------------
gamma <- function(x){
  return(factorial(x-1))
}

rbeta_ar <- function(n, a, b){
  k = 0
  y = numeric(n) 
  #density function of Beta(a,b)
  fbeta <- function(x){
    return(gamma(a+b)/gamma(a)/gamma(b)*x^(a-1)*(1-x)^(b-1))
  }
  x0 = as.numeric(optimize(fbeta,lower = 0, upper = 1, maximum = TRUE)[[1]])
  while(k < n){
    u = runif(1)
    x = runif(1)
    if((fbeta(x)/(fbeta(x0)+0.5)) > u){
      k = k+1
      y[k] = x
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
set.seed(234)
n = 1000
x = rbeta_ar(n, 3, 2)
y = seq(0, 1, 0.01)
hist(x, probability = TRUE, main = "Beta(3,2)", ylim = c(0,2))
lines(y, dbeta(y, 3, 2))

## -----------------------------------------------------------------------------
set.seed(345)
n = 1000
r = 4 
beta = 2
lambda = rgamma(n, r, beta)
x = rexp(n, lambda)
hist(x, probability = TRUE, main="Exponential-Gamma mixture")

## -----------------------------------------------------------------------------
set.seed(456)
n = 1000
r = 4 
beta = 2
lambda = rgamma(n, r, beta)
x = rexp(n, lambda)
y = seq(0, 10, 0.1)
hist(x, probability = TRUE, main="Exponential-Gamma mixture")
lines(y, r*beta^r/(beta+y)^(r+1))

## -----------------------------------------------------------------------------
# define a quicksort funciton
quicksort <- function(x){
  num = length(x)
  if(num < 2){return(x)}
  else{
    z = x[1] #In Step 1, pick x[1] to compare
    y = x[-1]
    lower = y[y<z]
    upper = y[y>=z]
    return(c(quicksort(lower),z,quicksort(upper)))
    }
}
# define a function to compute time over 100 simulations
avgtime <- function(n){
  t = numeric(100)
  for(i in 1:100){
    x = sample(1:n)
    t[i] = system.time(quicksort(x))[1]
  }
  return(mean(t))
}

## -----------------------------------------------------------------------------
set.seed(123)
a = numeric(5)
t = numeric(5)
n = c(1e4, 2*1e4, 4*1e4, 6*1e4, 8*1e4)
for (i in 1:5) {
  a[i] = avgtime(n[i])
  t[i] = n[i]*log(n[i])
}
fit = lm(a~t)
summary(fit)
plot(t,a,pch=16,xlab="nlog(n)",ylab="time")
abline(fit$coef[[1]],fit$coef[[2]],col="cyan")

## -----------------------------------------------------------------------------
MC_anti <- function(n, anti=FALSE){
  u = runif(n/2)
  if(anti) {
    v = 1-u
    return((mean(exp(u))+mean(exp(v)))/2)
    }
  else {
    v = runif(n/2)
    u = c(u,v)
    return(mean(exp(u)))
    }
}
set.seed(0)
m = 1e4
v1 <- v2 <- numeric(1000)
for(i in 1:1000){
 v1[i] = MC_anti(m)
 v2[i] = MC_anti(m,TRUE)
}
var1 = var(v1)
var2 = var(v2)
c(var1, var2, (var1-var2)/var1)

## -----------------------------------------------------------------------------
u = c(1,2,3)
v = 1-u
exp(c(u,v))
mean(exp(c(u,v)))

## -----------------------------------------------------------------------------
g <- function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>=1)
}
x = seq(1,8,0.01)
gs <- c(expression(g(x)),expression(f1(x)),expression(f2(x)))
par(mfrow=c(1,2))
# figure of g, f1, f2
plot(x, g(x), type="l", ylab="", ylim=c(0,0.6), lwd = 2, col=1)
lines(x, dnorm(x), lwd=2, col=2)
lines(x, dgamma(x,3,2), lwd=2, col=3)
legend("topright", legend = gs, lty=1, lwd=2, inset = 0.02,col=1:3)

# figure of g/f1, g/f2
plot(x, g(x)/dnorm(x), type="l", ylab="", ylim=c(0,5), lwd = 2, col=2)
lines(x, g(x)/dgamma(x,3,2), lwd=2, col=3)
legend("topright", legend = gs[-1], lty=1, lwd=2, inset = 0.02,col=2:3)

## -----------------------------------------------------------------------------
set.seed(0)
m = 1e4
theta <- se <- numeric(2)
# using f1
x <- rnorm(m) 
fg <- g(x) / dnorm(x)
theta[1] <- mean(fg)
se[1] <- sd(fg)
# using f2
x <- rgamma(m,3,2) 
fg <- g(x) / dgamma(x,3,2)
theta[2] <- mean(fg)
se[2] <- sd(fg)
rbind(theta, se)

## -----------------------------------------------------------------------------
se^2

## -----------------------------------------------------------------------------
set.seed(0)
m <- 1e6
est1 <- sd1 <- 0
g <- function(x){
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u <- runif(m) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
est1 <- mean(fg)
sd1 <- sd(fg)

## -----------------------------------------------------------------------------
M <- 10000
k <- 5 
m <- M/k # replicates per stratum
T1 <- T2 <- numeric(k)
est2 <- numeric(2)
fj <- function(x)5*exp(-x)/(1-exp(-1)) # fj(x)
g <- function(x)exp(-x)/(1+x^2)

## -----------------------------------------------------------------------------
set.seed(0)
for(j in 1:k){
  u = runif(m)
  x = -log(1-(1-exp(-1))*(u+j-1)/5) # inverse transform method
  T1[j] = mean(g(x)/fj(x))
  T2[j] = var(g(x)/fj(x))
}
est2[1] = sum(T1) 
est2[2] = sum(T2)
round(c(est1,est2[1]),4)
round(c(sd1,sqrt(est2[2])),5)

## -----------------------------------------------------------------------------
set.seed(0)
m = 1000 # number of replicates
n = 20 # sample size
alpha = 0.05
t = qt(1-alpha/2, n-1)
lci = uci = numeric(m)
ecp = 0
for(i in 1:m){
  x = rnorm(n)
  xbar = mean(x)
  s = sd(x)
  uci[i] = xbar+s*t/sqrt(n)
  lci[i] = xbar-s*t/sqrt(n)
  if(lci[i]<=0 & 0<=uci[i]){ecp = ecp+1}
}
ecp/m

## -----------------------------------------------------------------------------
count5test <-  function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

## -----------------------------------------------------------------------------
set.seed(1)
sigma1 = 1
sigma2 = 1.5
alpha = 0.055
m = 10000
n = c(10,100,1000)
power5 = powerF = numeric(3)
repl1 = repl2 = numeric(m) 
for(i in 1:3){
  for(j in 1:m){
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    repl1[j] = count5test(x, y)
    repl2[j] = var.test(x,y,conf.level=alpha)$p.value <= alpha
  }
  power5[i] = mean(repl1)
  powerF[i] = mean(repl2)
}
rbind(power5,powerF)

## -----------------------------------------------------------------------------
set.seed(0)
data("aircondit",package = "boot")
n = 12
B = 1000
x = aircondit$hours
lambda = n/sum(x)
lambdastar = numeric(B)
for(i in 1:B){
  xstar = sample(x,replace=TRUE)
  lambdastar[i] = n/sum(xstar)
}

## -----------------------------------------------------------------------------
round(c(MLE=lambda,bias=mean(lambdastar)-lambda,se=sd(lambdastar)),4)

## -----------------------------------------------------------------------------
library(boot)
set.seed(1)
boot.mean <- function(x,i) mean(x[i])
x = aircondit$hours
meantime <- boot(data=x,statistic=boot.mean,R = 1e3)
ci <- boot.ci(meantime,type=c("norm","basic","perc","bca"))

## -----------------------------------------------------------------------------
library(dplyr)
table1 = data.frame(Methods=c("basic","normal","percentile","BCa"),lower=round(c(ci$basic[4],ci$normal[2],ci$percent[4],ci$bca[4]),2),upper=round(c(ci$basic[5],ci$normal[3],ci$percent[5],ci$bca[5]),2))
knitr::kable(table1, format = "html", align=rep('c',3)) %>% 
  kableExtra::kable_styling()

## -----------------------------------------------------------------------------
len_norm = ci$normal[3]-ci$normal[2]
len_basic = ci$basic[5]-ci$basic[4]
len_perc = ci$percent[5]-ci$percent[4]
len_bca = ci$bca[5]-ci$bca[4]
cat('norm =',len_norm,'basic =',len_basic,'perc =',len_perc,'BCa =',len_bca)

## -----------------------------------------------------------------------------
theta = mean(x)
thetastar = meantime$t
hist(thetastar);abline(v=theta,col='red',lwd=2)
thetasort = sort(thetastar)
abline(v=thetasort[25],col="cyan",lwd=2);abline(v=thetasort[975],col="green",lwd=2)

## -----------------------------------------------------------------------------
mu = 1
sigma = 1
n = 30 # sample size
m = 1000 # times to replicate
set.seed(0)
ci.norm <- ci.basic <- ci.perc <- matrix(NA,m,2)
for(i in 1:m){
  normdata = rnorm(n,mu,sigma)
  de <- boot(data=normdata, statistic=boot.mean, R = 1000)
  ci <- boot.ci(de,type=c("norm","basic","perc"))
  ci.norm[i,] <- ci$norm[2:3]
  ci.basic[i,] <- ci$basic[4:5]
  ci.perc[i,] <- ci$percent[4:5]
}

## -----------------------------------------------------------------------------
miss_left = round(c(mean(ci.norm[,2]<=mu),mean(ci.basic[,2]<=mu),mean(ci.perc[,2]<=mu)),3)
cover_rate = c(mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu))
miss_right = c(mean(ci.norm[,1]>=mu),mean(ci.basic[,1]>=mu),mean(ci.perc[,1]>=mu))

table2 = data.frame(Methods=c("normal","basic","percentile"),miss_left, cover_rate, miss_right)
knitr::kable(table2, format = "html", align=rep('c',3)) %>%
  kableExtra::kable_styling()

## -----------------------------------------------------------------------------
data("scor",package="bootstrap")
data1 = scor
n = nrow(data1)
sigma.hat = matrix(0,5,5)
lambda.hat = numeric(5)
theta.jack = numeric(n)
for(i in 1:n){
  data_jack = data1[-i,]  # leave-one-out
  sigma.hat = (n-2)*cov(data_jack)/(n-1) # MLE of Sigma
  lambda.hat = eigen(sigma.hat)$values
  theta.jack[i] = lambda.hat[1]/sum(lambda.hat)
}

## -----------------------------------------------------------------------------
sigma = (n-1)*cov(data1)/n
lambda.hat = eigen(sigma)$values
theta.hat = lambda.hat[1]/sum(lambda.hat)
bias.jack = (n-1)*(mean(theta.jack)-theta.hat)
se.jack = sqrt((n-1)*mean((theta.jack-mean(theta.hat))^2))
round(c(bias.jack=bias.jack,se.jack=se.jack),3)

## -----------------------------------------------------------------------------
rm(list=ls()) #clear memory

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic) 
e1 <- e2 <- e3 <- e4 <- matrix(0,n,n) # store the squared prediction errors

# fit models on leave-two-out samples
for (i in 1:(n-1)){
  for (j in (i+1):n){
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2]*chemical[c(i,j)]
    e1[i,j] <- mean((magnetic[c(i,j)] - yhat1)^2)
  
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1]+J2$coef[2]*chemical[c(i,j)]+J2$coef[3]*chemical[c(i,j)]^2
    e2[i,j] <- mean((magnetic[c(i,j)] - yhat2)^2)
  
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e3[i,j] <- mean((magnetic[c(i,j)] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e4[i,j] <- mean((magnetic[c(i,j)] - yhat4)^2)
  }
}

## -----------------------------------------------------------------------------
2*n/(n-1)*c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
library(MASS)
n = 20
mu = c(0,0)
Sigma = matrix(c(1,0.2,0.2,1),2,2)
set.seed(0)
X = mvrnorm(n, mu, Sigma)
x = X[,1]; y = X[,2] # two samples
R = 1000
reps <- numeric(R)
r0 <- cor(x, y, method="spearman")
for(i in 1:R){
  x_per <- sample(x, replace=FALSE)
  reps[i] = cor(x_per, y, method="spearman")
}
p = mean(abs(reps)>=abs(r0)) # empirical p-value
round(c(p,cor.test(x,y,method="spearman")$p.value),4)

## -----------------------------------------------------------------------------
rm(list=ls())

## -----------------------------------------------------------------------------
f <- function(x){
  return(exp(-abs(x))/2)
}
invF <- function(x){
  return(-log(2*(1-x)))
}

rw.Laplace <- function(sigma, x0, N){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y)/f(x[i-1])))
    x[i] = y else {
    x[i] = x[i-1]
    k = k + 1
  }
}
return(list(x=x, k=k))
}

Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
  psi = as.matrix(psi)
  n = ncol(psi)
  k = nrow(psi)
  psi.means = rowMeans(psi) #row means
  B = n * var(psi.means) #between variance est.
  psi.w = apply(psi, 1, "var") #within variances
  W = mean(psi.w) #within est.
  v.hat = W*(n-1)/n + (B/n) #upper variance est.
  r.hat = v.hat / W #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
set.seed(1234)
n = 12000 #length of chains
k = 4 #number of chains to generate
x0 = c(-10,-5,5,10) #different initial values
sigma = c(0.5, 2, 5) #different variance
b = 1000 #burn-in length
X1 <- X2 <- X3 <- matrix(0,k,n)
rej = matrix(0,3,4)
refline = c(-invF(0.975), invF(0.975))

for(i in 1:k){
  rw1 = rw.Laplace(sigma[1], x0[i], n)
  rw2 = rw.Laplace(sigma[2], x0[i], n)
  rw3 = rw.Laplace(sigma[3], x0[i], n)
  X1[i,]=rw1$x;X2[i,]=rw2$x;X3[i,]=rw3$x
  rej[,i] = c(rw1$k,rw2$k,rw3$k)
}

# compute diagnostic statistics
psi1 <- t(apply(X1, 1, cumsum))
psi2 <- t(apply(X2, 1, cumsum))
psi3 <- t(apply(X3, 1, cumsum))

for (i in 1:k){
  psi1[i,] = psi1[i,] / (1:ncol(psi1))
  psi2[i,] = psi2[i,] / (1:ncol(psi2))
  psi3[i,] = psi3[i,] / (1:ncol(psi3))
}
print(c(Gelman.Rubin(psi1),Gelman.Rubin(psi2),Gelman.Rubin(psi3)))

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
plot(X1[1,1:2000], type="l", xlab=bquote(sigma==0.5), ylab="X", ylim=range(X1[1,]))
abline(h=refline)
plot(X2[1,1:2000], type="l", xlab=bquote(sigma==2), ylab="X", ylim=range(X2[1,]))
abline(h=refline)
plot(X3[1,1:2000], type="l", xlab=bquote(sigma==5), ylab="X", ylim=range(X3[1,]))
abline(h=refline)

## -----------------------------------------------------------------------------
rhat1 <- rhat2 <- rhat3 <- rep(0, n)
for (j in (b+1):n){
  rhat1[j] <- Gelman.Rubin(psi1[,1:j])
  rhat2[j] <- Gelman.Rubin(psi2[,1:j])
  rhat3[j] <- Gelman.Rubin(psi3[,1:j])
}
#par(mfrow=c(2,2))
plot(rhat1[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)
plot(rhat2[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)
plot(rhat3[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
1-apply(rej, 1, mean)/n

## -----------------------------------------------------------------------------
Gibbs.binorm <- function(N,mu1=0,mu2=0,sigma1=1,sigma2=1,rho=0.9,x0=c(mu1,mu2)){
  X = matrix(0, N, 2) #the chain, a bivariate sample
  s1 = sqrt(1-rho^2)*sigma1
  s2 = sqrt(1-rho^2)*sigma2
  X[1,] = x0
  for (i in 2:N) {
    y1 = X[i-1, 2]
    m1 = mu1 + rho * (y1 - mu2) * sigma1/sigma2
    X[i, 1] = rnorm(1, m1, s1)
    x1 = X[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] = rnorm(1, m2, s2)
  }
  return(X)
}

## -----------------------------------------------------------------------------
set.seed(1234)
k = 4 #number of chains to generate
N = 10000 #length of chain
X <- matrix(0,N,2*k)
Xt <- Yt <- matrix(0,k,N)
x0 = matrix(rep(c(-1,-0.5,0.5,1),2),4) #initial values
psix = psiy = matrix(0,k,N)

for(i in 1:k){
  X[,(2*i-1):(2*i)] = Gibbs.binorm(N,x0=x0[i,])
  Xt[i,] = t(X[,2*i-1])
  Yt[i,] = t(X[,2*i])
}
psix = t(apply(Xt, 1, cumsum))
psiy = t(apply(Yt, 1, cumsum))
for (i in 1:k){
  psix[i,] = psix[i,] / (1:ncol(psix))
  psiy[i,] = psiy[i,] / (1:ncol(psiy))
}
print(c(Gelman.Rubin(psix),Gelman.Rubin(psiy)))

## -----------------------------------------------------------------------------
b = 1000 #burn-in length
rhatx = rhaty = rep(0,N)
for (j in (b+1):N){
  rhatx[j] <- Gelman.Rubin(psix[,1:j])
  rhaty[j] <- Gelman.Rubin(psiy[,1:j])
}
#par(mfrow=c(1,2))
plot(rhatx[(b+1):N], type="l", xlab=bquote(X[t]), ylab="R")
abline(h=1.2, lty=2)
plot(rhaty[(b+1):N], type="l", xlab=bquote(Y[t]), ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
x <- X[(b+1):N, 1:2]
plot(x, main="", cex=.5, xlab=bquote(X[t]),
ylab=bquote(Y[t]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
fit = lm(x[,2] ~ x[,1]) #linear model Y=a+bX
summary(fit)

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
T <- function(x,m,y){
  fit1 = lm(m ~ x)
  fit2 = lm(y ~ m+x)
  alpha.hat = fit1$coef[[2]]
  beta.hat = fit2$coef[[2]]
  se1 = summary(fit1)$coef[2,2]
  se2 = summary(fit2)$coef[2,2]
  seab = sqrt(alpha.hat^2*se2+beta.hat^2*se1)
  return(alpha.hat*beta.hat/seab)
}

## -----------------------------------------------------------------------------
R = 100 #times to repeat in permutation
n = 20 #sample size
N = 100 
alpha = c(0,0,1); beta = c(0,1,0); gamma=1
t1err = numeric(3)

## -----------------------------------------------------------------------------
set.seed(1234)
Ts = numeric(R)
p = numeric(N)
for (i in 1:N){
  x1 = rexp(n,2)
  m1 = 1+alpha[1]*x1+rnorm(n)
  y1 = 1+gamma*x1+beta[1]*m1+rnorm(n)
  T0 = T(x1,m1,y1)
  for(j in 1:R){
    m_per = sample(m1, replace=FALSE)
    Ts[j] = T(x1,m_per,y1)
  }
  p[i] = mean(abs(Ts)>=abs(T0))
}
t1err[1] = mean(p<0.05)

## -----------------------------------------------------------------------------
set.seed(1234)
Ts = numeric(R)
p = numeric(N)
for (i in 1:N){
  x2 = rexp(n,2)
  m2 = 1+alpha[2]*x2+rnorm(n)
  y2 = 1+gamma*x2+beta[2]*m2+rnorm(n)
  T0 = T(x2,m2,y2)
  for(j in 1:R){
    x_per = sample(x2, replace=FALSE)
    Ts[j] = T(x_per,m2,y2)
  }
  p[i] = mean(abs(Ts)>=abs(T0))
}
t1err[2] = mean(p<0.05)

## -----------------------------------------------------------------------------
set.seed(1234)
Ts = numeric(R)
p = numeric(N)
for (i in 1:N){
  x3 = rexp(n,2)
  m3 = 1+alpha[3]*x3+rnorm(n)
  y3 = 1+gamma*x3+beta[3]*m3+rnorm(n)
  T0 = T(x3,m3,y3)
  for(j in 1:R){
    y_per = sample(y3, replace=FALSE)
    Ts[j] = T(x3,m3,y_per)
  }
  p[i] = mean(abs(Ts)>=abs(T0))
}
t1err[3] = mean(p<0.05)

## -----------------------------------------------------------------------------
table = data.frame(parameters=c(1,2,3),typeI.error=t1err)
knitr::kable(table, align=rep('c',2))

## -----------------------------------------------------------------------------
solve.alpha <- function(N, b1, b2, b3, f0){
  x1 = rpois(N,1)
  x2 = rexp(N)
  x3 = sample(0:1,N,replace=TRUE)
  g <- function(alpha){
    tmp = exp(-alpha-b1*x1-b2*x2-b3*x3)
    p = 1/(1+tmp)
    mean(p) - f0
  }
  solution = uniroot(g,c(-20,0))
  return(solution$root)
}

## -----------------------------------------------------------------------------
set.seed(0)
N = 1e6; b1=0; b2=1; b3=-1
f0 = c(1e-1,1e-2,1e-3,1e-4)
alpha = numeric(4)
for(i in 1:4){
  alpha[i] = solve.alpha(N, b1, b2, b3, f0[i])
}
alpha

## -----------------------------------------------------------------------------
plot(alpha, f0, pch=16)

## -----------------------------------------------------------------------------
data = c(11,12,8,9,27,28,13,14,16,17,0,1,23,24,10,11,24,25,2,3)
X = matrix(data, nrow=10, byrow=TRUE)
u = X[,1]; v = X[,2]
logL <- function(lambda){
  s = sum(log(exp(-lambda*u)-exp(-lambda*v)))
  return(s)
}
res = optimize(logL, lower=0, upper=5, maximum=TRUE)
lambdahat = res[[1]]

## -----------------------------------------------------------------------------
n = 10
delta = 1
lambda = numeric(100) 
lambda[1] = 1 #initialize 
i = 1
# E-step
Econd <- function(lambda,u,v){ 
  (u*exp(-lambda*u)-v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v))+1/lambda
}
while(delta > 1e-5){
  lambda[i+1] =  n/sum(Econd(lambda[i],u,v)) # M-step
  delta = abs(lambda[i+1]-lambda[i])
  i = i+1
}
lambdahat.EM = lambda[i]
round(c(lambdahat,lambdahat.EM),5)

## -----------------------------------------------------------------------------
x = data.frame(c())
nrow(x);ncol(x)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
head(data.frame(lapply(iris,function(x) if(is.numeric(x)) scale01(x) else x)))

## -----------------------------------------------------------------------------
head(data.frame(lapply(iris[sapply(iris,is.numeric)],scale01)))

## -----------------------------------------------------------------------------
vapply(cars, sd, numeric(1))

## -----------------------------------------------------------------------------
vapply(iris[vapply(iris, is.numeric, logical(1))],
       sd, 
       numeric(1))

## -----------------------------------------------------------------------------
GibbsR <- function(N,mu1=0,mu2=0,sigma1=1,sigma2=1,rho=0.9){
  x0 = c(0,0)
  X = matrix(0, N, 2) #the chain, a bivariate sample
  s1 = sqrt(1-rho^2)*sigma1
  s2 = sqrt(1-rho^2)*sigma2
  X[1,] = x0
  for (i in 2:N) {
    y1 = X[i-1, 2]
    m1 = mu1 + rho * (y1 - mu2) * sigma1/sigma2
    X[i, 1] = rnorm(1, m1, s1)
    x1 = X[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] = rnorm(1, m2, s2)
  }
  return(X)
}

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
#sourceCpp('GibbsC.cpp')
set.seed(0)
N = 5000
b = 1000 #burn-in length
gc = GibbsC(N)
gr = GibbsR(N)
qqplot(gc[(b+1):N,1], gr[(b+1):N,1],xlab=bquote(X[C]),ylab=bquote(X[R]))
qqplot(gc[(b+1):N,2], gr[(b+1):N,2],xlab=bquote(Y[C]),ylab=bquote(Y[R]))

## -----------------------------------------------------------------------------
ts <- microbenchmark(gibbR=GibbsR(N),
gibbC=GibbsC(N))
summary(ts)[,c(1,3,5,6)]

