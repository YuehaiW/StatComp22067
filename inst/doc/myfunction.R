## -----------------------------------------------------------------------------
library(StatComp22067)

## ----eval=FALSE---------------------------------------------------------------
#  pp <- function(x,kap,omiga){
#     n = nrow(x)
#     y = x
#     y[omiga] = 0
#     vecy = as.vector(y)
#     num1 = abs(vecy)
#     num_s = sort.int(num1,decreasing = TRUE,index.return = TRUE)
#     mat0 <- matrix(0,n,n)
#     for (s in 1:kap){
#       l = num_s$ix[s]
#       l1 = ceiling(l/n)
#       l0 = l-(l1-1)*n
#       mat0[l0,l1] = num_s$x[s]
#     }
#     return(mat0)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  cgd <- function(g,H,wk){
#    k = 1
#    b = as.vector(t(g))
#    n = length(b)
#    b[-wk] = 0
#    v0 = rep(0, n)
#    x <- list(); r <- list(); p <- list(); y <- list()
#    x[[1]] = v0
#    r[[1]] = b
#    p[[1]] = b
#    while (k <= n+1 & max(r[[k]]) > 0.05){
#      if(identical(r[[k]],v0) == TRUE){break}
#      y[[k]] = H %*% p[[k]]
#      y[[k]][-wk] = 0
#      alpha = as.numeric((t(r[[k]])%*%r[[k]])/(t(p[[k]])%*%y[[k]]))
#      x[[k+1]] = x[[k]] + alpha*p[[k]]
#      r[[k+1]] = r[[k]] - alpha*y[[k]]
#      beta = as.numeric((t(r[[k+1]])%*%r[[k+1]])/(t(r[[k]])%*%r[[k]]))
#      p[[k+1]] = r[[k+1]] + beta*p[[k]]
#      k = k+1
#    }
#    D = matrix(x[[k]],nrow(g))
#    return(D)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  arma::mat conjugate(arma::mat A, arma::vec b) {
#    int k = 1;
#    int n = b.n_elem;
#    arma::mat v0 = arma::zeros(n,1);
#    arma::mat x0 = v0;
#    arma::mat x1 = v0;
#    arma::mat p1 = v0;
#    arma::mat r0 = v0;
#    arma::mat r1 = v0;
#    arma::mat p0 = v0;
#    r0.col(0) = -b;
#    p0.col(0) = b;
#    while(k <= n+1 && arma::norm(r0,"inf") > 0.05){
#      arma::mat alpha = r0.t() * r0 / (p0.t() * A * p0);
#      x1 = x0 + alpha(0,0) * p0;
#      r1 = r0 + alpha(0,0) * A * p0;
#      arma::mat beta = r1.t() * r1 / (r0.t() * r0);
#      p1 = -r1 + beta(0,0) * p0;
#      p0 = p1;
#      r0 = r1;
#      x0 = x1;
#      k += 1;
#    }
#    return(x0);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  ANA <- function(stop.coef,L0,kk,omiga,lambda,sigma,elta,alpha_max,alpha_min,alpha_min_newton,epson,l,p,q,s){
#    #stepsize alpha_bb
#    alpha_bb <- function(x0,x1,lambda,s){
#      a = t(f_gra1(x1,lambda,s) - f_gra1(x0,lambda,s))
#      alp_b = sum(diag(a %*% (x1-x0)))/(norm(x1-x0,"F")^2)
#      return(alp_b)
#    }
#  
#    k = 1
#    x1 = t(L0)%*%L0
#    x0 = x1 + (10^(-l))*f_gra1(x1,lambda,s)
#    L1 = L0
#    L0 = chol(x0)
#    o1 = 1; o2 = 0
#    grad_step = TRUE
#    while (o1 >= o2){
#      if (k%%q!=p) {
#        if (grad_step==TRUE){
#           wk = which(x1!=0 | abs(x0-x1)>=epson)
#        }else{
#           Y = pp(x1 - alpha * f_gra1(x1,lambda,s),kk,omiga)
#           wk = which(x1!=0 | abs(Y-x1)>=epson)
#        }
#        D = cgd(f_gra1(x1,lambda,s),f_gra2(x1,lambda,s),wk)
#  
#        alpha = 1
#        alpha_0 = 1
#        j = 0
#        ls = TRUE
#        while (ls==TRUE & alpha >= alpha_min_newton) {
#          alpha = alpha_0*sigma^j
#          x_tr = pp(x1-alpha*D,kk,omiga)
#          fit <- try(chol(x_tr),silent = TRUE)
#          if ("try-error"%in%class(fit)){
#            j = j+1
#          }else{
#            L_tr = chol(x_tr)
#            if (f(x_tr,lambda,s)>max(f(x0,lambda,s),f(x1,lambda,s))-elta/2*sum(diag(t(x_tr-x1)%*%(x_tr-x1)))){
#             j = j+1
#             if (alpha<alpha_min){
#               break
#             }
#            }else{
#              ls = FALSE
#              L2 = L_tr
#              L0 = L1
#              L1 = L2
#              x0 = t(L0)%*%L0
#              x1 = t(L1)%*%L1
#            }
#           }
#        }
#      }
#      grad_step = FALSE
#      alpha = 1
#      if (k%%q==p | alpha<alpha_min_newton){
#        G = f_gra1(x1,lambda,s)
#        alpha_0 = min(alpha_max,max(alpha_min,alpha_bb(x0,x1,lambda,s)))
#        j = 0
#        ls = TRUE
#        while (ls == TRUE) {
#          alpha = alpha_0*sigma^j
#          x_tr = pp(x1-alpha*G,kk,omiga)
#          fit <- try(chol(x_tr),silent = TRUE)
#          if ("try-error"%in%class(fit)){
#            j = j+1
#          }else{
#            L_tr = chol(x_tr)
#            if (f(x_tr,lambda,s)>max(f(x0,lambda,s),f(x1,lambda,s))-elta/2*sum(diag(t(x_tr-x1)%*%(x_tr-x1)))){
#             j = j+1
#             if (alpha<alpha_min){
#               break
#             }
#            }else{
#              ls = FALSE
#              L2 = L_tr
#              L0 = L1
#              L1 = L2
#              x0 = t(L0)%*%L0
#              x1 = t(L1)%*%L1
#            }
#          }
#        }
#        grad_step = TRUE
#      }
#       k = k+1
#       o1 = abs(f(x1,lambda,s)-f(x0,lambda,s))
#       o2 = max(1,abs(f(x0,lambda,s)))*stop.coef
#    }
#      return(L1)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  HANA <- function(n,kkk,omiga,s,epson = 1e-2,lambda_L=1e-5,lambda_0=1e-1,sigma=0.5,elta = 1e-6,alpha_min=1e-15,alpha_max =1e15,alpha_min_newton=1e-5,L=10,p=1,q=5){
#    t1 = Sys.time()
#    stop.coef = 1e-2
#    x0 = solve(diag(diag(s))+0.01*diag(n))
#    L0 = chol(x0)
#    kk = seq(from=2*kkk,to=kkk,length.out=L+1)
#    lambda_L_log = log10(lambda_L)
#    lambda_0_log = log10(lambda_0)
#    lambda_log = seq(from=lambda_0_log,to = lambda_L_log,length.out = L+1)
#    lambda = 10^lambda_log
#    for (l in 1:L+1) {
#      if (l==L+1){
#        stop.coef=1e-5
#      }
#     L1 = ANA(stop.coef,L0,kk[l],omiga,lambda[l],sigma,elta,alpha_max,alpha_min,alpha_min_newton,epson,l,p,q,s)
#     L0 = L1
#     }
#     t2 = Sys.time()
#     print(t2-t1)
#     return(L1)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  TNR <- function(s,x){
#    isclose0 <- function(x){
#      return(isTRUE(all.equal(x,0)))
#    }
#    k1 = sum(s == 0 & apply(x,c(1,2),isclose0))
#    k2 = sum(s == 0)
#    return(k1/k2)
#  }
#  
#  TPR <- function(s,x){
#    isclose0 <- function(x){
#      return(isTRUE(all.equal(x,0)))
#    }
#    k1 = sum(s != 0 & !apply(x,c(1,2),isclose0))
#    k2 = sum(s != 0)
#    return(k1/k2)
#  }
#  
#  likelihood <- function(s,x){
#    d = log(det(x)) - sum(diag(s %*% x))
#    return(d)
#  }
#  
#  LossE <- function(s,x,n){
#    d = (1/n)*(sum(diag(s %*% x))-log(det(s %*% x))-n)
#    return(d)
#  }
#  
#  LossQ <- function(s,x,n){
#    d = 1/n*norm(s %*% x - diag(n),"F")
#    return(d)
#  }

## -----------------------------------------------------------------------------
AR2 <- function(n){
  A = diag(n)
  for (i in 1:n){
    for (j in 1:n){
      if (abs(i-j)==1) {
        A[i,j]=0.5
      }else{
        if (abs(i-j)==2){
          A[i,j] = 0.25
        }
      }
    }
  }
  return(A)
}

## ----eval=TRUE----------------------------------------------------------------
library(MASS)
set.seed(1)
n = 100
A = AR2(n)
m = rep(0,n)
data <- mvrnorm(2*n,m,solve(A))
s = var(data)
omiga = integer(0)
L = HANA(n,494,omiga,s)
x = t(L) %*% L

## ----eval=TRUE----------------------------------------------------------------
TNR(A,x)
TPR(A,x)
LossQ(solve(A),x,n)
LossE(solve(A),x,n)
likelihood(s,x)

