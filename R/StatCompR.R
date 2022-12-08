#' @title The objective function we want to minimize
#' @description The objective function we want to minimize
#' @param x a positive definite symmetric matrix
#' @param lambda a regularization parameter
#' @param s the sample covariance
#' @return d value of the objective function
#' @examples
#' \dontrun{
#' f(diag(3),1,diag(3))
#' }
#' @export
f <- function(x,lambda,s){
  d = sum(diag(s %*% x)) - log(det(x)) + lambda/2 * norm(x,"F")
  return(d)
}

#' @title Compute first-order gradient
#' @description Compute first-order gradient of the objective function
#' @param x a positive definite symmetric matrix
#' @param lambda a regularization parameter
#' @param s the sample covariance
#' @return d the first-order gradient
#' @examples
#' \dontrun{
#' f_gra1(diag(3),1,diag(3))
#' }
#' @export
f_gra1 <- function(x,lambda,s){
  d = s - solve(x) + lambda * x
  return(d)
}

#' @title Compute second-order gradient
#' @description Compute second-order gradient of the objective function
#' @param x a positive definite symmetric matrix
#' @param lambda a regularization parameter
#' @param s the sample covariance
#' @return d the second-order gradient
#' @examples
#' \dontrun{
#' f_gra2(diag(3),1,diag(3))
#' }
#' @export
f_gra2 <- function(x,lambda,s){
  len = dim(x)[1]
  d = kronecker(solve(x),solve(x)) + lambda * diag(len^2)
  return(d)
} 

#' @title True negative rate
#' @description Compute true negative rate of the estimation
#' @param s the ground truth inverse of the covariance
#' @param x the estimated inverse of the covariance
#' @return k1/k2 true negative rate
#' @examples
#' \dontrun{
#' s = diag(3)
#' x = matrix(c(1:9),3)
#' TNR(s,x)
#' }
#' @export
TNR <- function(s,x){
  isclose0 <- function(x){
    return(isTRUE(all.equal(x,0)))
  }
  k1 = sum(s == 0 & apply(x,c(1,2),isclose0))
  k2 = sum(s == 0)
  return(k1/k2)
}

#' @title True positive rate
#' @description Compute true positive rate of the estimation
#' @param s the ground truth inverse of the covariance
#' @param x the estimated inverse of the covariance
#' @return k1/k2 true positive rate
#' @examples
#' \dontrun{
#' s = diag(3)
#' x = matrix(c(1:9),3)
#' TPR(s,x)
#' }
#' @export
TPR <- function(s,x){
  isclose0 <- function(x){
    return(isTRUE(all.equal(x,0)))
  }
  k1 = sum(s != 0 & !apply(x,c(1,2),isclose0))
  k2 = sum(s != 0) 
  return(k1/k2)
}
  
#' @title The log-likelihood performance metric
#' @description Compute the log-likelihood performance metric
#' @param s the ground truth covariance
#' @param x the estimated inverse of the covariance
#' @return d the log-likelihood performance metric
#' @examples
#' \dontrun{
#' s = diag(3)
#' x = matrix(c(1:9),3)
#' likelihood(s,x)
#' }
#' @export
likelihood <- function(s,x){
  d = -(sum(diag(s %*% x)) - log(det(x)))
  return(d)
}

#' @title The normalized entropy loss
#' @description Compute the normalized entropy loss
#' @param s the ground truth covariance
#' @param x the estimated inverse of the covariance
#' @param n number of variables
#' @return d the normalized entropy loss
#' @examples
#' \dontrun{
#' s = diag(3)
#' x = matrix(c(1:9),3)
#' LossE(s,x,3)
#' }
#' @export
LossE <- function(s,x,n){
  d = (1/n)*(sum(diag(s %*% x))-log(det(s %*% x))-n)
  return(d)
}

#' @title The quadratic loss
#' @description Compute the quadratic loss
#' @param s the ground truth covariance
#' @param x the estimated inverse of the covariance
#' @param n number of variables
#' @return d the quadratic loss
#' @examples
#' \dontrun{
#' s = diag(3)
#' x = matrix(c(1:9),3)
#' LossQ(s,x,3)
#' }
#' @export
LossQ <- function(s,x,n){
  d = 1/n*norm(s %*% x - diag(n),"F")
  return(d)
}

#' @title Generate a synthetic data
#' @description Generate a synthetic data. Specificly, A[i,i]=1, A[i,i-1]=A[i-1,i]=0.5, A[i,i-2]=A[i-2,i]=0.25 and zero otherwise.
#' @param n number of variables
#' @return A the synthetic data
#' @examples
#' \dontrun{
#' AR2(30)
#' }
#' @export
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

#' @title A projection operator
#' @description Project a matrix into a space Omiga
#' @param x the operated matrix
#' @param kap the l0 constraint
#' @param omiga the space
#' @return a projection of the operated matrix
#' @examples
#' \dontrun{
#' x = matrix(c(1:9),3)
#' omiga = c(3,6,7,8)
#' pp(x,5,omiga)
#' }
#' @export
pp <- function(x,kap,omiga){
  n = nrow(x)
  y = x
  y[omiga] = 0
  vecy = as.vector(y)
  num1 = abs(vecy)
  num_s = sort.int(num1,decreasing = TRUE,index.return = TRUE)
  mat0 <- matrix(0,n,n)
  for (s in 1:kap){
    l = num_s$ix[s]
    l1 = ceiling(l/n)
    l0 = l-(l1-1)*n
    mat0[l0,l1] = num_s$x[s] 
  }
  return(mat0)
}


#' @title Obtain an approximate Newton direction
#' @description Using a conjugate gradient method to compute approxiamte Newton direction.
#' @param g the first-order gradient of the objective function 
#' @param H the second-order gradient of the objective function or we can say Hessian matrix
#' @param wk the working set, which all but a working set of the matrix entries will not be fixed to zero
#' @return an approximate Newton direction
#' @examples
#' \dontrun{
#' x= matrix(c(2,1,1,1,3,2,1,2,4),3)
#' cgd(f_gra1(x,1,diag(3)),f_gra2(x,1,diag(3)),c(1,5,6,8,9))
#' }
#' @export
cgd <- function(g,H,wk){
  k = 1
  b = as.vector(t(g))
  n = length(b)
  b[-wk] = 0
  v0 = rep(0, n)
  x <- list(); r <- list(); p <- list(); y <- list()
  x[[1]] = v0
  r[[1]] = b
  p[[1]] = b 
  while (k <= n+1 & max(r[[k]]) > 0.05){
    if(identical(r[[k]],v0) == TRUE){break} 
    y[[k]] = H %*% p[[k]]
    y[[k]][-wk] = 0
    alpha = as.numeric((t(r[[k]])%*%r[[k]])/(t(p[[k]])%*%y[[k]]))
    x[[k+1]] = x[[k]] + alpha*p[[k]]
    r[[k+1]] = r[[k]] - alpha*y[[k]]
    beta = as.numeric((t(r[[k+1]])%*%r[[k+1]])/(t(r[[k]])%*%r[[k]]))
    p[[k+1]] = r[[k+1]] + beta*p[[k]]
    k = k+1
  }
  D = matrix(x[[k]],nrow(g))
  return(D)
}

#' @title Approximate Newton Algorithm
#' @description an approximate Newton algorithm to compute the minimum point of the objective function
#' @param stop.coef a stopping criteria for the algorithm to terminate while loop
#' @param L0 the cholesky decomposition of the initial point 
#' @param kk the l0 constraint, i.e. a maximally allowable number of nonzeros
#' @param omiga a sparsity constraint set which includes the indexes of variables known to be conditionally independent
#' @param lambda a regularization parameter
#' @param sigma a parameter takes value in (0,1) to ensure the stepsize is sufficiently small
#' @param elta a positive parameter to ensure the objective function value is sufficiently decreased
#' @param alpha_max the maximum stepsize for Newton's method
#' @param alpha_min the minimum stepsize for Newton's method
#' @param alpha_min_newton the first-order gradient of the objective function
#' @param epson the first-order gradient of the objective function
#' @param l the index of a decreased sequence 
#' @param p the remainder of the iteration k mod q
#' @param q every q − 1 iterations of the Newton steps we perform a single iteration of the gradient projection step
#' @param s the sample covariance
#' @return the cholesky decompositon of the minimum point
#' @export
ANA <- function(stop.coef,L0,kk,omiga,lambda,sigma,elta,alpha_max,alpha_min,alpha_min_newton,epson,l,p,q,s){
  #stepsize alpha_bb
  alpha_bb <- function(x0,x1,lambda,s){
    a = t(f_gra1(x1,lambda,s) - f_gra1(x0,lambda,s))
    alp_b = sum(diag(a %*% (x1-x0)))/(norm(x1-x0,"F")^2)
    return(alp_b)
  }
  
  k = 1
  x1 = t(L0)%*%L0
  x0 = x1 + (10^(-l))*f_gra1(x1,lambda,s)
  L1 = L0
  L0 = chol(x0)
  o1 = 1; o2 = 0
  grad_step = TRUE
  while (o1 >= o2){ 
    if (k%%q!=p) {
      if (grad_step==TRUE){
        wk = which(x1!=0 | abs(x0-x1)>=epson)
      }else{
        Y = pp(x1 - alpha * f_gra1(x1,lambda,s),kk,omiga)
        wk = which(x1!=0 | abs(Y-x1)>=epson)
      } 
      D = cgd(f_gra1(x1,lambda,s),f_gra2(x1,lambda,s),wk)
      
      alpha = 1
      alpha_0 = 1
      j = 0
      ls = TRUE
      while (ls==TRUE & alpha >= alpha_min_newton) {
        alpha = alpha_0*sigma^j
        x_tr = pp(x1-alpha*D,kk,omiga)
        fit <- try(chol(x_tr),silent = TRUE)
        if ("try-error"%in%class(fit)){
          j = j+1
        }else{
          L_tr = chol(x_tr)
          if (f(x_tr,lambda,s)>max(f(x0,lambda,s),f(x1,lambda,s))-elta/2*sum(diag(t(x_tr-x1)%*%(x_tr-x1)))){
            j = j+1
            if (alpha<alpha_min){
              break
            }
          }else{
            ls = FALSE
            L2 = L_tr
            L0 = L1
            L1 = L2
            x0 = t(L0)%*%L0
            x1 = t(L1)%*%L1
          }
        }
      }
    }
    grad_step = FALSE
    alpha = 1
    if (k%%q==p | alpha<alpha_min_newton){
      G = f_gra1(x1,lambda,s)
      alpha_0 = min(alpha_max,max(alpha_min,alpha_bb(x0,x1,lambda,s)))
      j = 0
      ls = TRUE
      while (ls == TRUE) {
        alpha = alpha_0*sigma^j
        x_tr = pp(x1-alpha*G,kk,omiga)
        fit <- try(chol(x_tr),silent = TRUE)
        if ("try-error"%in%class(fit)){
          j = j+1
        }else{
          L_tr = chol(x_tr)
          if (f(x_tr,lambda,s)>max(f(x0,lambda,s),f(x1,lambda,s))-elta/2*sum(diag(t(x_tr-x1)%*%(x_tr-x1)))){
            j = j+1
            if (alpha<alpha_min){
              break
            }
          }else{
            ls = FALSE
            L2 = L_tr
            L0 = L1
            L1 = L2
            x0 = t(L0)%*%L0
            x1 = t(L1)%*%L1
          }
        }
      }
      grad_step = TRUE
    }
    k = k+1
    o1 = abs(f(x1,lambda,s)-f(x0,lambda,s))
    o2 = max(1,abs(f(x0,lambda,s)))*stop.coef
  }
  return(L1)
}

#' @title Homotopy Approximate Newton Algorithm
#' @description a homotopy approximate Newton algorithm to compute MLE of the inverse covariance
#' @param n the number of variables for a multinormal distribution
#' @param kkk the l0 constraint, i.e. a maximally allowable number of nonzeros
#' @param omiga a sparsity constraint set which includes the indexes of variables known to be conditionally independent
#' @param s the sample covariance
#' @param epson the first-order gradient of the objective function
#' @param lambda_L the last point of the decreased sequence
#' @param lambda_0 the first point of the decreased sequence
#' @param sigma a parameter takes value in (0,1) to ensure the stepsize is sufficiently small
#' @param elta a positive parameter to ensure the objective function value is sufficiently decreased
#' @param alpha_min the minimum stepsize for Newton's method
#' @param alpha_max the maximum stepsize for Newton's method
#' @param alpha_min_newton the first-order gradient of the objective function
#' @param L length of a decreased sequence 
#' @param p the remainder of the iteration k mod q
#' @param q every q − 1 iterations of the Newton steps we perform a single iteration of the gradient projection step
#' @return the cholesky decompositon of the minimum point
#' @examples
#' \dontrun{
#' library(MASS)
#' n = 100
#' A = AR2(n)
#' m = rep(0,n)
#' data <- mvrnorm(2*n,m,solve(A))
#' s = var(data)
#' omiga = which(A==0)
#' x = HANA(n,494,omiga,s)
#' t(x) %*% x
#' }
#' @export
HANA <- function(n,kkk,omiga,s,epson = 1e-2,lambda_L=1e-5,lambda_0=1e-1,sigma=0.5,elta = 1e-6,alpha_min=1e-15,alpha_max =1e15,alpha_min_newton=1e-5,L=10,p=1,q=5){
  t1 = Sys.time()
  stop.coef = 1e-2
  x0 = solve(diag(diag(s))+0.01*diag(n)) 
  L0 = chol(x0)
  kk = seq(from=2*kkk,to=kkk,length.out=L+1)
  lambda_L_log = log10(lambda_L)
  lambda_0_log = log10(lambda_0)
  lambda_log = seq(from=lambda_0_log,to = lambda_L_log,length.out = L+1)
  lambda = 10^lambda_log
  for (l in 1:L+1) {
    if (l==L+1){
      stop.coef=1e-5
    }
    L1 = ANA(stop.coef,L0,kk[l],omiga,lambda[l],sigma,elta,alpha_max,alpha_min,alpha_min_newton,epson,l,p,q,s)
    L0 = L1
  }
  t2 = Sys.time()
  print(t2-t1)
  return(L1)
}

