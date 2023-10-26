
# Sampling ----------------------------------------------------------------

# C's
auxiliar <- function (p, alpha)
{
  if (any(sort(p) != p)) {
    stop(paste("p must be crecient", "\n", ""))
  }
  if (any(alpha <= 1)) {
    stop(paste("alpha must be > 1", "\n", ""))
  }
  
  aux <- c(1, rep(NA, length(p) - 1), 0)
  for (i in 2:length(p)) {
    aux[i] <- prod((p[i] / p[i - 1]) ^ (1 - alpha[i - 1]), aux[(i - 1)])
  }
  
  return (aux)
}



# density distribution function
dpwpowerlaw <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar(p, alpha)
  j <- max(which(p < x)) #antes pusimos min
  d <- ( (alpha[j] -1)/p[j] )*(x / p[j])^(-alpha[j]) * C[j]

  return(d)
}



# cumulative distribution function
ppwpowerlaw <- function (x, p, alpha)
{
  C <- auxiliar(p, alpha)
  j <- max(which(p <= x)) #antes pusimos min
  u <- 1 - (x / p[j]) ^ (1 - alpha[j]) * C[j]
  
  return(u)
}

# survival function
spwpowerlaw <- function (x, p, alpha)
{
  s <- 1 - ppwpowerlaw(x, p, alpha)
  
  return(s)
}

# hazard function
hpwpowerlaw <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar(p, alpha)
  j <- max(which(p <= x))
  (alpha[j] - 1) / x
}

# quantile function
qpwpowerlaw <- function (u, p, alpha)
{
  C <- auxiliar(p, alpha)
  aux <- u < (1 - C)
  j <- min(which(aux == TRUE)) - 1
  x <- p[j] * ((1 - u) / C[j]) ^ (1 / (1 - alpha[j]))
  
  return (x)
}

# sampling piece-wise power law
rpwpowerlaw <- function(n, p, alpha)
{
  flag <- 0
  while(flag == 0)
	{
    flag <- 0
    x <- vector(length = n)
	  C <- auxiliar(p, alpha)
	  for (k in 1:n) 
	  {
	    u <- runif(1)
	    aux <- u < (1 - C)
	    j <- min(which(aux == TRUE)) - 1
	    x[k] <- p[j] * ((1 - u) / C[j]) ^ (1 / (1 - alpha[j]))
	  }
	  nj <- n_each_interval(x, p)
	  if (all(nj > 2)) flag <- 1  # all(nj) != 0
	}
  
  return (x)
}

## Simulation under type I censure
censI_pwpowerlaw <- function (n, p, alpha, p_cens)
{
  x <- rpwpowerlaw (n, p, alpha)
  t.c <- qpwpowerlaw (1 - p_cens, p, alpha)
  cens <- rep(1, n)
  for (i in 1:n) {
    if (x[i] >= t.c) {
      x[i] <- t.c
      cens[i] <- 0
    }
  }
  
  return (list(x = x, cens = cens))
}

# identifying the observations in each interval
index_each_interval <- function (x, p)
{
  k <- length(p)
  nj <- vector("list", length = k)
  for (j in 1:k)
  {
    if (j < k)
    {
      nj[[j]] <- which(x >= p[j] & x < p[j + 1])
    } else {
      nj[[j]] <- which(x >= p[j])
    }
  }
  
  return (nj)
}

# counting the observations in each interval
n_each_interval <- function (x, p)
{
  nj <- lengths(index_each_interval(x, p))
  return (nj)
}

mle_pwpowerlaw <- function (x, p, bias.correction = TRUE)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x, p)
  nj <- n_each_interval(x, p)
  if (any(nj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum(log(x[index[[j]]] / p[j]))
    b <- sum(nj[(j + 1):k] * log(p[(j + 1)] / p[j]))
    if (bias.correction) {
      alpha[j] <- 1 + (nj[j] - 1) * (a + b) ^ (-1)
    } else {
      alpha[j] <- 1 + nj[j] * (a + b) ^ (-1)
    }
  }
  
  if (bias.correction) {
    alpha[k] <- 1 + (nj[k] - 1) * sum(log(x[index[[k]]] / p[k])) ^ (-1)
  } else {
    alpha[k] <- 1 + nj[k] * sum(log(x[index[[k]]] / p[k])) ^ (-1)
  }
  
  return(alpha)
}

## logC[j]
log_auxiliar <- function (p, alpha)
{
  k <- length(p) - 1
  logC <- vector(length = k)
  logC[1] <- 0
  for (j in 1:k)
  {
    logC[j + 1] <-
      sum((1 - alpha[1:j]) * (log(p[2:(j + 1)]) - log(p[1:j])))
  }
  
  return (logC)
}

## sum(logdata) in each I[j]
sum_logdata <- function (x, p)
{
  k <- length(p)
  index <- index_each_interval(x, p)
  logdata <- vector(length = k)
  for (j in 1:k)
  {
    logdata[j] <- sum(log(x[index[[j]]]))
  }
  
  return (logdata)
}

## loglikelihood
loglik_pwpowerlaw <- function(theta, x, p)
{
  logC <- log_auxiliar(p, theta)
  logdata <- sum_logdata(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (log(theta - 1) - (1 - theta) * log(p) + logC)) - sum(theta * logdata)
  
  return (l)
}

## loglikelihood reparametrized
loglik_pwpowerlaw_reparametrized <- function(theta, x, p)
{
  arg <- exp(theta) + 1
  logC <- log_auxiliar(p, arg)
  logdata <- sum_logdata(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (log(arg - 1) - (1 - arg) * log(p) + logC)) - sum(arg * logdata)
  
  return (l)
}

r_moment<-function(r, alpha, p)
{
	k<-length(p); C<-auxiliar(p,alpha)[1:length(p)]
	## First k-1 integrals
	Ij<-vector(length=k-1)
	for(j in 1:(k-1))
	{
		if(alpha[j]!=r+1){ Ij[j]<-C[j]*(alpha[j]-1)/(p[j])^(1-alpha[j])*(  1/(r+1-alpha[j])*( p[j+1]^(r+1-alpha[j]) -p[j]^(r+1-alpha[j]) )  ) }
		if(alpha[j]==r+1){ Ij[j]<-C[j]*(alpha[j]-1)/(p[j])^(1-alpha[j])*( log(p[j+1]) -log(p[j]) ) }
	}
	
	## k-ith integral
	Ik<- -(alpha[k] -1)*C[k]*p[k]^r/(r+1-alpha[k])
	
	## r-moment
	sum(Ij) +Ik
}

# Simulation under random censoring
censrand_pwpowerlaw <- function(n, p, alpha, PI)
{
  flag <- 0
  lambda <- r_moment(1, alpha, p)/PI
  x <- vector(length = n)
  delta <- vector(length = n)
  C <- auxiliar(p, alpha)
  
  while(flag == 0)
  {

  for (k in 1:n) 
  {
    u <- runif(1)
    aux <- u < (1 - C)
    j <- min(which(aux == TRUE)) - 1
    aux <- p[j] * ((1 - u) / C[j]) ^ (1 / (1 - alpha[j]))
    c <- runif(1, p[1], lambda) #runif(1, 0, lambda)
    if (aux < c) {
      delta[k] <- 1
      x[k] <- aux
    } else {
      delta[k] <- 0
      x[k] <- c
    }
    
  }
  
  out <- data.frame(time = x, status = delta)
  dj <- n_each_interval(out$time[out$status == 1], p)
  if (all(dj > 2)) flag <- 1  # all(nj) != 0 
  }
  
  return (out)
}

diff_c <- function(x) c(x[1] - x[2], x[2] - x[3], x[3] - x[4])

# explicit mle censored
mle_cens_pwpowerlaw <- function (x, p, bias.correction = T)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  if (any(dj <= 2)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum(log(x$time[index[[j]]] / p[j]))
    b <- sum(nj[(j + 1):k] * log(p[(j + 1)] / p[j]))
    alpha[j] <- 1 + dj[j] * (a + b) ^ (-1)
  }
  
  alpha[k] <- 1 + dj[k] * sum(log(x$time[index[[k]]] / p[k])) ^ (-1)
  
  if (bias.correction) {
    h_ij <- diag(-sum(nj)*diff_c(auxiliar(p, alpha)) / (alpha - 1)^2) # K
    h_ij1 <- diag(c(2*sum(nj) * diff_c(auxiliar(p, alpha))[1] / (alpha[1] - 1)^3, 0, 0))
    h_ij2 <- diag(c(0, 2*sum(nj)*diff_c(auxiliar(p, alpha))[2] / (alpha[2] - 1)^3, 0))
    h_ij3 <- diag(c(0, 0, 2*sum(nj)*diff_c(auxiliar(p, alpha))[3] / (alpha[3] - 1)^3))
    
    f_temp <- function(x) -sum(nj)*diff_c(auxiliar(p, x)) / (x - 1)^2
    temp <- jacobian(f_temp, alpha)
    h_ij_1 <- diag(temp[,1])
    h_ij_2 <- diag(temp[,2])
    h_ij_3 <- diag(temp[,3])
    
    A1 <- h_ij_1 - 0.5 * h_ij1
    A2 <- h_ij_2 - 0.5 * h_ij2
    A3 <- h_ij_3 - 0.5 * h_ij3
    A <- cbind(A1, A2, A3)
    
    Bias <- c(solve(h_ij) %*% A %*% fBasics::vec(solve(h_ij)))
    alpha <- alpha - Bias
  }
  
  return(alpha)
}

#profile_log_likelihood to estimate the change points
profile_loglik_pwpowerlaw <- function(p, x)
{
  theta <- mle_pwpowerlaw(x, p)
  logC <- log_auxiliar(p, theta)
  logdata <- sum_logdata(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (log(theta - 1) - (1 - theta) * log(p) + logC)) - sum(theta * logdata)
  
  return (l)
}

loglik_cens_pwpowerlaw <- function(theta, x, p)
{
  logC <- log(auxiliar(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * log(theta - 1) +
          nj * logC - nj * (1 - theta) * log(p)) +
    sum((1 - theta) * logdata) - sum(log(x$time[x$status==1]))
  
  return (l)
}
