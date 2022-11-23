
# HMM functions -----------------------------------------------------------

zgamma_scalar = function(x, mu, sigma, p){
  if(x == 0){
    return(p)
  }
  else{
    return((1-p)*dgamma(x, shape = mu^2/sigma^2, scale = sigma^2/mu))
  }
}

zgamma = function(x, mu, sigma, p){
  return(sapply(x, zgamma_scalar, mu = mu, sigma = sigma, p = p))
}

flip_gamma = function(x, mu, sigma){
  return(dgamma(-x, shape = mu^2/sigma^2, scale = sigma^2/-mu))
}


zbeta_scalar = function(x, shape1, shape2, p){
  if(x == 0){
    return(p)
  }
  else{
    return((1-p)*dbeta(x, shape1 = shape1, shape2 = shape2))
  }
}

zbeta = function(x, shape1, shape2, p){
  return(sapply(x, zbeta_scalar, shape1 = shape1, shape2 = shape2, p = p))
}




# likelihood for 3D-time series assuming contemporaneous conditional independence
mllk = function(theta.star, X, N){
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) # stationary, smaller tolerance to avoid numerical problems
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  p.b = plogis(theta.star[(N-1)*N+4*N+1:N])
  
  # normal distribution: Height first difference
  mu = theta.star[(N-1)*N+5*N+1:N] # means of normal distributions
  sigma = exp(theta.star[(N-1)*N+6*N+1:N]) # sds of normal distributions
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind,j] =
      dgamma(X$step[ind], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      zbeta(X$angle[ind], shape1 = alpha[j], shape2 = beta[j], p = p.b[j])* # beta
      dnorm(X$height.fd[ind], mean = mu[j], sd = sigma[j]) # normal
  }
  
  # forward algorithm to compute the log-likelihood
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo/sum(foo)
  for (t in 2:nrow(X)){
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l+log(sum(foo))
    phi = foo/sum(foo)
  }
  
  return(-l) 
}


mllk_cov = function(theta.star, X, N){
  # avoid zero beta
  X$angle[which(X$angle == 0)] = runif(length(X$angle[which(X$angle == 0)]), 0, 0.03)
  
  # get betas
  coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # normal distribution: Height first difference
  mu = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
  sigma = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
  
  delta = c(1, exp(theta.star[2*(N-1)*N+6*N+1:(N-1)])) # initial distribution
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind,j] =
      dgamma(X$step[ind], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind], shape1 = alpha[j], shape2 = beta[j])* # beta
      dnorm(X$height.fd[ind], mean = mu[j], sd = sigma[j]) # normal
  }
  
  # forward algorithm to compute the log-likelihood
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo/sum(foo)
  for (t in 2:nrow(X)){
    eta = coef[,1] + coef[,2]*X$temp[t]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l+log(sum(foo))
    phi = foo/sum(foo)
  }
  
  return(-l) 
}



get_probs = function(theta.star, X, N){
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-20) # stationary, smaller tolerance to avoid numerical problems
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  p.b = plogis(theta.star[(N-1)*N+4*N+1:N])
  
  # normal distribution: Height first difference
  mu = theta.star[(N-1)*N+5*N+1:N] # means of normal distributions
  sigma = exp(theta.star[(N-1)*N+6*N+1:N]) # sds of normal distributions
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind,j] =
      dgamma(X$step[ind], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      zbeta(X$angle[ind], shape1 = alpha[j], shape2 = beta[j], p = p.b[j])* # beta
      dnorm(X$height.fd[ind], mean = mu[j], sd = sigma[j]) # normal
  }
  
  probs = matrix(data = NA, nrow = nrow(X), N)
  
  # forward algorithm to compute the log-likelihood
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo/sum(foo)
  probs[1,] = phi/sum(phi)
  for (t in 2:nrow(X)){
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l+log(sum(foo))
    phi = foo/sum(foo)
    probs[t,] = phi/sum(phi)
  }
  
  return(probs) 
}


viterbi = function(theta.star, X, N){ 
  n = nrow(X)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-20) # stationary, smaller tolerance to avoid numerical problems
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  p.b = plogis(theta.star[(N-1)*N+4*N+1:N])
  
  # normal distribution: Height first difference
  mu = theta.star[(N-1)*N+5*N+1:N] # means of normal distributions
  sigma = exp(theta.star[(N-1)*N+6*N+1:N]) # sds of normal distributions
  
  allprobs = matrix(1, n, N)
  ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind,j] =
      dgamma(X$step[ind], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      zbeta(X$angle[ind], shape1 = alpha[j], shape2 = beta[j], p = p.b[j])* # beta
      dnorm(X$height.fd[ind], mean = mu[j], sd = sigma[j]) # normal
  }
  
  xi = matrix(0, n, ncol = N) 
  foo = delta * allprobs[1, ]
  xi[1, ] = foo / sum(foo)
  for (t in 2:n){
    foo = apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ] = foo / sum(foo) 
    }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv)
}



# Specifying likelihood with wild distribution for height.fd --------------

mllk_ = function(theta.star, X, N){
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) # stationary, smaller tolerance to avoid numerical problems
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  p.b = plogis(theta.star[(N-1)*N+4*N+1:N])
  
  # crazy distribution: Height first difference
  mu1 = exp(theta.star[22])
  mu2 = theta.star[23]
  mu3 = theta.star[24]
  mu = c(mu1, mu2, mu3)
  sigma = exp(theta.star[(N-1)*N+6*N+1:N]) # sds of normal distributions
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  
  allprobs[ind,1] =
    dgamma(X$step[ind], shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])* #gamma
    zbeta(X$angle[ind], shape1 = alpha[1], shape2 = beta[1], p = p.b[1])* # beta
    dgamma(X$height.fd[ind], shape = mu[1]^2/sigma[1]^2, scale = sigma[1]^2/mu[1]) # normal
  allprobs[ind,2] =
    dgamma(X$step[ind], shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])* #gamma
    zbeta(X$angle[ind], shape1 = alpha[2], shape2 = beta[2], p = p.b[2])* # beta
    flip_gamma(X$height.fd[ind], mu[2], sigma[2]) # flipped gamma
  allprobs[ind,3] =
    dgamma(X$step[ind], shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])* #gamma
    zbeta(X$angle[ind], shape1 = alpha[3], shape2 = beta[3], p = p.b[3])* # beta
    dnorm(X$height.fd[ind], mu[3], sigma[3]) # gamma 
  
  # forward algorithm to compute the log-likelihood
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo/sum(foo)
  for (t in 2:nrow(X)){
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l+log(sum(foo))
    phi = foo/sum(foo)
  }
  
  return(-l) 
}
