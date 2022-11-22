
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
