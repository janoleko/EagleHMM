library(sn)


# HMM functions -----------------------------------------------------------

zgamma_scalar = function(x, mu, sigma, p){
  if(x == 0){return(p)}
  else{return((1-p)*dgamma(x, shape = mu^2/sigma^2, scale = sigma^2/mu))}
}
zgamma = function(x, mu, sigma, p){
  return(sapply(x, zgamma_scalar, mu = mu, sigma = sigma, p = p))
}


# Likelihoods -------------------------------------------------------------


mllk = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) # stationary, smaller tolerance to avoid numerical problems
  
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  xi = theta.star[(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      sn::dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      sn::dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  l = numeric(numdays)
  for(i in 1:numdays){
    index = which(X$day == days[i])
    
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}


viterbi = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) # stationary, smaller tolerance to avoid numerical problems
  
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  xi = theta.star[(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  l = numeric(numdays)
  n = nrow(X)
  xi = matrix(0, n, ncol = N)
  
  for(i in 1:numdays){
    index = which(X$day == days[i])
    
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    for (t in 2:length(index)){
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    iv[t] = which.max(Gamma[, iv[t + 1]] * xi[t, ]) 
  }
  return(iv)  
}




# likelihood for 3D-time series assuming contemporaneous conditional independence
# mllk = function(theta.star, X, N){
#   Gamma = diag(N)
#   Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
#   Gamma = Gamma/rowSums(Gamma)
#   delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) # stationary, smaller tolerance to avoid numerical problems
#   
#   # gamma distribution: Step length
#   mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
#   sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
#   
#   # beta distribution: Turning angle/ pi
#   alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
#   beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
#   p.b = plogis(theta.star[(N-1)*N+4*N+1:N])
#   
#   # normal distribution: Height first difference
#   mu = theta.star[(N-1)*N+5*N+1:N] # means of normal distributions
#   sigma = exp(theta.star[(N-1)*N+6*N+1:N]) # sds of normal distributions
#   
#   allprobs = matrix(1, nrow(X), N)
#   ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
#   
#   for (j in 1:N){ # allprobs matrix
#     allprobs[ind,j] =
#       dgamma(X$step[ind], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
#       zbeta(X$angle[ind], shape1 = alpha[j], shape2 = beta[j], p = p.b[j])* # beta
#       dnorm(X$height.fd[ind], mean = mu[j], sd = sigma[j]) # normal
#   }
#   
#   # forward algorithm to compute the log-likelihood
#   foo = delta%*%diag(allprobs[1,])
#   l = log(sum(foo))
#   phi = foo/sum(foo)
#   for (t in 2:nrow(X)){
#     foo = phi%*%Gamma%*%diag(allprobs[t,])
#     l = l+log(sum(foo))
#     phi = foo/sum(foo)
#   }
#   
#   return(-l) 
# }

mllk_na = function(theta.star, X, N){
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
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      zbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j], p = p.b[j])* # beta
      dnorm(X$height.fd[ind1], mean = mu[j], sd = sigma[j]) # normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dnorm(X$height.fd[ind2], mean = mu[j], sd = sigma[j]) # normal
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


mllk_sn = function(theta.star, X, N){
  # avoid zero beta
  # X$angle[which(X$angle == 0)] = runif(length(X$angle[which(X$angle == 0)]), 0, 0.03)
  
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
  Gamma = Gamma/rowSums(Gamma)
  delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-40) # stationary, smaller tolerance to avoid numerical problems
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions

  # normal distribution: Height first difference
  xi = theta.star[(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind,j] =
      dgamma(X$step[ind], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
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


mllk_na_sn = function(theta.star, X, N){
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
  
  # normal distribution: Height first difference
  xi = theta.star[(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  t1 = Sys.time()
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo/sum(foo)
  for (t in 2:nrow(X)){
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l+log(sum(foo))
    phi = foo/sum(foo)
  }
  # Doppelschleife über Tage
  # am Ende die log Likelihoods pro Tag addieren!
  # oder parallelisierbar (mclapply())
  Sys.time()-t1
  return(-l) 
}


mllk_na_sn = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  
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
  
  # normal distribution: Height first difference
  xi = theta.star[(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  for(i in 1:numdays){
    index = which(X$day == i)
    
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}


# mllk_na_sn_cov = function(theta.star, X, N){
#   coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
#   
#   # gamma distribution: Step length
#   mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
#   sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
#   
#   # beta distribution: Turning angle/ pi
#   alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
#   beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
#   
#   # normal distribution: Height first difference
#   xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
#   omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
#   al = theta.star[2*(N-1)*N+6*N+1:N]
#   
#   delta = c(1, exp(theta.star[2*(N-1)*N+7*N+1:(N-1)]))
#   delta = delta/sum(delta)
#   
#   allprobs = matrix(1, nrow(X), N)
#   ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
#   ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
#   
#   for (j in 1:N){ # allprobs matrix
#     allprobs[ind1,j] =
#       dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
#       dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
#       dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
#     
#     allprobs[ind2,j] =
#       dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
#       dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
#   }
#   
#   # forward algorithm to compute the log-likelihood
#   foo = delta%*%diag(allprobs[1,])
#   l = log(sum(foo))
#   phi = foo/sum(foo)
#   for (t in 2:nrow(X)){
#     eta = coef[,1] + coef[,2]*X$temp[t]
#     Gamma = diag(N)
#     Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
#     Gamma = Gamma/rowSums(Gamma)
#     
#     foo = phi%*%Gamma%*%diag(allprobs[t,])
#     l = l+log(sum(foo))
#     phi = foo/sum(foo)
#   }
#   return(-l) 
# }


# mllk_na_sn_cov2 = function(theta.star, X, N){
#   days = unique(X$day)
#   numdays = length(days)
#   
#   coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
#   
#   # gamma distribution: Step length
#   mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
#   sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
#   
#   # beta distribution: Turning angle/ pi
#   alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
#   beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
#   
#   # normal distribution: Height first difference
#   xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
#   omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
#   al = theta.star[2*(N-1)*N+6*N+1:N]
#   
#   delta = c(1, exp(theta.star[2*(N-1)*N+7*N+1:(N-1)]))
#   delta = delta/sum(delta)
#   
#   allprobs = matrix(1, nrow(X), N)
#   ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
#   ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
#   
#   for (j in 1:N){ # allprobs matrix
#     allprobs[ind1,j] =
#       dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
#       dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
#       dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
#     
#     allprobs[ind2,j] =
#       dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
#       dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
#   }
#   
#   # forward algorithm to compute the log-likelihood
#   l = numeric(numdays)
#   # compute the daywise likelihoods seperately and sum up in the end
#   for(i in 1:numdays){
#     index = which(X$day == i)
#     
#     # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
#     foo = delta%*%diag(allprobs[index[1],])
#     l[i] = log(sum(foo))
#     phi = foo/sum(foo)
#     
#     for (t in 2:length(index)){
#       eta = coef[,1] + coef[,2]*X$temp[index[t]]
#       Gamma = diag(N)
#       Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
#       Gamma = Gamma/rowSums(Gamma)
#       
#       foo = phi%*%Gamma%*%diag(allprobs[index[t],])
#       l[i] = l[i]+log(sum(foo))
#       phi = foo/sum(foo)
#     }
#   }
#   l = sum(l)
#   return(-l) 
# }


mllk_tod = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  coef = matrix(theta.star[1:(3*(N-1)*N)], (N-1)*N, 3)
  
  mu.g = exp(theta.star[3*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[3*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  alpha = exp(theta.star[3*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[3*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  xi = theta.star[3*(N-1)*N+4*N+1:N]
  omega = exp(theta.star[3*(N-1)*N+5*N+1:N])
  al = theta.star[3*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[3*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd))
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])*
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])*
      sn::dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j])
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])*
      sn::dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j])
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  for(i in 1:numdays){ # days independent
    index = which(X$day == days[i])
    
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[t]]/24) + coef[,3]*cos(2*pi*X$time[index[t]]/24)
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}


# Likelihood Funktion with:
# - NA handling for turning angle
# - Skew normal distribution for height.fd
# - Seperate daily computation of the likelihood 
# - Covariates: Temperature and Time of Day

mllk_tt = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[4*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[4*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[4*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[4*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # skew normal distribution: Height first difference
  xi = theta.star[4*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[4*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[4*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[4*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  # compute the daywise likelihoods seperately and sum up in the end
  for(i in 1:numdays){
    index = which(X$day == days[i])
    
    eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[1]]/24) + coef[,3]*cos(2*pi*X$time[index[1]]/24) + 
      coef[,4]*X$temp[index[1]]
    
    # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[t]]/24) + coef[,3]*cos(2*pi*X$time[index[t]]/24) +
        coef[,4]*X$temp[index[t]]
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}


mllk_landform = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(6*(N-1)*N)], (N-1)*N, 6)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[6*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[6*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[6*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[6*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # skew normal distribution: Height first difference
  xi = theta.star[6*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[6*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[6*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[6*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  # compute the daywise likelihoods seperately and sum up in the end
  for(i in 1:numdays){
    index = which(X$day == i)
    
    eta = coef[,1] + coef[,2]*X$lower_slope[index[1]] + coef[,3]*X$mountain_divide[index[1]] + 
      coef[,4]*X$peak_ridge[index[1]] + coef[,5]*X$upper_slope[index[1]] + 
      coef[,6]*X$valley[index[1]]
    
    # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$lower_slope[index[t]] + coef[,3]*X$mountain_divide[index[t]] + 
        coef[,4]*X$peak_ridge[index[t]] + coef[,5]*X$upper_slope[index[t]] + 
        coef[,6]*X$valley[index[t]]
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}



mllk_ttl = function(theta.star, X, N){
  numdays = max(X$day)
  
  coef = matrix(theta.star[1:(9*(N-1)*N)], (N-1)*N, 9)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[9*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[9*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[9*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[9*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # skew normal distribution: Height first difference
  xi = theta.star[9*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[9*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[9*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[9*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  # compute the daywise likelihoods seperately and sum up in the end
  for(i in 1:numdays){
    index = which(X$day == i)
    
    eta = coef[,1] + coef[,2]*X$lower_slope[index[1]] + coef[,3]*X$mountain_divide[index[1]] + coef[,4]*X$peak_ridge[index[1]] + coef[,5]*X$upper_slope[index[1]] + coef[,6]*X$valley[index[1]] +
      coef[,7]*X$temp[index[1]] + coef[,8]*sin(2*pi*X$time[index[1]]/24) + coef[,9]*cos(2*pi*X$time[index[1]]/24)
    
    # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$lower_slope[index[t]] + coef[,3]*X$mountain_divide[index[t]] + coef[,4]*X$peak_ridge[index[t]] + coef[,5]*X$upper_slope[index[t]] + coef[,6]*X$valley[index[t]] +
        coef[,7]*X$temp[index[t]] + coef[,8]*sin(2*pi*X$time[index[t]]/24) + coef[,9]*cos(2*pi*X$time[index[t]]/24)
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}


mllk_ttm = function(theta.star, X, N){
  days = unique(X$day)
  numdays = length(days)
  coef = matrix(theta.star[1:(5*(N-1)*N)], (N-1)*N, 5)
  
  mu.g = exp(theta.star[5*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[5*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  alpha = exp(theta.star[5*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[5*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  xi = theta.star[5*(N-1)*N+4*N+1:N]
  omega = exp(theta.star[5*(N-1)*N+5*N+1:N])
  al = theta.star[5*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[5*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd))
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])*
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])*
      sn::dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j])
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])*
      sn::dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j])
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  for(i in 1:numdays){ # days independent
    index = which(X$day == days[i])
    
    foo = delta%*%diag(allprobs[index[1],])
    l[i] = log(sum(foo))
    phi = foo/sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$temp[index[t]] + coef[,3]*sin(2*pi*X$time[index[t]]/24) + coef[,4]*cos(2*pi*X$time[index[t]]/24) + coef[,5]*X$mTPI[index[t]]
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = phi%*%Gamma%*%diag(allprobs[index[t],])
      l[i] = l[i]+log(sum(foo))
      phi = foo/sum(foo)
    }
  }
  l = sum(l)
  return(-l) 
}



# Viterbi functions -------------------------------------------------------

viterbi_na_sn_cov2 = function(theta.star, X, N){
  n = nrow(X)
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # normal distribution: Height first difference
  xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[2*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[2*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  xi = matrix(0, n, ncol = N)
  # compute the daywise likelihoods seperately and sum up in the end
  for(i in 1:numdays){
    index = which(X$day == i)
    
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$temp[index[t]]
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*X$temp[t]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv) 
}


viterbi_tod = function(theta.star, X, N){
  n = nrow(X)
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(3*(N-1)*N)], (N-1)*N, 3)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[3*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[3*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[3*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[3*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # normal distribution: Height first difference
  xi = theta.star[3*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[3*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[3*(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  xi = matrix(0, n, ncol = N)

  for(i in 1:numdays){
    index = which(X$day == i)
    
    eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[1]]/24) + coef[,3]*cos(2*pi*X$time[index[1]]/24)
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18)
    
    # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[t]]/24) + coef[,3]*cos(2*pi*X$time[index[t]]/24)
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*sin(2*pi*X$time[t]/24) + coef[,3]*cos(2*pi*X$time[t]/24)
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] = which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv)  
}


viterbi_tt = function(theta.star, X, N){
  n = nrow(X)
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[4*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[4*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[4*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[4*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # normal distribution: Height first difference
  xi = theta.star[4*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[4*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[4*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[4*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  xi = matrix(0, n, ncol = N)
  
  for(i in 1:numdays){
    index = which(X$day == days[i])
    
    eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[1]]/24) + coef[,3]*cos(2*pi*X$time[index[1]]/24) + 
      coef[,4]*X$temp[index[1]]
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)

    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*sin(2*pi*X$time[index[t]]/24) + coef[,3]*cos(2*pi*X$time[index[t]]/24) + 
        coef[,4]*X$temp[index[t]]
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*sin(2*pi*X$time[t]/24) + coef[,3]*cos(2*pi*X$time[t]/24) + 
      coef[,4]*X$temp[t]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] = which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv) 
}


viterbi_landform = function(theta.star, X, N){
  n = nrow(X)
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(6*(N-1)*N)], (N-1)*N, 6)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[6*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[6*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[6*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[6*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # skew normal distribution: Height first difference
  xi = theta.star[6*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[6*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[6*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[6*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  xi = matrix(0, n, ncol = N)
  
  for(i in 1:numdays){
    index = which(X$day == i)
    
    eta = coef[,1] + coef[,2]*X$lower_slope[index[1]] + coef[,3]*X$mountain_divide[index[1]] + 
      coef[,4]*X$peak_ridge[index[1]] + coef[,5]*X$upper_slope[index[1]] + 
      coef[,6]*X$valley[index[1]]
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)

    # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$lower_slope[index[t]] + coef[,3]*X$mountain_divide[index[t]] + 
        coef[,4]*X$peak_ridge[index[t]] + coef[,5]*X$upper_slope[index[t]] + 
        coef[,6]*X$valley[index[t]]
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*X$lower_slope[t] + coef[,3]*X$mountain_divide[t] + 
      coef[,4]*X$peak_ridge[t] + coef[,5]*X$upper_slope[t] + 
      coef[,6]*X$valley[t]
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] = which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv)
}


viterbi_ttl = function(theta.star, X, N){
  n = nrow(X)
  numdays = max(X$day)
  
  coef = matrix(theta.star[1:(9*(N-1)*N)], (N-1)*N, 9)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[9*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[9*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[9*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[9*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # skew normal distribution: Height first difference
  xi = theta.star[9*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[9*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[9*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[9*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  xi = matrix(0, n, ncol = N)
  
  for(i in 1:numdays){
    index = which(X$day == i)
    
    eta = coef[,1] + coef[,2]*X$lower_slope[index[1]] + coef[,3]*X$mountain_divide[index[1]] + coef[,4]*X$peak_ridge[index[1]] + coef[,5]*X$upper_slope[index[1]] + coef[,6]*X$valley[index[1]] +
      coef[,7]*X$temp[index[1]] + coef[,8]*sin(2*pi*X$time[index[1]]/24) + coef[,9]*cos(2*pi*X$time[index[1]]/24)
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    # vllt statt delta die hypothetische stationary zu Tagesbeginn oder gleichverteilt
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$lower_slope[index[t]] + coef[,3]*X$mountain_divide[index[t]] + coef[,4]*X$peak_ridge[index[t]] + coef[,5]*X$upper_slope[index[t]] + coef[,6]*X$valley[index[t]] +
        coef[,7]*X$temp[index[t]] + coef[,8]*sin(2*pi*X$time[index[t]]/24) + coef[,9]*cos(2*pi*X$time[index[t]]/24)
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*X$lower_slope[t] + coef[,3]*X$mountain_divide[t] + coef[,4]*X$peak_ridge[t] + coef[,5]*X$upper_slope[t] + coef[,6]*X$valley[t] +
      coef[,7]*X$temp[t] + coef[,8]*sin(2*pi*X$time[t]/24) + coef[,9]*cos(2*pi*X$time[t]/24)
    
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] = which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv) 
}


viterbi_na_sn_cov3 = function(theta.star, X, N){
  n = nrow(X)
  days = unique(X$day)
  numdays = length(days)
  
  coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # normal distribution: Height first difference
  xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[2*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[2*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  xi = matrix(0, n, ncol = N)
  # compute the daywise likelihoods seperately and sum up in the end
  for(i in 1:numdays){
    index = which(X$day == i)
    
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$temp[index[t]]
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*X$elevation[t]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv) 
}


viterbi_na_sn_cov = function(theta.star, X, N){
  n = nrow(X)
  
  coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
  
  # gamma distribution: Step length
  mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  # beta distribution: Turning angle/ pi
  alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  # normal distribution: Height first difference
  xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[2*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[2*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
  }

  xi = matrix(0, n, ncol = N) 
  foo = delta * allprobs[1, ]
  xi[1, ] = foo / sum(foo)
  for (t in 2:n){
    eta = coef[,1] + coef[,2]*X$temp[t]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    foo = apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ] = foo / sum(foo) 
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv) 
}


viterbi_ttm = function(theta.star, X, N){
  n = nrow(X)
  days = unique(X$day)
  numdays = length(days)
  coef = matrix(theta.star[1:(5*(N-1)*N)], (N-1)*N, 5)
  
  mu.g = exp(theta.star[5*(N-1)*N+1:N]) # means of gamma distributions
  sigma.g = exp(theta.star[5*(N-1)*N+N+1:N]) # sds of gamma distributions
  
  alpha = exp(theta.star[5*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
  beta = exp(theta.star[5*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
  
  xi = theta.star[5*(N-1)*N+4*N+1:N]
  omega = exp(theta.star[5*(N-1)*N+5*N+1:N])
  al = theta.star[5*(N-1)*N+6*N+1:N]
  
  delta = c(1, exp(theta.star[5*(N-1)*N+7*N+1:(N-1)]))
  delta = delta/sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd))
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])*
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])*
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j])
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])*
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j])
  }
  
  # forward algorithm to compute the log-likelihood
  l = numeric(numdays)
  xi = matrix(0, n, ncol = N)
  
  for(i in 1:numdays){
    index = which(X$day == days[i])
    
    foo = delta * allprobs[index[1], ]
    xi[index[1], ] = foo / sum(foo)
    
    for (t in 2:length(index)){
      eta = coef[,1] + coef[,2]*X$temp[index[t]] + coef[,3]*sin(2*pi*X$time[index[t]]/24) + coef[,4]*cos(2*pi*X$time[index[t]]/24) + coef[,5]*X$mTPI[index[t]]
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      foo = apply(xi[index[t - 1], ] * Gamma, 2, max) * allprobs[index[t], ]
      xi[index[t], ] = foo / sum(foo) 
    }
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    eta = coef[,1] + coef[,2]*X$temp[t] + coef[,3]*sin(2*pi*X$time[t]/24) + coef[,4]*cos(2*pi*X$time[t]/24) + coef[,5]*X$mTPI[t]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    iv[t] = which.max(Gamma[, iv[t + 1]] * xi[t, ]) }
  return(iv) 
}



# Functions to get transprobs and delta -----------------------------------

solve_gamma_na_sn_cov = function(theta.star, temp, N){
  coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
  delta = matrix(data = NA, nrow = length(temp), ncol = N)
  
  for (i in 1:length(temp)){
    eta = coef[,1] + coef[,2]*temp[i]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18)
  }
  return(delta)
}


solve_gamma_tod = function(theta.star, tod, N){
  coef = matrix(theta.star[1:(3*(N-1)*N)], (N-1)*N, 3)
  delta = matrix(data = NA, nrow = length(tod), ncol = N)
  
  for (i in 1:length(tod)){
    eta = coef[,1] + coef[,2]*sin(2*pi*tod[i]/24) + coef[,3]*cos(2*pi*tod[i]/24)
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18)
  }
  return(delta)
}

solve_gamma_tt1 = function(theta.star, tod, tempmean, N){
  coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
  delta = matrix(data = NA, nrow = length(tod), ncol = N)
  
  for (i in 1:length(tod)){
    eta = coef[,1] + coef[,2]*sin(2*pi*tod[i]/24) + coef[,3]*cos(2*pi*tod[i]/24) + coef[,4]*tempmean
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
  }
  return(delta)
}

solve_gamma_tt2 = function(theta.star, temp, todmean, N){
  coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
  delta = matrix(data = NA, nrow = length(temp), ncol = N)
  
  for (i in 1:length(temp)){
    eta = coef[,1] + coef[,2]*sin(2*pi*todmean/24) + coef[,3]*cos(2*pi*todmean/24) + coef[,4]*temp[i]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
  }
  return(delta)
}


solve_gamma_landform = function(theta.star, N){
  landform = c("C", "LS", "MD", "PR", "US", "V")
  coef = matrix(theta.star[1:(6*(N-1)*N)], (N-1)*N, 6)
  delta = matrix(data = NA, nrow = length(landform), ncol = N)
  
  eta = coef[,1] 
  Gamma = diag(N)
  Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
  Gamma = Gamma/rowSums(Gamma)
  
  delta[1,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18)
  
  for (i in 2:length(landform)){
    eta = coef[,1] + coef[,i]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18)
  }
  return(list(landform = as.factor(landform), delta = delta))
}


solve_gamma_ttl1 = function(theta.star, tod, tempmean, N){
  coef = matrix(theta.star[1:(9*(N-1)*N)], (N-1)*N, 9)
  delta = matrix(data = NA, nrow = length(tod), ncol = N*6)
  
  for (i in 1:length(tod)){
    eta = coef[,1] + 
      coef[,7]*tempmean + coef[,8]*sin(2*pi*tod[i]/24) + coef[,9]*cos(2*pi*tod[i]/24)
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,1:4] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
    
    for (j in 2:6){
      eta = coef[,1] + coef[,j] +
        coef[,7]*tempmean + coef[,8]*sin(2*pi*tod[i]/24) + coef[,9]*cos(2*pi*tod[i]/24)
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      delta[i,(j-1)*4+1:4] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
    }
  }
  return(delta)
}

solve_gamma_ttl2 = function(theta.star, temp, todmean, N){
  coef = matrix(theta.star[1:(9*(N-1)*N)], (N-1)*N, 9)
  delta = matrix(data = NA, nrow = length(temp), ncol = N*6)
  
  for (i in 1:length(temp)){
    eta = coef[,1] + 
      coef[,7]*temp[i] + coef[,8]*sin(2*pi*todmean/24) + coef[,9]*cos(2*pi*todmean/24)
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,1:4] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
    
    for (j in 2:6){
      eta = coef[,1] + coef[,j] +
      coef[,7]*temp[i] + coef[,8]*sin(2*pi*todmean/24) + coef[,9]*cos(2*pi*todmean/24)
      
      Gamma = diag(N)
      Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
      Gamma = Gamma/rowSums(Gamma)
      
      delta[i,(j-1)*4+1:4] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
    }
  }
  return(delta)
}

solve_gamma_mTPI = function(theta.star, mTPI, tempmean, todmean, N){
  coef = matrix(theta.star[1:(5*(N-1)*N)], (N-1)*N, 5)
  delta = matrix(data = NA, nrow = length(mTPI), ncol = N)
  
  for (i in 1:length(mTPI)){
    eta = coef[,1] + coef[,2]*tempmean + coef[,3]*sin(2*pi*todmean/24) + coef[,4]*cos(2*pi*todmean/24) + coef[,5]*mTPI[i]
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta)
    Gamma = Gamma/rowSums(Gamma)
    
    delta[i,] = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-30)
  }
  return(delta)
}


get_transprobs = function(theta.star, temp){
  N = 4
  coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
  transprobs = matrix(data = NA, nrow = length(temp), ncol = N^2)
  
  for (i in 1:length(temp)){
    eta = coef[,1] + coef[,2]*temp[i]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    transprobs[i,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
  }
  transprobs = as.data.frame(transprobs)
  colnames(transprobs) = c("1 -> 1","1 -> 2", "1 -> 3", "1 -> 4",
                          "2 -> 1", "2 -> 2", "2 -> 3", "2 -> 4",
                          "3 -> 1", "3 -> 2", "3 -> 3", "3 -> 4",
                          "4 -> 1", "4 -> 2", "4 -> 3", "4 -> 4")
  return(transprobs)
}


get_transprobs_tod = function(theta.star, tod){
  N = 4
  coef = matrix(theta.star[1:(3*(N-1)*N)], (N-1)*N, 3)
  transprobs = matrix(data = NA, nrow = length(tod), ncol = N^2)
  
  for (i in 1:length(tod)){
    eta = coef[,1] + coef[,2]*sin(2*pi*tod[i]/24) + coef[,3]*cos(2*pi*tod[i]/24)
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    transprobs[i,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
  }
  transprobs = as.data.frame(transprobs)
  colnames(transprobs) = c("1 -> 1","1 -> 2", "1 -> 3", "1 -> 4",
                           "2 -> 1", "2 -> 2", "2 -> 3", "2 -> 4",
                           "3 -> 1", "3 -> 2", "3 -> 3", "3 -> 4",
                           "4 -> 1", "4 -> 2", "4 -> 3", "4 -> 4")
  return(transprobs)
}

get_transprobs_tt1 = function(theta.star, tod, tempmean){
  N = 4
  coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
  transprobs = matrix(data = NA, nrow = length(tod), ncol = N^2)
  
  for (i in 1:length(tod)){
    eta = coef[,1] + coef[,2]*sin(2*pi*tod[i]/24) + coef[,3]*cos(2*pi*tod[i]/24) + coef[,4]*tempmean
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    transprobs[i,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
  }
  transprobs = as.data.frame(transprobs)
  colnames(transprobs) = c("1 -> 1","1 -> 2", "1 -> 3", "1 -> 4",
                           "2 -> 1", "2 -> 2", "2 -> 3", "2 -> 4",
                           "3 -> 1", "3 -> 2", "3 -> 3", "3 -> 4",
                           "4 -> 1", "4 -> 2", "4 -> 3", "4 -> 4")
  return(transprobs)
}

get_transprobs_tt2 = function(theta.star, temp, todmean){
  N = 4
  coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
  transprobs = matrix(data = NA, nrow = length(temp), ncol = N^2)
  
  for (i in 1:length(temp)){
    eta = coef[,1] + coef[,2]*sin(2*pi*todmean/24) + coef[,3]*cos(2*pi*todmean/24) + coef[,4]*temp[i]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    transprobs[i,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
  }
  transprobs = as.data.frame(transprobs)
  colnames(transprobs) = c("1 -> 1","1 -> 2", "1 -> 3", "1 -> 4",
                           "2 -> 1", "2 -> 2", "2 -> 3", "2 -> 4",
                           "3 -> 1", "3 -> 2", "3 -> 3", "3 -> 4",
                           "4 -> 1", "4 -> 2", "4 -> 3", "4 -> 4")
  return(transprobs)
}


get_transprobs_landform = function(theta.star){
  N = 4
  landform = c("C", "LS", "MD", "PR", "US", "V")
  coef = matrix(theta.star[1:(6*(N-1)*N)], (N-1)*N, 6)
  transprobs = matrix(data = NA, nrow = length(landform), ncol = N^2)
  
  eta = coef[,1]
  Gamma = diag(N)
  Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
  Gamma = Gamma/rowSums(Gamma)
  transprobs[1,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
  
  for (i in 2:length(landform)){
    eta = coef[,1] + coef[,i]
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    transprobs[i,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
  }
  transprobs = as.data.frame(transprobs)
  colnames(transprobs) = c("1 -> 1","1 -> 2", "1 -> 3", "1 -> 4",
                           "2 -> 1", "2 -> 2", "2 -> 3", "2 -> 4",
                           "3 -> 1", "3 -> 2", "3 -> 3", "3 -> 4",
                           "4 -> 1", "4 -> 2", "4 -> 3", "4 -> 4")
  return(list(landform = as.factor(landform), transprobs = transprobs))
}


get_transprobs_ttl1 = function(theta.star, tod, tempmean){
  N = 4
  coef = matrix(theta.star[1:(9*(N-1)*N)], (N-1)*N, 9)
  transprobs = matrix(data = NA, nrow = 6*length(tod), ncol = N^2)
  
  for (i in 1:length(tod)){
    eta = coef[,1] + 
      coef[,7]*tempmean + coef[,8]*sin(2*pi*tod[i]/24) + coef[,9]*cos(2*pi*tod[i]/24)    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
    Gamma = Gamma/rowSums(Gamma)
    
    transprobs[(i-1)*6+1,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
    
    for (j in 2:6){
      eta = coef[,1] + coef[,j] +
        coef[,7]*tempmean + coef[,8]*sin(2*pi*tod[i]/24) + coef[,9]*cos(2*pi*tod[i]/24)
        Gamma = diag(N)
        Gamma[!Gamma] = exp(eta) # dynamically changing Gamma-Matrix
        Gamma = Gamma/rowSums(Gamma)
        
        transprobs[(i-1)*6+j,] = c(Gamma[1,], Gamma[2,], Gamma[3,], Gamma[4,])
    }
  }
  transprobs = as.data.frame(transprobs)
  colnames(transprobs) = c("1 -> 1","1 -> 2", "1 -> 3", "1 -> 4",
                           "2 -> 1", "2 -> 2", "2 -> 3", "2 -> 4",
                           "3 -> 1", "3 -> 2", "3 -> 3", "3 -> 4",
                           "4 -> 1", "4 -> 2", "4 -> 3", "4 -> 4")
  return(transprobs)
}



get_transprobs_ttm = function(theta.star, variable = "temp", X, N = 4, tempmean = mean(X$temp, na.rm = T),
                              todmean = 12, mTPImean = mean(X$mTPI, na.rm = T)){

  sequence = seq(min(X[variable], na.rm = T), max(X[variable], na.rm = T), length.out = 100)
  
  coef = matrix(theta.star[1:(5*(N-1)*N)], (N-1)*N, 5)
  transprobs = matrix(data = NA, nrow = length(sequence), ncol = N^2)
  
  for (i in 1:length(sequence)){
    if(variable == "temp"){
      eta = coef[,1] + coef[,2]*sequence[i] + coef[,3]*sin(2*pi*todmean/24) + coef[,4]*cos(2*pi*todmean/24) + coef[,5]*mTPImean
    }
    if(variable == "time"){
      eta = coef[,1] + coef[,2]*tempmean + coef[,3]*sin(2*pi*sequence[i]/24) + coef[,4]*cos(2*pi*sequence[i]/24) + coef[,5]*mTPImean
    }
    if(variable == "mTPI"){
      eta = coef[,1] + coef[,2]*tempmean + coef[,3]*sin(2*pi*todmean/24) + coef[,4]*cos(2*pi*todmean/24) + coef[,5]*sequence[i]
    }
    
    Gamma = diag(N)
    Gamma[!Gamma] = exp(eta)
    Gamma = Gamma/rowSums(Gamma)
    
    for (j in 1:N){
      transprobs[i,(j-1)*N+1:N] = Gamma[j,]
    }
  }
  
  names = rep(1,N)
  for (j in 2:N){
    names = c(names, rep(j,N))
  }
  colnames(transprobs) = names
  
  return(transprobs)
}






viterbi_na_sn = function(theta.star, X, N){
  n = nrow(X)
  
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
  
  # normal distribution: Height first difference
  xi = theta.star[(N-1)*N+4*N+1:N] # means of normal distributions
  omega = exp(theta.star[(N-1)*N+5*N+1:N]) # sds of normal distributions
  al = theta.star[(N-1)*N+6*N+1:N]
  
  allprobs = matrix(1, nrow(X), N)
  ind1 = which(!is.na(X$step) & !is.na(X$angle) & !is.na(X$height.fd)) # Fragen wie wir das machen --> sonst immer gar keine Information nur weil in angle Zeitreihe so viele NAs
  ind2 = which(!is.na(X$step) & is.na(X$angle) & !is.na(X$height.fd))
  
  for (j in 1:N){ # allprobs matrix
    allprobs[ind1,j] =
      dgamma(X$step[ind1], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dbeta(X$angle[ind1], shape1 = alpha[j], shape2 = beta[j])* # beta
      dsn(X$height.fd[ind1], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
    
    allprobs[ind2,j] =
      dgamma(X$step[ind2], shape = mu.g[j]^2/sigma.g[j]^2, scale = sigma.g[j]^2/mu.g[j])* #gamma
      dsn(X$height.fd[ind2], xi = xi[j], omega = omega[j], alpha = al[j]) # skew normal
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
