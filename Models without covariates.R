
# 3 state model -----------------------------------------------------------

# hist(data4$angle, prob = T, breaks = 50)
# curve(dbeta(x, 10, 55), add = T, lwd = 2)
# curve(dbeta(x, 1, 50), add = T, lwd = 2)
# curve(dbeta(x, 0.5, 2), add = T, lwd = 2)
# 
# hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))
# library(sn)
# curve(dsn(x, 0, 1, 3), add = T, lwd = 2)
# curve(dsn(x, 0, 1, -3), add = T, lwd = 2)
# curve(dsn(x, 0, 0.1, 0), add = T, lwd = 2)

theta0 = c(rep(0.1, 6),
             8, 15, 2, # mu.gamma
             4, 10, 2, # sigma.gamma
             10, 1, 1, # alphas
             55, 40, 10, # betas
             0.2, -0.2, 0.4, # xi
             0.5, 0.5, 0.5, # omega
             5, -5, 0) 

theta.star0 = c(log(theta0[1:18]), theta0[19:21], log(theta0[22:24]), theta0[25:27])

data4 = data2[1:5000,]
data4 = data4[-which(data4$day >= 41 & data4$day <= 48),]

mod1 = nlm(mllk, theta.star0, X = data4, N = 3, print.level = 2, hessian = T, iterlim = 1000)
theta.star = mod1$estimate
sds = sqrt(diag(solve(mod1$hessian)))
N = 3

# setwd("/Users/jan-ole/R/EagleHMM")
# saveRDS(mod1, "mod_no_cov_3st.rds")

plotrix::plotCI(x = 1:length(theta.star), y = theta.star, ui=theta.star+1.96*sds, li=theta.star-1.96*sds, 
       pch = 20, bty = "n", xlab = "index of theta*", ylab = "theta*")
abline(h = 0)


Gamma = diag(N)
Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
Gamma = Gamma/rowSums(Gamma)
delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) 
mu.g = exp(theta.star[(N-1)*N+1:N])
sigma.g = exp(theta.star[(N-1)*N+N+1:N])
alpha = exp(theta.star[(N-1)*N+2*N+1:N])
beta = exp(theta.star[(N-1)*N+3*N+1:N]) 
xi = theta.star[(N-1)*N+4*N+1:N]
omega = exp(theta.star[(N-1)*N+5*N+1:N])
al = theta.star[(N-1)*N+6*N+1:N]


# Plotting results
color = c("orange", "springgreen3", "deepskyblue")
par(mfrow = c(1,1))
hist(data4$step, prob = T, breaks = 50)
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$angle, prob = T, breaks = 50)
curve(delta[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[1]*dbeta(x, alpha[1], beta[1])+
        delta[2]*dbeta(x, alpha[2], beta[2])+
        delta[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))
curve(delta[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[1]*dsn(x, xi[1], omega[1], al[1])+
        delta[2]*dsn(x, xi[2], omega[2], al[2])+
        delta[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, lty = "dashed", n = 500)

states = viterbi(theta.star, data4, N = 3)

# Decoded states

par(mfrow = c(3,1))
plot(data4$step, type = "h", bty = "n", col = color[states], ylab = "step length")
plot(data4$angle, type = "h", bty = "n", col = color[states], ylab = "turning angle")
plot(data4$height.fd, type = "h", bty = "n", col = color[states], ylab = "height difference")



# 4 state model -----------------------------------------------------------


theta0.4 = c(rep(0.1, 12),
              1, 8, 15, 2, # mu.gamma
              1, 4, 10, 2, # sigma.gamma
              0.7, 10, 1, 1, # alphas
              2, 55, 40, 10, # betas
              0, 0.2, -0.2, 0.4, # xi
              0.05, 0.5, 0.5, 0.5, # omega
              0, 5, -5, 5) 

theta.star0.4 = c(log(theta0.4[1:28]),
                   theta0.4[29:32], 
                   log(theta0.4[33:36]),
                   theta0.4[37:40])

mod = nlm(mllk, theta.star0.4, X = data4, N = 4, print.level = 2, hessian = T, iterlim = 1000)
theta.star = mod$estimate
sds = sqrt(diag(solve(mod$hessian)))

# saveRDS(mod, "mod_no_cov_4st.rds")

plotrix::plotCI(x = 1:length(theta.star), y = theta.star, ui=theta.star+1.96*sds, li=theta.star-1.96*sds, 
                pch = 20, bty = "n", xlab = "index of theta*", ylab = "theta*")
abline(h = 0)

N = 4

Gamma = diag(N)
Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
Gamma = Gamma/rowSums(Gamma)
delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-18) 
mu.g = exp(theta.star[(N-1)*N+1:N])
sigma.g = exp(theta.star[(N-1)*N+N+1:N])
alpha = exp(theta.star[(N-1)*N+2*N+1:N])
beta = exp(theta.star[(N-1)*N+3*N+1:N]) 
xi = theta.star[(N-1)*N+4*N+1:N]
omega = exp(theta.star[(N-1)*N+5*N+1:N])
al = theta.star[(N-1)*N+6*N+1:N]

# Plotting results
color = c("deepskyblue", "orange", "springgreen3", "dodgerblue3")

hist(data4$step, prob = T, breaks = 50)
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$angle, prob = T, breaks = 30)
curve(delta[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta[1]*dbeta(x, alpha[1], beta[1])+
        delta[2]*dbeta(x, alpha[2], beta[2])+
        delta[3]*dbeta(x, alpha[3], beta[3])+
        delta[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))
curve(delta[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta[1]*dsn(x, xi[1], omega[1], al[1])+
        delta[2]*dsn(x, xi[2], omega[2], al[2])+
        delta[3]*dsn(x, xi[3], omega[3], al[3])+
        delta[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, lty = "dashed", n = 500)

states = viterbi(theta.star, data4, N = 4)


# Plotting the decoded states ---------------------------------------------

par(mfrow = c(3,1))
plot(data4$step, type = "h", bty = "n", col = color[states], ylab = "step length")
plot(data4$angle, type = "h", bty = "n", col = color[states], ylab = "turning angle")
plot(data4$height.fd, type = "h", bty = "n", col = color[states], ylab = "height difference")


       