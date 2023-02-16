
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
# noch nicht data4 = data4[-which(data4$day >= 41 & data4$day <= 48),]

mod1 = nlm(mllk, theta.star0, X = data4, N = 3, print.level = 2, hessian = T, iterlim = 1000)

mod1 = readRDS("mod_no_cov_3st.rds")

theta.star = mod1$estimate
sds = sqrt(diag(solve(mod1$hessian)))


(AIC3 = 2*mod1$minimum + 2*length(mod1$estimate))
(BIC3 = 2*mod1$minimum + log(nrow(data4))*length(mod1$estimate))


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

par(mfrow = c(1,3))

color = c("orange", "springgreen3", "deepskyblue")
par(mfrow = c(1,3))
hist(data4$step, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "step length", ylim = c(0,0.7))
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, lty = "dashed", n = 500)
legend("topright", paste("state", c(1,2,3)), col = c("deepskyblue", "orange", "springgreen3"), lwd = 2, bty = "n")

hist(data4$angle, prob = T, breaks = 50, border = "white", col = "grey70", main = NULL, xlab = "turning angle", ylim = c(0,12))
curve(delta[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[1]*dbeta(x, alpha[1], beta[1])+
        delta[2]*dbeta(x, alpha[2], beta[2])+
        delta[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5), border = "white", col = "grey70", main = NULL, xlab = "height difference", ylim = c(0,1))
curve(delta[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[1]*dsn(x, xi[1], omega[1], al[1])+
        delta[2]*dsn(x, xi[2], omega[2], al[2])+
        delta[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, lty = "dashed", n = 500)

states = viterbi(theta.star, data4, N = 3)

# Decoded states

par(mfrow = c(3,1))
plot(data4$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
plot(data4$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
plot(data4$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference")



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

# parallel optimization
library(optimParallel)
cl = makeCluster(8); setDefaultCluster(cl=cl)
t1 = Sys.time()
mod = optimParallel(par = theta.star0.4, fn = mllk, X = data4, N = 4, hessian = T, control = list(trace = 5, maxit = 10000))
Sys.time()-t1
stopCluster(cl)

saveRDS(mod, "mod_no_cov_4_st_new.rds")

setwd("/Users/jan-ole/R/EagleHMM")
mod = readRDS("mod_no_cov_4st.rds")

theta.star = mod$par
sds = sqrt(diag(solve(mod$hessian)))

# saveRDS(mod, "mod_no_cov_4st.rds")

par(mfrow = c(1,1))
plotrix::plotCI(x = 1:length(theta.star), y = theta.star, ui=theta.star+1.96*sds, li=theta.star-1.96*sds, 
                pch = 20, bty = "n", xlab = "index of theta*", ylab = "theta*")
abline(h = 0)

(AIC = 2*mod$value + 2*length(theta.star))
(BIC = 2*mod$value + log(nrow(data4))*length(theta.star))


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

par(mfrow = c(1,3))
hist(data4$step, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "step length", ylim = c(0,0.7))
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed", n = 500)
legend("topright", paste("state", c(1,2,3,4)), col = color, lwd = 2, bty = "n")


hist(data4$angle, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "turning angle", ylim = c(0,12))
curve(delta[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta[1]*dbeta(x, alpha[1], beta[1])+
        delta[2]*dbeta(x, alpha[2], beta[2])+
        delta[3]*dbeta(x, alpha[3], beta[3])+
        delta[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5), border = "white", main = NULL, col = "grey70", xlab = "height difference", ylim = c(0,1))
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
plot(data4$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
plot(data4$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
plot(data4$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference")


       