library(optimParallel)

# days 1:110
data4 = data2[which(data2$day <= 110),]
data4 = data4[-which(data4$day >= 41 & data4$day <= 48),] # exclude weird part

# days 187-324
data6 = data2[which(data2$day>=187 & data2$day <= 320),]
data7 = rbind(data4, data6)



# Model only including ToD ------------------------------------------------

setwd("/Users/jan-ole/R/EagleHMM")
mod = readRDS("mod_no_cov_4_st_new.rds")

theta.star = mod$par
N = 4
Gamma = diag(N)
Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
Gamma = Gamma/rowSums(Gamma)

theta0.tod = c(Gamma[-1,1], Gamma[-2,2], Gamma[-3,3], Gamma[-4,4], rep(0,24),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.7, 10, 1, 3, # alphas
               2, 55, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5, # al
               -2, -2, -2) # delta

theta.star0.tod = c(theta0.tod[1:36],
                    log(theta0.tod[37:52]),
                    theta0.tod[53:56],
                    log(theta0.tod[57:60]),
                    theta0.tod[61:67])

# Parallelization

# cl = makeCluster(8); setDefaultCluster(cl=cl)
# t1 = Sys.time()
# mod_tod_bothways = optimParallel(par = theta.star0.tod, fn = mllk_tod, X = data7, N = 4, hessian = T, control = list(trace = 5, maxit = 10000))
# Sys.time()-t1
# stopCluster(cl)

# saveRDS(mod_tod_bothways, "mod_tod_bothways.rds")

theta.star = mod_tod_bothways$par
N = 4

(AIC = 2*mod_tod_bothways$value + 2*length(theta.star))
(BIC = 2*mod_tod_bothways$value + log(nrow(data7))*length(theta.star))

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

color = c("deepskyblue", "orange", "springgreen3", "dodgerblue3")

states = viterbi_tod(theta.star, data7, N = 4)
delta_star = c(sum(states == 1), sum(states == 2), sum(states == 3), sum(states == 4))/length(states)

# Plotting results
color = c("deepskyblue", "orange", "springgreen3", "dodgerblue3")

par(mfrow = c(1,3))
hist(data7$step, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "step length", ylim = c(0,0.7))
curve(delta_star[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_star[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta_star[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta_star[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        delta_star[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed", n = 500)
legend("topright", paste("state", c(1,2,3,4)), col = color, lwd = 2, bty = "n")


hist(data7$angle, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "turning angle", ylim = c(0,12))
curve(delta_star[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_star[1]*dbeta(x, alpha[1], beta[1])+
        delta_star[2]*dbeta(x, alpha[2], beta[2])+
        delta_star[3]*dbeta(x, alpha[3], beta[3])+
        delta_star[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data7$height.fd, prob = T, breaks = 200, xlim = c(-5,5), border = "white", main = NULL, col = "grey70", xlab = "height difference", ylim = c(0,1))
curve(delta_star[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_star[1]*dsn(x, xi[1], omega[1], al[1])+
        delta_star[2]*dsn(x, xi[2], omega[2], al[2])+
        delta_star[3]*dsn(x, xi[3], omega[3], al[3])+
        delta_star[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, lty = "dashed", n = 500)


par(mfrow = c(3,1))
plot(data7$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
plot(data7$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
plot(data7$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference")


## Dwell times

par(mfrow = c(2,2))
for (st in 1:4){
  dwt = get_dwell_times(states, st)
  dwf = get_dwell_freq(dwt)
  
  plot(ddwell(1:(24*4), state = st, N = 4, coef = coef), pch = 16, col = color[st], ylim = c(0,1), xlim = c(1,40),
       ylab = paste0("d",st,"(r)"), xlab = "r", main = paste("Dwell times in state",st), bty = "n")
  points(1:max(dwt), dwf, type = "h")
  legend("topright", c("implied dwell times", "decoded"), bty = "n", 
         pch = c(16,NA), lty = c(NA,1), col = c(color[st], "black"))
}


## Tranition probs
todseq = seq(min(data7$time), max(data7$time), length.out = 500)
transprobs_tod = get_transprobs_tod(theta.star, todseq)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(todseq, transprobs_tod[,i], main = NULL, 
       ylab = colnames(transprobs_tod)[i], xlab = "time of day", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}




# Model with ToD and temp -------------------------------------------------

theta0.tt = c(Gamma[-1,1], Gamma[-2,2], Gamma[-3,3], Gamma[-4,4], rep(0,36),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.7, 10, 1, 3, # alphas
               2, 55, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5, # al
               -2, -2, -2) # delta

theta.star0.tt = c(theta0.tt[1:48],
                    log(theta0.tt[49:64]),
                    theta0.tt[65:68],
                    log(theta0.tt[69:72]),
                    theta0.tt[73:79])

cl = makeCluster(8); setDefaultCluster(cl=cl)
t1 = Sys.time()
mod_tt_bothways = optimParallel(par = theta.star0.tt, fn = mllk_tt, X = data7, N = 4, hessian = T, control = list(trace = 5, maxit = 10000))
Sys.time()-t1
stopCluster(cl)

saveRDS(mod_tt_bothways, "mod_tt_bothways.rds")

theta.star = mod_tt_bothways$par

(AIC = 2*mod_tt_bothways$value + 2*length(theta.star))
(BIC = 2*mod_tt_bothways$value + log(nrow(data7))*length(theta.star))

N = 4
coef = matrix(theta.star[1:(4*(N-1)*N)], (N-1)*N, 4)
mu.g = exp(theta.star[4*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[4*(N-1)*N+N+1:N]) # sds of gamma distributions
alpha = exp(theta.star[4*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[4*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
xi = theta.star[4*(N-1)*N+4*N+1:N]
omega = exp(theta.star[4*(N-1)*N+5*N+1:N])
al = theta.star[4*(N-1)*N+6*N+1:N]
delta = c(1, exp(theta.star[4*(N-1)*N+7*N+1:(N-1)]))
delta = delta/sum(delta)

states = viterbi_tt(theta.star, data7, N = 4)
delta_star = c(sum(states == 1), sum(states == 2), sum(states == 3), sum(states == 4))/length(states)


par(mfrow = c(1,3))
hist(data7$step, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "step length", ylim = c(0,0.7))
curve(delta_star[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_star[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta_star[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta_star[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        delta_star[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed", n = 500)
legend("topright", paste("state", c(1,2,3,4)), col = color, lwd = 2, bty = "n")


hist(data7$angle, prob = T, breaks = 50, border = "white", main = NULL, col = "grey70", xlab = "turning angle", ylim = c(0,12))
curve(delta_star[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_star[1]*dbeta(x, alpha[1], beta[1])+
        delta_star[2]*dbeta(x, alpha[2], beta[2])+
        delta_star[3]*dbeta(x, alpha[3], beta[3])+
        delta_star[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data7$height.fd, prob = T, breaks = 200, xlim = c(-5,5), border = "white", main = NULL, col = "grey70", xlab = "height difference", ylim = c(0,1))
curve(delta_star[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_star[1]*dsn(x, xi[1], omega[1], al[1])+
        delta_star[2]*dsn(x, xi[2], omega[2], al[2])+
        delta_star[3]*dsn(x, xi[3], omega[3], al[3])+
        delta_star[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, lty = "dashed", n = 500)

# state dependent distributions unter dem Modell sehen echt kacke aus

par(mfrow = c(3,1))
plot(data7$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
plot(data7$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
plot(data7$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference")

