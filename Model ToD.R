theta0.tod = c(Gamma[2:4,1], Gamma[c(1,3:4),2], Gamma[c(1:2,4),3], Gamma[1:3,4], rep(0,24),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.7, 10, 1, 3, # alphas
               2, 55, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5) # al

theta.star0.tod = c(theta0.tod[1:36],
                    log(theta0.tod[37:52]),
                    theta0.tod[53:56],
                    log(theta0.tod[57:60]),
                    theta0.tod[61:64])

t1 = Sys.time()
mod8 = nlm(f = mllk_tod, p = theta.star0.tod, X = data4, N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod8$estimate
states = viterbi_tod(mod8$estimate, X = data2[1:5000,], N = 4)

# AIC
2*mllk_tod(mod8$estimate, X = data2[1:5000,], N = 4) + 2*length(mod8$estimate)
# BIC
2*mllk_tod(mod8$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod8$estimate)


# Results -----------------------------------------------------------------

color = c("deepskyblue", "orange", "springgreen4", "dodgerblue3")

# Scatterplot

data4 = data2[1:5000,]
par(mfrow = c(3,1))
plot(data4$step[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Step lenght")
legend(x = 4700, y = 30, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1)
plot(data4$angle[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Turning angle")
plot(data4$height.fd[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Height (fd.)")


# Hypothetical stationary

todseq = seq(min(data4$time), max(data4$time), length.out = 500)
h_delta = solve_gamma_tod(mod8$estimate, todseq, 4)

par(mfrow = c(1,1))
plot(todseq, h_delta[,1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
lines(todseq, h_delta[,2], lwd = 2, col = color[2])
lines(todseq, h_delta[,3], lwd = 2, col = color[3])
lines(todseq, h_delta[,4], lwd = 2, col = color[4])
legend(15, 0.8, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1, box.lwd = 0)


# Transition probabilities

transprobs_tod = get_transprobs_tod(mod8$estimate, todseq)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(todseq, transprobs_tod[,i], main = NULL, 
       ylab = colnames(transprobs_tod)[i], xlab = "time of day", 
       type = "l", lwd = 2, col = color[st])
}


# Hypothetical marginals

todseq2 = seq(min(data4$time), max(data4$time), length.out = 4)

N = 4
mu.g = exp(theta.star[3*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[3*(N-1)*N+N+1:N]) # sds of gamma distributions
alpha = exp(theta.star[3*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[3*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
xi = theta.star[3*(N-1)*N+4*N+1:N] # means of normal distributions
omega = exp(theta.star[3*(N-1)*N+5*N+1:N]) # sds of normal distributions
al = theta.star[3*(N-1)*N+6*N+1:N]

par(mfrow = c(4,3))
par(mar = c(4, 4, 1.5, 1.5))
for (i in 1:length(todseq2)){
  delta = solve_gamma_tod(mod8$estimate, todseq2[i], N = 4)
  curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), col = color[1], lwd = 1, xlim = c(0,30), ylim = c(0,0.35), ylab = "density", xlab = "step length", n = 300, main = paste("ToD =", round(todseq2[i], 0)))
  curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), col = color[2], add = T, lwd = 1, n = 300)
  curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), col = color[3], add = T, lwd = 1, n = 300)
  curve(delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), col = color[4], add = T, lwd = 1, n = 300)
  curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
          delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
          delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
          delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), lty = "dashed", lwd = 2, add = T, n = 300)
  
  curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), col = color[1], lwd = 1, xlim = c(0,1), ylim = c(0,6.5), ylab = "density", xlab = "turning angle", n = 300, main = paste("ToD =", round(todseq2[i], 0)))
  curve(delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), col = color[2], lwd = 1, add = T, n = 300)
  curve(delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), col = color[3], lwd = 1, add = T, n = 300)
  curve(delta[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), col = color[4], lwd = 1, add = T, n = 300)
  curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
          delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
          delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3])+
          delta[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), lty = "dashed", lwd = 2, add = T, n = 300)
  
  curve(delta[1]*dsn(x, xi = xi[1], omega = omega[1], alpha = al[1]), col = color[1], lwd = 1, xlim = c(-10,10), ylim = c(0,0.25), ylab = "density", xlab = "height.fd", n = 300, main = paste("ToD =", round(todseq2[i], 0)))
  curve(delta[2]*dsn(x, xi = xi[2], omega = omega[2], alpha = al[2]), col = color[2], lwd = 1, add = T, n = 300)
  curve(delta[3]*dsn(x, xi = xi[3], omega = omega[3], alpha = al[3]), col = color[3], lwd = 1, add = T, n = 300)
  curve(delta[4]*dsn(x, xi = xi[4], omega = omega[4], alpha = al[4]), col = color[4], lwd = 1, add = T, n = 300)
  curve(delta[1]*dsn(x, xi = xi[1], omega = omega[1], alpha = al[1])+
          delta[2]*dsn(x, xi = xi[2], omega = omega[2], alpha = al[2])+
          delta[3]*dsn(x, xi = xi[3], omega = omega[3], alpha = al[3])+
          delta[4]*dsn(x, xi = xi[4], omega = omega[4], alpha = al[4]), lty = "dashed", lwd = 2, add = T, n = 300)
}

