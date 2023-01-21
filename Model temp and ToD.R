theta0.tt = c(rep(log(0.05), 12), rep(0,36),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.3, 9, 1, 3, # alphas
               2, 50, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5, # al
              0, 0, 0) # delta

theta.star0.tt = c(theta0.tt[1:48],
                    log(theta0.tt[49:64]),
                    theta0.tt[65:68],
                    log(theta0.tt[69:72]),
                    theta0.tt[73:79])

data4 = data2[1:5000,]
data4 = data4[-which(data4$day >= 41 & data4$day <= 48),]

mllk_tt(theta.star0.tt, X = data4, N = 4)

t1 = Sys.time()
mod9 = nlm(f = mllk_tt, p = theta.star0.tt, X = data4, N = 4, print.level = 2, iterlim = 2000, steptol = 1e-15)
Sys.time()-t1


#mod9 = readRDS("mod9.rds")
theta.star = mod9$estimate
states = viterbi_tt(theta.star, X = data4, N = 4)

# AIC
2*mod9$minimum + 2*length(mod9$estimate)
# BIC
2*mod9$minimum + log(nrow(data4))*length(mod9$estimate)
# BIC better than model only with ToD

# Results -----------------------------------------------------------------

color = c("deepskyblue", "orange", "springgreen4", "dodgerblue3")

# Scatterplot

par(mfrow = c(3,1))
plot(data4$step, pch = 20, col = color[states], ylab = "Step length")
legend(x = 4200, y = 30, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1)
plot(data4$angle, pch = 20, col = color[states], ylab = "Turning angle")
plot(data4$height.fd, pch = 20, col = color[states], ylab = "Height (fd.)")

par(mfrow = c(1,1))

# weirder Zeitraum ist 2000-2350

data4$day[2000:2350]
# Tage 41-48 ausschließen

# Hypothetical stationary

todseq = seq(min(data4$time), max(data4$time), length.out = 500)
tempseq = seq(min(data4$temp), max(data4$temp), length.out = 500)

h_delta1 = solve_gamma_tt1(mod9$estimate, todseq, mean(data4$temp), 4)
h_delta2 = solve_gamma_tt2(mod9$estimate, tempseq, mean(data4$time), 4)

# Plot for time of day
dev.off()
par(mfrow = c(1,2))
plot(todseq, h_delta1[,1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
lines(todseq, h_delta1[,2], lwd = 2, col = color[2])
lines(todseq, h_delta1[,3], lwd = 2, col = color[3])
lines(todseq, h_delta1[,4], lwd = 2, col = color[4])
legend(8.5, 1, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1, box.lwd = 0)

# Plot for temperature
plot(tempseq, h_delta2[,1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "temperature")
lines(tempseq, h_delta2[,2], lwd = 2, col = color[2])
lines(tempseq, h_delta2[,3], lwd = 2, col = color[3])
lines(tempseq, h_delta2[,4], lwd = 2, col = color[4])
legend(15, 1, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1, box.lwd = 0)


# Transition probabilities

transprobs_tt1 = get_transprobs_tt1(mod9$estimate, todseq, mean(data4$temp))
transprobs_tt2 = get_transprobs_tt2(mod9$estimate, tempseq, mean(data4$time))

# Plot for time of day
par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(todseq, transprobs_tt1[,i], main = NULL, 
       ylab = colnames(transprobs_tt1)[i], xlab = "time of day", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}

# Plot for temperature
par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(tempseq, transprobs_tt2[,i], main = NULL, 
       ylab = colnames(transprobs_tt2)[i], xlab = "temperature", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}


# Plot marginals ----------------------------------------------------------
N = 4

mu.g = exp(theta.star[4*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[4*(N-1)*N+N+1:N]) # sds of gamma distributions
# beta distribution: Turning angle/ pi
alpha = exp(theta.star[4*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[4*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
# skew normal distribution: Height first difference
xi = theta.star[4*(N-1)*N+4*N+1:N] # means of normal distributions
omega = exp(theta.star[4*(N-1)*N+5*N+1:N]) # sds of normal distributions
al = theta.star[4*(N-1)*N+6*N+1:N]

delta_hat = c(sum(states==1), sum(states==2), sum(states==3), sum(states==4))/length(states)

par(mfrow = c(1,1))
hist(data4$step, prob = T, breaks = 50)
curve(delta_hat[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_hat[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_hat[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_hat[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(delta_hat[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        delta_hat[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        delta_hat[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        delta_hat[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed")

hist(data4$angle, prob = T, breaks = 100)
curve(delta_hat[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), add = T, col = color[1], lwd = 2)
curve(delta_hat[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), add = T, col = color[2], lwd = 2)
curve(delta_hat[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), add = T, col = color[3], lwd = 2)
curve(delta_hat[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), add = T, col = color[4], lwd = 2)
curve(delta_hat[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
        delta_hat[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
        delta_hat[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3])+
        delta_hat[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), add = T, lwd = 2, lty = "dashed")

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))
curve(delta_hat[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1])
curve(delta_hat[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2])
curve(delta_hat[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3])
curve(delta_hat[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4])
curve(delta_hat[1]*dsn(x, xi[1], omega[1], al[1])+
        delta_hat[2]*dsn(x, xi[2], omega[2], al[2])+
        delta_hat[3]*dsn(x, xi[3], omega[3], al[3])+
        delta_hat[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, lty = "dashed")



hist(data4$angle, prob = T, breaks = 100)
curve(delta_hat[1]*dbeta(x, shape1 = 0.3, shape2 = 2), add = T, col = color[2], lwd = 2)

