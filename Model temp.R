# theta0.5 = c(Gamma[2:4,1], Gamma[c(1,3:4),2], Gamma[c(1:2,4),3], Gamma[1:3,4], rep(0,12),
#              1, 8, 15, 2, # mu.gamma
#              1, 4, 10, 2, # sigma.gamma
#              0.7, 10, 1, 3, # alphas
#              2, 55, 40, 40, # betas
#              0, 0.2, -0.2, 0.4, # xi
#              0.05, 0.5, 0.5, 0.5, # omega
#              0, 5, -5, 5, # al
#              1,1,1) # delta
# 
# theta.star0.5 = c(theta0.5[1:24],
#                   log(theta0.5[25:40]),
#                   theta0.5[41:44],
#                   log(theta0.5[45:48]),
#                   theta0.5[49:55])

# t1 = Sys.time()
# mod6 = nlm(f = mllk_na_sn_cov2, p = theta.star0.5, X = data2[1:5000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8, hessian = T)
# Sys.time()-t1
mod6 = readRDS("mod6.rds")
theta.star = mod6$estimate
states = viterbi_na_sn_cov2(mod6$estimate, X = data2[1:5000,], N = 4)

# AIC
2*mllk_na_sn_cov2(mod6$estimate, X = data2[1:5000,], N = 4) + 2*length(mod6$estimate)
# BIC
2*mllk_na_sn_cov2(mod6$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod6$estimate)


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

tempseq = seq(min(data4$temp), max(data4$temp), length.out = 500)
h_delta = solve_gamma_na_sn_cov(mod6$estimate, tempseq, 4)

par(mfrow = c(1,1))
plot(tempseq, h_delta[,1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
lines(tempseq, h_delta[,2], lwd = 2, col = color[2])
lines(tempseq, h_delta[,3], lwd = 2, col = color[3])
lines(tempseq, h_delta[,4], lwd = 2, col = color[4])
legend(15, 1, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1, box.lwd = 0)


# Transition probabilities

transprobs_temp = get_transprobs(mod6$estimate, tempseq)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(tempseq, transprobs_temp[,i], main = NULL, 
       ylab = colnames(transprobs_temp)[i], xlab = "time of day", 
       type = "l", lwd = 2, col = color[st])
}


# Marginal distribution with relative state frequencies

delta_star = c(sum(states == 1), sum(states == 2), sum(states == 3), sum(states == 4))
delta_star = delta_star/sum(delta_star)

# getting the parameters:
N = 4
coef = matrix(theta.star[1:(2*(N-1)*N)], (N-1)*N, 2)
mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions
alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
al = theta.star[2*(N-1)*N+6*N+1:N]

par(mfrow = c(1,1))

# Step length
hist(data4$step, prob = T, breaks = 100, main = "Histogram of Step length", xlab = "Step length")
curve(delta_star[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)

curve(
  delta_star[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta_star[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta_star[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
    delta_star[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]),
  add = T, lty = "dashed", lwd = 2, n = 500
)

# Turning angle
hist(data4$angle, prob = T, breaks = 50, main = "Histogram of Turning angle", xlab = "Turning angle", xlim = c(0,1))
curve(delta_star[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta_star[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta_star[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(delta_star[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), add = T, lwd = 2, col = color[4], n = 500)

curve(
  delta_star[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
    delta_star[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
    delta_star[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3])+
    delta_star[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]),
  add = T, lty = "dashed", lwd = 2
)
# vllt noch bessere Startwerte für State 2?

# Height (fd.)
hist(data4$height.fd, prob = T, breaks = 190, main = "Histogram of Height (fd.)", xlab = "Height (fd.)", xlim = c(-5,5), ylim = c(0, 1.8))
curve(delta_star[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 1000)
curve(delta_star[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2])
curve(delta_star[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3])
curve(delta_star[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4])

curve(
  delta_star[1]*dsn(x, xi[1], omega[1], al[1])+
    delta_star[2]*dsn(x, xi[2], omega[2], al[2])+
    delta_star[3]*dsn(x, xi[3], omega[3], al[3])+
    delta_star[4]*dsn(x, xi[4], omega[4], al[4]),
  add = T, lty = "dashed", lwd = 2, n = 1000
)

