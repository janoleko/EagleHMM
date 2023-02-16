# days 1:110
data4 = data2[which(data2$day <= 110),]
data4 = data4[-which(data4$day >= 41 & data4$day <= 48),] # exclude weird part

# 3 state -----------------------------------------------------------------

theta0 = c(log(rep(0.1, 6)), rep(0, 24), # intercepts and slopes for 6 etas
           8, 15, 2, # mu.gamma
           4, 10, 2, # sigma.gamma
           10, 1, 1, # alphas
           55, 40, 10, # betas
           0.2, -0.2, 0.4, # xi
           0.5, 0.5, 0.5, # omega
           5, -5, 0,
           0.5, 0.5) # delta

theta.star0 = c(theta0[1:30], log(theta0[31:42]), theta0[43:45], log(theta0[46:48]), theta0[49:51], log(theta0[52:53]))

mod_ttm3 = nlm(mllk_ttm, theta.star0, X = data4, N = 3, print.level = 2, hessian = T, iterlim = 10000)

mod_ttm3 = readRDS("mod_ttm3.rds")
theta.star = mod_ttm3$estimate
sds = sqrt(diag(solve(mod_ttm3$hessian)))

AIC = 2*mod_ttm3$minimum + 2*length(mod_ttm3$estimate)
BIC = 2*mod_ttm3$minimum + log(nrow(data4))*length(mod_ttm3$estimate)

# saveRDS(mod_ttm3, "mod_ttm3.rds")

N = 3
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

states = viterbi_ttm(theta.star, data4, N = 3)
stationary = c(sum(states == 1), sum(states == 2), sum(states == 3))/length(states)

# Plotting results

color = c("orange", "springgreen3", "deepskyblue")
par(mfrow = c(1,1))
hist(data4$step, prob = T, breaks = 50)
curve(stationary[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        stationary[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        stationary[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$angle, prob = T, breaks = 50)
curve(stationary[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[1]*dbeta(x, alpha[1], beta[1])+
        stationary[2]*dbeta(x, alpha[2], beta[2])+
        stationary[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))
curve(stationary[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[1]*dsn(x, xi[1], omega[1], al[1])+
        stationary[2]*dsn(x, xi[2], omega[2], al[2])+
        stationary[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, lty = "dashed", n = 500)

# time series
par(mfrow = c(3,1))
plot(data4$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
plot(data4$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
plot(data4$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference")



# 4 state model -----------------------------------------------------------

theta0.4 = c(log(rep(0.1, 12)), rep(0, 12*4),
             1, 8, 15, 2, # mu.gamma
             1, 4, 10, 2, # sigma.gamma
             0.2, 10, 1, 1, # alphas
             2, 55, 40, 10, # betas
             0, 0.2, -0.2, 0.4, # xi
             0.05, 0.5, 0.5, 0.5, # omega
             0, 5, -5, 5,
             0.5, 0.5, 0.5) # delta

theta.star0.4 = c(theta0.4[1:60], log(theta0.4[61:76]), theta0.4[77:80], log(theta0.4[81:84]), theta0.4[85:88], log(theta0.4[89:91]))

t1 = Sys.time()
mod_ttm4 = nlm(mllk_ttm, theta.star0.4, X = data4, N = 4, print.level = 2, hessian = T, iterlim = 10000)
Sys.time()-t1

mod_ttm4 = readRDS("mod_ttm4.rds")

theta.star = mod_ttm4$estimate
sds = sqrt(diag(solve(mod_ttm4$hessian)))

AIC = 2*mod_ttm4$minimum + 2*length(mod_ttm4$estimate)
BIC = 2*mod_ttm4$minimum + log(nrow(data4))*length(mod_ttm4$estimate)

states = viterbi_ttm(theta.star, X = data4, N = 4)
stationary = c(sum(states == 1), sum(states == 2), sum(states == 3), sum(states == 4))/length(states)

color = c("deepskyblue", "orange", "springgreen3", "dodgerblue3")
par(mfrow = c(3,1))
plot(data4$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
plot(data4$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
plot(data4$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference")


N = 4
mu.g = exp(theta.star[5*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[5*(N-1)*N+N+1:N]) # sds of gamma distributions
alpha = exp(theta.star[5*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[5*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
xi = theta.star[5*(N-1)*N+4*N+1:N]
omega = exp(theta.star[5*(N-1)*N+5*N+1:N])
al = theta.star[5*(N-1)*N+6*N+1:N]

par(mfrow = c(1,1))
hist(data4$step, prob = T, breaks = 50)
curve(stationary[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(stationary[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        stationary[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        stationary[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        stationary[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$angle, prob = T, breaks = 50)
curve(stationary[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(stationary[1]*dbeta(x, alpha[1], beta[1])+
        stationary[2]*dbeta(x, alpha[2], beta[2])+
        stationary[3]*dbeta(x, alpha[3], beta[3])+
        stationary[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, lty = "dashed", n = 500)

hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))
curve(stationary[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(stationary[1]*dsn(x, xi[1], omega[1], al[1])+
        stationary[2]*dsn(x, xi[2], omega[2], al[2])+
        stationary[3]*dsn(x, xi[3], omega[3], al[3])+
        stationary[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, lty = "dashed", n = 500)



# Plot transprobs ---------------------------------------------------------

transprobs_temp = get_transprobs_ttm(theta.star, variable = "temp", X = data4, N = 4)
tempseq = seq(min(data4$temp, na.rm =T),  max(data4$temp, na.rm = T), length.out = 100)

# CIs
#transprobs_temp_CI1 = get_transprobs_ttm(theta.star-1.96*sds, variable = "temp", X = data4, N = 4)
#transprobs_temp_CI2 = get_transprobs_ttm(theta.star+1.96*sds, variable = "temp", X = data4, N = 4)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(tempseq, transprobs_temp[,i], main = NULL, 
       ylab = paste(st, "->", i%%4), xlab = "temperature", 
       type = "l", lwd = 2, col = color[st])
  #lines(tempseq, transprobs_temp_CI1[,i], lwd = 1, col = color[st], linetype = "dashed")
  #lines(tempseq, transprobs_temp_CI2[,i], lwd = 1, col = color[st], linetype = "dashed")
  
}



transprobs_mTPI = get_transprobs_ttm(theta.star, variable = "mTPI", X = data4, N = 4)
mTPIseq = seq(min(data4$mTPI, na.rm =T),  max(data4$mTPI, na.rm = T), length.out = 100)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(mTPIseq, transprobs_mTPI[,i], main = NULL, 
       ylab = paste(st, "->", i%%4), xlab = "mTPI", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}


transprobs_tod = get_transprobs_ttm(theta.star, variable = "time", X = data4, N = 4)
todseq = seq(min(data4$time, na.rm =T),  max(data4$time, na.rm = T), length.out = 100)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(todseq, transprobs_tod[,i], main = NULL, 
       ylab = paste(st, "->", i%%4), xlab = "time of day", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}



###########################################################################
# Estimating model also for return trip -----------------------------------
###########################################################################
# days 187-324

data6 = data2[which(data2$day>=187 & data2$day <= 320),]

data7 = rbind(data4, data6)

# 3 state -----------------------------------------------------------------

theta.star0 = mod_ttm3$estimate

t1 = Sys.time()
mod_ttm3_return = nlm(mllk_ttm, theta.star0, X = data6, N = 3, print.level = 2, hessian = T, iterlim = 10000)
Sys.time()-t1

mod_ttm3_return = readRDS("mod_ttm3.rds")

N = 3
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

states = viterbi_ttm(theta.star, data7, N = 3)
stationary = c(sum(states == 1), sum(states == 2), sum(states == 3))/length(states)

st = 3

dwt = get_dwell_times(states, st)
dwf = get_dwell_freq(dwt)

time = seq(0.25, 24, length.out = 24*4)
plot(ddwell(1:48, state = st, N = 3, coef = coef), pch = 16, col = color[st], ylim = c(0,1),
     ylab = paste0("d",st,"(r)"), xlab = "r", main = paste("Dwell times in state",st), bty = "n")
points(1:max(dwt), dwf, type = "h")
legend(x = 1, y = 1, c("simulated", "analytically derived"), bty = "n", pch = 16, col = c(color[st], "black"))




# 4 state -----------------------------------------------------------------

curve(dsn(x, 0, 0.5, 0), add = T)

theta0.4 = c(mod_ttm4$estimate[1:60],
             0.3, 8, 13, 3, # mu.gamma
             0.5, 4, 10, 2, # sigma.gamma
             0.5, 7, 1, 1, # alphas
             3, 45, 100, 20, # betas
             0, 0, 0, 0, # xi
             0.05, 0.5, 0.5, 0.5, # omega
             0, 5, -5, 0,
             0.1, 0.2, 0.2) # delta

theta.star0.4 = c(theta0.4[1:60], log(theta0.4[61:76]), theta0.4[77:80], log(theta0.4[81:84]), theta0.4[85:88], log(theta0.4[89:91]))
# theta.star0.4 = mod_ttm4$estimate

# Parallelization
library(optimParallel)
cl = makeCluster(8)
setDefaultCluster(cl=cl)

t1 = Sys.time()
mod_bothways = optimParallel(par = mod_ttm4$estimate, fn = mllk_ttm, X = data7, N = 4, hessian = T, control = list(trace = 5, maxit = 10000))
Sys.time()-t1
theta.star = mod_bothways$par


# saveRDS(mod_bothways, "mod_bothways.rds")
mod_bothways = readRDS("mod_bothways.rds")
theta.star = mod_bothways$par

t1 = Sys.time()
mod_ttm4_return = nlm(mllk_ttm, theta.star0.4, X = data6[1:1000,], N = 4, print.level = 2, hessian = T, iterlim = 10000)
Sys.time()-t1
theta.star = mod_ttm4_return$estimate

sds = sqrt(diag(solve(mod_ttm4_return$hessian+diag(rep(1e-8, length(theta.star))))))

AIC = 2*mod_parallel$value + 2*length(theta.star)
BIC = 2*mod_parallel$value + log(nrow(data7))*length(theta.star)


# State decoding
states = viterbi_ttm(theta.star, X = data7, N = 4)
states1 = states[1:4658]
states2 = states[4659:length(states)]

stationary = c(sum(states == 1), sum(states == 2), sum(states == 3), sum(states == 4))/length(states)
stationary1 = c(sum(states1 == 1), sum(states1 == 2), sum(states1 == 3), sum(states1 == 4))/length(states1)
stationary2 = c(sum(states2 == 1), sum(states2 == 2), sum(states2 == 3), sum(states2 == 4))/length(states2)
names(stationary1) = c("state 1", "state 2", "state 3", "state 4")
stat = t(cbind(stationary1, stationary2))


## Plotting the time spent in each state for both trips
par(mfrow = c(1,1))
barplot(stat, beside = T, col = c("deepskyblue", "deepskyblue3", "orange", "orange3", "springgreen3", "springgreen4", "dodgerblue", "dodgerblue3"), 
        border = "white", ylim = c(0,0.7), bty = "n",
        ylab = "Proportion of time in state", main = "Time spent in each state for first and second trip")
legend(9, 0.65, c("first trip", "second trip"), col = c("lightgrey", "darkgrey"), pch = 15, bty = "n")

## different layout
stat = cbind(stationary1, stationary2)
barplot(stat, col = c("deepskyblue", "orange", "springgreen3", "dodgerblue"), xlim = c(0,2), width = 0.5, space = 0.7,
        border = "white", ylim = c(0,1), bty = "n", names.arg = c("first trip", "second trip"),
        ylab = "Proportion of time in state", main = "Time spent in each state for first and second trip")


dev.off()
par(mfrow = c(1,1))
# spent time in states at different times in the day
data7$states = states

states_sub = data7 %>% filter(time < 5) %>% select(states)
stat1 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 5 & time < 6) %>% select(states)
stat2 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 6 & time < 7) %>% select(states)
stat3 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 7 & time < 8) %>% select(states)
stat4 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 8 & time < 9) %>% select(states)
stat5 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 9 & time < 10) %>% select(states)
stat6 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 10 & time < 11) %>% select(states)
stat7 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 11 & time < 12) %>% select(states)
stat8 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 12 & time < 13) %>% select(states)
stat9 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 13 & time < 14) %>% select(states)
stat10 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 14 & time < 15) %>% select(states)
stat11 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 15 & time < 16) %>% select(states)
stat12 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 16 & time < 17) %>% select(states)
stat13 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 17 & time < 18) %>% select(states)
stat14 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 18 & time < 19) %>% select(states)
stat15 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(time >= 19) %>% select(states)
stat16 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
stat_time = cbind(stat1, stat2, stat3, stat4, stat5, stat6, stat7, stat8, stat9, stat10, stat11, stat12, stat13, stat14, stat15, stat16)

barplot(stat_time, col = c("deepskyblue", "orange", "springgreen3", "dodgerblue"), space = 0.5,
        border = "white", ylim = c(0,1), bty = "n", names.arg = c("<5", 5:18, ">=19"),
        ylab = "Proportion of time in state", main = "Time spent in each state for different times", xlab = "Time")
legend(x = 0.23, col = c("deepskyblue", "orange", "springgreen3", "dodgerblue"), legend = c("state 1", "state 2", "state 3", "state 4"), pch = 16, box.col = "white")


## for different temperatures
seq(from = min(data7$temp), to = max(data7$temp), length.out = 5)

states_sub = data7 %>% filter(temp < 23) %>% select(states)
stat1 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(temp >= 23 & temp < 33) %>% select(states)
stat2 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(temp >= 33 & temp < 43) %>% select(states)
stat3 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(temp >= 43) %>% select(states)
stat4 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
stat_time = cbind(stat1, stat2, stat3, stat4)

barplot(stat_time, col = c("deepskyblue", "orange", "springgreen3", "dodgerblue"), space = 0.5,
        border = "white", ylim = c(0,1), bty = "n", names.arg = c("< 23°", "23° to 33°", "33° to 43°", "> 43°"),
        ylab = "Proportion of time in state", main = "Time spent in each state for different temperatures")

## for different mTPI
par(mfrow = c(1,1))

seq(from = min(data7$mTPI), to = max(data7$mTPI), length.out = 6)

states_sub = data7 %>% filter(mTPI < -170) %>% select(states)
stat1 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(mTPI >= -170 & mTPI < -65) %>% select(states)
stat2 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(mTPI >= -65 & mTPI < 35) %>% select(states)
stat3 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(mTPI >= 35 & mTPI < 140) %>% select(states)
stat4 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
states_sub = data7 %>% filter(mTPI >= 140) %>% select(states)
stat5 = c(sum(states_sub == 1), sum(states_sub == 2), sum(states_sub == 3), sum(states_sub == 4))/nrow(states_sub)
stat_time = cbind(stat1, stat2, stat3, stat4, stat5)

barplot(stat_time, col = c("deepskyblue", "orange", "springgreen3", "dodgerblue"), space = 0.4,
        border = "white", ylim = c(0,1), bty = "n", 
        names.arg = c("< -170", "-170° to -65", "-65 to 35", "35 to 140", "> 140"),
        ylab = "Proportion of time in state", main = "Time spent in each state for different mTPI")




# Plotting the time series with decoded states
color = c("deepskyblue", "orange", "springgreen3", "dodgerblue3")
par(mfrow = c(3,1))
plot(data7$step, pch = 20, bty = "n", col = color[states], ylab = "step length")
abline(v = 4658, col = "darkgrey", lty = "dashed", lwd = 1)
plot(data7$angle, pch = 20, bty = "n", col = color[states], ylab = "turning angle")
abline(v = 4658, col = "darkgrey", lty = "dashed", lwd = 1)
plot(data7$height.fd, pch = 20, bty = "n", col = color[states], ylab = "height difference", ylim = c(-7,7))
abline(v = 4658, col = "darkgrey", lty = "dashed", lwd = 1)


N = 4
mu.g = exp(theta.star[5*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[5*(N-1)*N+N+1:N]) # sds of gamma distributions
alpha = exp(theta.star[5*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[5*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
xi = theta.star[5*(N-1)*N+4*N+1:N]
omega = exp(theta.star[5*(N-1)*N+5*N+1:N])
al = theta.star[5*(N-1)*N+6*N+1:N]

par(mfrow = c(1,1))
hist(data7$step, prob = T, breaks = 50)
curve(stationary[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(stationary[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
        stationary[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
        stationary[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
        stationary[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, lty = "dashed", n = 500)
hist(data4$step, prob = T, breaks = 50)

hist(data7$angle, prob = T, breaks = 50)
curve(stationary[1]*dbeta(x, alpha[1], beta[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dbeta(x, alpha[2], beta[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dbeta(x, alpha[3], beta[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(stationary[1]*dbeta(x, alpha[1], beta[1])+
        stationary[2]*dbeta(x, alpha[2], beta[2])+
        stationary[3]*dbeta(x, alpha[3], beta[3])+
        stationary[4]*dbeta(x, alpha[4], beta[4]), add = T, lwd = 2, lty = "dashed", n = 500)
hist(data4$angle, prob = T, breaks = 50)

hist(data7$height.fd, prob = T, breaks = 200, xlim = c(-5,5))
curve(stationary[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(stationary[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(stationary[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3], n = 500)
curve(stationary[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4], n = 500)
curve(stationary[1]*dsn(x, xi[1], omega[1], al[1])+
        stationary[2]*dsn(x, xi[2], omega[2], al[2])+
        stationary[3]*dsn(x, xi[3], omega[3], al[3])+
        stationary[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, lty = "dashed", n = 500)
hist(data4$height.fd, prob = T, breaks = 100, xlim = c(-5,5))




## Transprobs

transprobs_temp = get_transprobs_ttm(theta.star, variable = "temp", X = data7, N = 4, todmean = 9, mTPImean = 150)
tempseq = seq(min(data7$temp, na.rm =T),  max(data7$temp, na.rm = T), length.out = 100)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(tempseq, transprobs_temp[,i], main = NULL, 
       ylab = paste(st, "->", i%%4), xlab = "temperature", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}


transprobs_time = get_transprobs_ttm(theta.star, variable = "time", X = data7, N = 4)
timeseq = seq(min(data7$time, na.rm =T),  max(data7$time, na.rm = T), length.out = 100)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(timeseq, transprobs_time[,i], main = NULL, 
       ylab = paste(st, "->", i%%4), xlab = "time", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}


transprobs_mTPI = get_transprobs_ttm(theta.star, variable = "mTPI", X = data7, N = 4)
mTPIseq = seq(min(data7$mTPI, na.rm =T),  max(data7$mTPI, na.rm = T), length.out = 100)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(mTPIseq, transprobs_mTPI[,i], main = NULL, 
       ylab = paste(st, "->", i%%4), xlab = "mTPI", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
}



# Adding decoded states to ider -------------------------------------------

data = as.data.frame(cbind(ider,d, measurement))

# aggregation with mean
data_aggr = data %>% 
  dplyr::select(step, angle, height = height.above.msl, height.fd, measurement, x, y, day, burst, mTPI, landform, landform.type, elevation = elevation.m, temp = external.temperature, time = solar.time) %>% 
  group_by(measurement) %>% 
  summarise(count_steps = sum(!is.na(step)),
            count_angles = sum(!is.na(angle)), # find intervals where there are barely any values
            count_heights = sum(!is.na(height)),
            step = mean(step, na.rm = T)*1000, # converting to m/s
            angle = abs(mean(angle, na.rm = T))/pi,
            height.fd = mean(diff(height), na.rm = T),
            height = mean(height, na.rm = T),
            x = mean(x, na.rm = T),
            y = mean(y, na.rm = T),
            day = floor(median(day, na.rm = T)),
            burst = median(burst),
            mTPI = mean(mTPI, na.rm = T),
            landform = as.integer(names(which.max(table(landform)))),
            landform.type = as.character(names(which.max(table(na.omit(landform.type))))),
            elevation = mean(elevation, na.rm = T),
            temp = mean(temp, na.rm = T),
            time = mean(time, na.rm = T)) %>% 
  ungroup()


data_aggr$states = NA



index1 = which(data2$day <= 110)
index1 = index1[-which(data4$day >= 41 & data4$day <= 48)]

index2 = which(data2$day>=187 & data2$day <= 320)

index = c(index1, index2)

data_aggr$states[index] = states

data_marielle = data_aggr %>% select(x,y,day, burst, states)
View(data_marielle)

write.csv(data_marielle, file = "decoded_states.csv")
