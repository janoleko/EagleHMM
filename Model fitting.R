library(moveHMM)
library(tidyverse)

# Load in the data --------------------------------------------------------
setwd("/Users/jan-ole/R/HMM Project")
ider = read.csv("Ider_annotated.csv")
colnames(ider)

ider$height.fd = c(NA, diff(ider$height.above.msl))

N = nrow(ider)/30

measurement = rep(NA, nrow(ider))
for (i in 1:N){
  measurement[(i-1)*30+1:30] = rep(i, 30)
}

d = ider[,5:6]
d$ID = "ID"
colnames(d) = c("x", "y", "ID")
d = prepData(d)


# deleting wrong step lengths and angles
for (i in 1:N){
  d$step[30*i] = NA
}
for (i in 1:(N-1)){
  d$angle[30*i+1] = NA
}

data = as.data.frame(cbind(ider,d, measurement))

# aggregation with mean
data2 = data %>% 
  dplyr::select(step, angle, height = height.above.msl, height.fd, measurement, x, y) %>% 
  group_by(measurement) %>% 
  summarise(step = mean(step, na.rm = T)*1000, # converting to m/s
            count_angles = sum(!is.na(angle)), # find intervals where there are barely any values
            angle = abs(mean(angle, na.rm = T))/pi,
            height.fd = mean(diff(height), na.rm = T),
            x = mean(x, na.rm = T),
            y = mean(y, na.rm = T))

data2$angle[which(data2$count_angles <= 2)] = NA # set interval with barely any values to NA

# Fitting a model with turning angle --------------------------------------

theta0 = c(rep(0.1, 6),
           9, 13, 0.5, # mu.gamma
           1.5, 4, 0.4, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 2, # betas
           0.01, 0.01, 0.01, # zero masses
           1, -0.5, 0, # mu
           1, 1, 0.1) # sigma

theta.star0 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24], 
                log(theta0[25:27]))
 
theta0.2 = c(rep(0.05, 6),
           9, 12, 0.3, # mu.gamma
           1.5, 6, 0.2, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 2, # betas
           0.01, 0.05, 0.05, # zero masses
           1, -1.8, 0, # mu
           1, 0.6, 0.1) # sigma

theta.star0.2 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24],
                log(theta0[25:27]))

### 4 state HMM? 4-th state is flapping flight -> directed movement but upwards

theta0.4 = c(rep(0.1, 12),
             9, 12, 0.3, 9, # mu.gamma
             1.5, 6, 0.2, 3, # sigma.gamma
             10, 1, 0.5, 1, # alphas
             55, 40, 2, 40, # betas
             0.01, 0.01, 0.01, 0.01, # zero masses
             1, -2, 0, 1, # mu
             0.5, 0.5, 0.1, 0.5) # sigma

theta.star0.4 = c(log(theta0.4[1:28]),
                  qlogis(theta0.4[29:32]),
                  theta0.4[33:36],
                  log(theta0.4[37:40]))


data2$height.fd[which(abs(data2$height.fd) > 20)] = NA
data2$angle[which(data2$angle == 1)] = runif(length(data2$angle[which(data2$angle == 1)]), 0.95, 1)
# Wir brauchen zero und one inflated beta Verteilung??

t1 = Sys.time()
mod = nlm(f = mllk, p = theta.star0, X = data2, N = 3, print.level = 2, iterlim = 1000)
Sys.time()-t1

t1 = Sys.time()
mod2 = nlm(f = mllk, p = theta.star0.2, X = data2, N = 3, print.level = 2, iterlim = 1000)
Sys.time()-t1

t1 = Sys.time()
mod4 = nlm(f = mllk, p = theta.star0.4, X = data2, N = 4, print.level = 2, iterlim = 1000)
Sys.time()-t1

-mllk(mod$estimate, X = data2, N = 3)
# -23738.48
-mllk(mod4$estimate, X = data2, N = 4)

theta.star = mod$estimate
# theta.star = mod2$estimate
theta.star = mod4$estimate

N = 3

Gamma = diag(N)
Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
Gamma = Gamma/rowSums(Gamma)
delta = solve(t(diag(N)-Gamma+1),rep(1,N)) # stationary

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



# Checking for global optimum ---------------------------------------------

# llk = rep(NA, 30)
# mods = vector("list")
# for (k in 1:30){
#   r_theta0 = theta0.3 = c(runif(6, 0, 0.5),
#                           c(9, 16) + runif(2,-5,5), runif(1,0,2),
#                           c(1.5, 4) + runif(2, -1.5, 1.5), runif(1,0,2),
#                           runif(1,5,10), runif(0,3), runif(0,2),
#                           c(55, 40) + runif(2, -20, 20), runif(1,0,4),
#                           runif(3,0,0.2),
#                           c(1, -0.9, 0) + runif(3, -2,2),
#                           runif(3,0,2))
#   r_theta.star0 = c(log(theta0[1:18]),
#                     qlogis(theta0[19:21]),
#                     theta0[22:24], 
#                     log(theta0[25:27]))
#   mods[[k]] = nlm(mllk, r_theta.star0, X = data2, N = 3, iterlim = 300)
#   llks[k] = -mods[[k]]$minimum
# }


# Plotting ----------------------------------------------------------------

color = c("deepskyblue", "orange", "forestgreen", "tomato3")

par(mfrow = c(1,1))
# total
hist(data2$step, prob = T, breaks = 300, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3],n = 500)

curve(
  delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]),
  add = T, lty = "dashed", lwd = 2, n = 500
)
# Mass around zero not captured accurately
# Maybe 4th state small movements??

hist(data2$angle, prob = T, breaks = 50, xlab = "Angle", xlim = c(0,1))
curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), add = T, lwd = 2, col = color[3])

curve(
  delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
    delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
    delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]),
  add = T, lty = "dashed", lwd = 2
)
# point masses on 0 and 1 for beta distribution!

hist(data2$height.fd, prob = T, breaks = 200, xlab = "Height.fd", xlim = c(-6,6))

curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dnorm(x, mu[3], sigma[3]), add = T, lwd = 2, col = color[3],n =1000)

curve(
  delta[1]*dnorm(x, mu[1], sigma[1])+
    delta[2]*dnorm(x, mu[2], sigma[2])+
    delta[3]*dnorm(x, mu[3], sigma[3]),
  add = T, lty = "dashed", lwd = 2, n=1000
)


# Plotting 4 state HMM ----------------------------------------------------

color = c("deepskyblue", "orange", "forestgreen", "tomato3")

par(mfrow = c(1,1))
# total
hist(data2$step, prob = T, breaks = 200, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3],n = 500)
curve(delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), add = T, lwd = 2, col = color[4],n = 500)

curve(
  delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
    delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]),
  add = T, lty = "dashed", lwd = 2, n = 500
)

hist(data2$angle, prob = T, breaks = 50, xlab = "Angle", xlim = c(0,1))
curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), add = T, lwd = 2, col = color[3])
curve(delta[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), add = T, lwd = 2, col = color[4], n = 500)

curve(
  delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
    delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
    delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3])+
    delta[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]),
  add = T, lty = "dashed", lwd = 2
)

hist(data2$height.fd, prob = T, breaks = 100, xlab = "Height.fd", xlim = c(-6,6))

curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dnorm(x, mu[3], sigma[3]), add = T, lwd = 2, col = color[3],n =1000)
curve(delta[4]*dnorm(x, mu[4], sigma[4]), add = T, lwd = 2, col = color[4],n =1000)

curve(
  delta[1]*dnorm(x, mu[1], sigma[1])+
    delta[2]*dnorm(x, mu[2], sigma[2])+
    delta[3]*dnorm(x, mu[3], sigma[3])+
    delta[4]*dnorm(x, mu[4], sigma[4]),
  add = T, lty = "dashed", lwd = 2, n=1000
)

states = viterbi(mod$estimate, data2, 3)

par(mfrow = c(3,1))
plot(data2$step[3001:6000], type = "h", col = color[states[3001:6000]])
plot(data2$angle[3001:6000], type = "h", col = color[states[3001:6000]])
plot(data2$height.fd[3001:6000], type = "h", col = color[states[3001:6000]])



# Looking at weird turning angles -----------------------------------------

pi_half_ind = which(data2$angle > 0.49 & data2$angle < 0.51)
pi_ind = which(data2$angle > 0.97 & data2$angle <= 1)

# Intervalle mit <= 2 Datenpunkten = NA setzen
for (value in pi_half_ind){
  print(d$angle[(value-30)+1:30])
  Sys.sleep(1)
}

for (value in pi_ind){
  print(d$angle[(value-30)+1:30])
  Sys.sleep(1)
}


# Getting better starting values for beta distribution --------------------

hist(data2$angle, prob = T, breaks = 50, xlab = "Angle", xlim = c(0,1))
curve(dbeta(x, shape1 = .5, shape2 = 2), add = T)

hist(data2$step, prob = T, breaks = 200, xlab = "Step length", xlim = c(0,30), ylim = c(0,1.4), main = "Resting")
curve(dgamma(x, shape = 0.5^2/0.4^2, scale = 0.4^2/0.5), add = T, n = 500)
