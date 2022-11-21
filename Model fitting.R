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


# deleting wrong step lenghts and angles
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
            angle = abs(mean(angle, na.rm = T))/pi,
            height.fd = mean(diff(height), na.rm = T),
            x = mean(x, na.rm = T),
            y = mean(y, na.rm = T))

# Fitting a model with turning angle --------------------------------------

theta0 = c(rep(0.05, 6),
           9, 16, 0.3, # mu.gamma
           1.5, 6, 0.2, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 1.5, # betas
           0.01, 0.05, 0.05, # zero masses
           1, -0.5, 0, # mu
           1, 1, 0.1) # sigma

theta.star0 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24], 
                log(theta0[25:27]))

theta0.2 = c(rep(0.05, 6),
           10, 16, 0.3, # mu.gamma
           4, 6, 0.2, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 1.5, # betas
           0.01, 0.05, 0.05, # zero masses
           1, -0.5, 0, # mu
           1, 1, 0.1) # sigma

theta.star0.2 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24], 
                log(theta0[25:27]))


data2$height.fd[which(abs(data2$height.fd) > 20)] = NA
data2$angle[which(data2$angle == 1)] = .99 # Wir brauchen zero und one inflated beta Verteilung??

t1 = Sys.time()
mod = nlm(f = mllk, p = theta.star0, X = data2, N = 3, print.level = 2, iterlim = 1000)
Sys.time()-t1

t1 = Sys.time()
mod2 = nlm(f = mllk, p = theta.star0.2, X = data2, N = 3, print.level = 2, iterlim = 1000)
Sys.time()-t1

-mllk(mod$estimate, X = data2, N = 3)
-mllk(mod2$estimate, X = data2, N = 3) # same Maximum

theta.star = mod$estimate

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

color = c("cornflowerblue", "orange", "forestgreen")

par(mfrow = c(1,1))
# total
hist(data2$step, prob = T, breaks = 200, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3])

curve(
  delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]),
  add = T, lty = "dashed", lwd = 2
)

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


hist(data2$height.fd, prob = T, breaks = 100, xlab = "Height.fd")

curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dnorm(x, mu[3], sigma[3]), add = T, lwd = 2, col = color[3],n =1000)

curve(
  delta[1]*dnorm(x, mu[1], sigma[1])+
    delta[2]*dnorm(x, mu[2], sigma[2])+
    delta[3]*dnorm(x, mu[3], sigma[3]),
  add = T, lty = "dashed", lwd = 2, n=1000
)

states = viterbi(mod$estimate, data2, 3)

par(mfrow = c(3,1))
plot(data2$step[4700:5000], type = "h", col = color[states[4700:5000]])
plot(data2$angle[4700:5000], type = "h", col = color[states[4700:5000]])
plot(data2$height.fd[4700:5000], type = "h", col = color[states[4700:5000]])

