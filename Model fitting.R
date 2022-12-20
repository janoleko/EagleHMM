library(moveHMM)
library(tidyverse)
library(sn)

# Load in the data --------------------------------------------------------
setwd("/Users/jan-ole/R/HMM Project")
ider = read.csv("Ider_annotated.csv")
colnames(ider)

ider$height.fd = c(NA, diff(ider$height.above.msl))
# ider$turning2 = c(NA, diff(ider$heading*pi/180))

N = nrow(ider)/30

measurement = rep(1:N, each=30)

d = ider[,5:6]
d$ID = "ID"
colnames(d) = c("x", "y", "ID")
d = prepData(d)

index_stepNA <- (1:N)*30
d$step[index_stepNA] <- NA

index_angleNA <- (1:(N-1))*30+1
d$angle[c(index_stepNA, index_angleNA)] <- NA

# deleting wrong step lengths and angles
# for (i in 1:N){
#   d$step[30*i] = NA
#   d$angle[30*i] = NA
# }
# for (i in 1:(N-1)){
#   d$angle[30*i+1] = NA
# }

data = as.data.frame(cbind(ider,d, measurement))

# aggregation with mean
data2 = data %>% 
  dplyr::select(step, angle, height = height.above.msl, height.fd, measurement, x, y, day, landform, landform.type, elevation = elevation.m, temp = external.temperature, time = solar.time) %>% 
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
            day = mean(day, na.rm = T),
            landform = as.integer(names(which.max(table(landform)))),
            landform.type = as.character(names(which.max(table(na.omit(landform.type))))),
            elevation = mean(elevation, na.rm = T),
            temp = mean(temp, na.rm = T),
            time = mean(time, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(remove_step_height = (count_steps <= 3 | count_heights <= 3),
         remove_angle = count_angles <= 3)
# seasonal dummy?

data2$step[which(data2$remove_step_height == T)] = NA # set interval with barely any values to NA
data2$height.fd[which(data2$remove_step_height == T)] = NA # set interval with barely any values to NA
data2$angle[which(data2$remove_angle == T)] = NA # set interval with barely any values to NA

data2$height.fd[which(data2$height.fd > 20)] = NA
data2$height.fd[which(data2$height.fd < -110)] = NA

# get rid of point masses
set.seed(123)
data2$angle[which(data2$angle == 1)] = runif(length(data2$angle[which(data2$angle == 1)]), 0.98, 1)
data2$angle[which(data2$angle == 0)] = runif(length(data2$angle[which(data2$angle == 0)]), 0, 0.02)
data2$step[which(data2$step == 0)] = runif(length(data2$step[which(data2$step == 0)]), 0, 0.02)


# Data prep done ----------------------------------------------------------
















# More EDA ----------------------------------------------------------------

coloor = rep(1, 500)
coloor[which(data3$height.fd>0.5)] = 2
coloor[which(data3$height.fd < (-0.5))] = 3
coloor[which(abs(data3$height.fd) < 0.5 & data3$step > 5)] = 4

par(mfrow = c(3,1))
data3 = data2[4000:4500,]
plot(data3$step, type = "h", col = coloor)
plot(data3$angle, type = "h", col = coloor)
plot(data3$height.fd, type = "h", col = coloor)

par(mfrow = c(4,1))
hist(data3$step[which(coloor == 1)], prob = T, breaks = 20, xlim = c(0,20))
hist(data3$step[which(coloor == 2)], prob = T, breaks = 20, xlim = c(0,20))
hist(data3$step[which(coloor == 3)], prob = T, breaks = 20, xlim = c(0,20))
hist(data3$step[which(coloor == 4)], prob = T, breaks = 20, xlim = c(0,20))

hist(data3$angle[which(coloor == 1)], prob = T, breaks = 10, xlim = c(0,1))
hist(data3$angle[which(coloor == 2)], prob = T, breaks = 10, xlim = c(0,1))
hist(data3$angle[which(coloor == 3)], prob = T, breaks = 10, xlim = c(0,1))
hist(data3$angle[which(coloor == 4)], prob = T, breaks = 10, xlim = c(0,1))

hist(data3$height.fd[which(coloor == 1)], prob = T, breaks = 10, xlim = c(-5,5))
hist(data3$height.fd[which(coloor == 2)], prob = T, breaks = 10, xlim = c(-5,5))
hist(data3$height.fd[which(coloor == 3)], prob = T, breaks = 10, xlim = c(-5,5))
hist(data3$height.fd[which(coloor == 4)], prob = T, breaks = 10, xlim = c(-5,5))



# Fitting a model with turning angle --------------------------------------

theta0 = c(rep(0.1, 6),
           9, 13, 0.5, # mu.gamma
           1.5, 4, 0.4, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 2, # betas
           0.01, 0.01, 0.01, # zero masses
           1, -1, 0, # mu
           2, 1, 0.1) # sigma

theta.star0 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24], 
                log(theta0[25:27]))
 
theta0.2 = c(rep(0.1, 6),
             1, 8, 15, # mu.gamma
             1, 4, 10, # sigma.gamma
             0.5, 10, 1, # alphas
             2, 55, 40, # betas
             0, 0.2, -0.2, # xi
             0.05, 0.6, 0.6, # omega
             0, 5, -5) 

theta.star0.2 = c(log(theta0.2[1:18]),
                theta0.2[19:21], 
                log(theta0.2[22:24]),
                theta0.2[25:27])


theta0.2.4 = c(rep(0.1, 12),
              1, 8, 15, 10, # mu.gamma
              1, 4, 10, 8, # sigma.gamma
              0.5, 10, 1, 2, # alphas
              2, 55, 40, 40, # betas
              0, 0.2, -0.2, 0, # xi
              0.05, 0.5, 0.5, 0.6, # omega
              0, 5, -5, 0) 

theta.star0.2.4 = c(log(theta0.2.4[1:28]),
                  theta0.2.4[29:32], 
                  log(theta0.2.4[33:36]),
                  theta0.2.4[37:40])


theta0.2.4_ = c(rep(0.1, 12),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.7, 10, 1, 3, # alphas
               2, 55, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5) 

theta.star0.2.4_ = c(log(theta0.2.4_[1:28]),
                    theta0.2.4_[29:32], 
                    log(theta0.2.4_[33:36]),
                    theta0.2.4_[37:40])


# theta0.3 = c(rep(0.1, 6),
#              9, 12, 0.3, # mu.gamma
#              1.5, 6, 0.2, # sigma.gamma
#              10, 1, 0.5, # alphas
#              55, 40, 2, # betas
#              0.01, 0.01, 0.01, # zero masses
#              2, -2, 0, # mu
#              1, 1, 0.1) # sigma
# 
# theta.star0.3 = c(log(theta0.3[1:18]),
#                   qlogis(theta0.3[19:21]),
#                   log(theta0.3[22]),
#                   theta0.3[23:24],
#                   log(theta0.3[25:27]))

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


### HMM with covariate temperature

theta0.5 = c(-1.2, -3.27, 0.77, -1.75, -1.15, -0.72, # beta0
             rep(0, 6), # beta1
             9, 13, 0.5, # mu.gamma
             1.5, 4, 0.4, # sigma.gamma
             10, 1, 0.5, # alphas
             55, 40, 2, # betas
             1, -1.8, 0, # mu
             1, 0.6, 0.1,  # sigma
             3, 5 # delta
             )
theta.star0.5 = c(theta0.5[1:12],
                  log(theta0.5[13:24]),
                  theta0.5[25:27],
                  log(theta0.5[28:32]))


theta0.5 = c(Gamma[2:4,1], Gamma[c(1,3:4),2], Gamma[c(1:2,4),3], Gamma[1:3,4], rep(0,12),
                1, 8, 15, 2, # mu.gamma
                1, 4, 10, 2, # sigma.gamma
                0.7, 10, 1, 3, # alphas
                2, 55, 40, 40, # betas
                0, 0.2, -0.2, 0.4, # xi
                0.05, 0.5, 0.5, 0.5, # omega
                0, 5, -5, 5, # al
                1,1,1) # delta

theta.star0.5 = c(theta0.5[1:24],
                  log(theta0.5[25:40]),
                  theta0.5[41:44],
                  log(theta0.5[45:48]),
                  theta0.5[49:55])

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
mod = nlm(f = mllk, p = theta.star0, X = data2, N = 3, print.level = 2, iterlim = 1000)
Sys.time()-t1

mllk_na_sn(theta.star0.2, X = data2, N = 3)
t1 = Sys.time()
mod2 = nlm(f = mllk_na_sn, p = theta.star0.2, X = data2, N = 3, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
mod2 = readRDS("mod2.rds")
theta.star = mod2$estimate
states = viterbi_na_sn(theta.star, data2, N = 3)

# 4 State HMM with na handling and skew normal distribution (part of data set)
t1 = Sys.time()
mod2.4 = nlm(f = mllk_na_sn, p = theta.star0.2.4, X = data2[500:4000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod2.4$estimate
data4 = data2[500:4000,]
states = viterbi_na_sn(theta.star, data4, N = 4)


# 4 State HMM with na and skew norm (full data set) -----------------------

t1 = Sys.time()
mod2.4_ = nlm(f = mllk_na_sn, p = theta.star0.2.4, X = data2, N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod2.4_$estimate
states = viterbi_na_sn(theta.star, data2, N = 4)

t1 = Sys.time()
mod2.4_2 = nlm(f = mllk_na_sn, p = theta.star0.2.4, X = data2, N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod2.4_2$estimate
states = viterbi_na_sn(theta.star, data2, N = 4)

t1 = Sys.time()
mod2.4_3 = nlm(f = mllk_na_sn, p = theta.star0.2.4_, X = data2[1:5000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
mod2.4_3 = readRDS("mod2.4_3.rds")
theta.star = mod2.4_3$estimate
data4 = data2[1:5000,]
states = viterbi_na_sn(theta.star, data4, N = 4)
2*mllk_na_sn(mod2.4_3$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod2.4_3$estimate)


t1 = Sys.time()
mod5 = nlm(f = mllk_na_sn_cov, p = theta.star0.5, X = data2[1:5000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod5$estimate
states = viterbi_na_sn_cov(mod5$estimate, X = data2[1:5000,], N = 4)
2*mllk_na_sn_cov(mod5$estimate, X = data2[1:4000,], N = 4) + log(nrow(data2[1:4000,]))*length(mod5$estimate)


t1 = Sys.time()
mod6 = nlm(f = mllk_na_sn_cov2, p = theta.star0.5, X = data2[1:5000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8, hessian = T)
mod6.1 = nlm(f = mllk_na_sn_cov2, p = theta.star0.5, X = data2[1:5000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8, hessian = T)
Sys.time()-t1
theta.star = mod6$estimate
states = viterbi_na_sn_cov2(mod6$estimate, X = data2[1:5000,], N = 4)

t1 = Sys.time()
mod7 = nlm(f = mllk_na_sn_cov3, p = theta.star0.5, X = data2[1:5000,], N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod7$estimate
states = viterbi_na_sn_cov3(mod7$estimate, X = data2[1:5000,], N = 4)
# mod 7 ist quatsch


t1 = Sys.time()
mod8 = nlm(f = mllk_tod, p = theta.star0.tod, X = data4, N = 4, print.level = 2, iterlim = 1000, steptol = 1e-8)
Sys.time()-t1
theta.star = mod8$estimate
states = viterbi_tod(mod8$estimate, X = data2[1:5000,], N = 4)





# AIC:
2*mllk_na_sn(mod2.4_3$estimate, X = data2[1:5000,], N = 4) + 2*length(mod2.4_3$estimate)
2*mllk_na_sn_cov(mod5$estimate, X = data2[1:5000,], N = 4) + 2*length(mod5$estimate)
2*mllk_na_sn_cov(mod6$estimate, X = data2[1:5000,], N = 4) + 2*length(mod6$estimate)
2*mllk_tod(mod8$estimate, X = data2[1:5000,], N = 4) + 2*length(mod8$estimate)


# BIC:
2*mllk_na_sn(mod2.4_3$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod2.4_3$estimate)
2*mllk_na_sn_cov(mod5$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod5$estimate)
2*mllk_na_sn_cov(mod6$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod6$estimate)
2*mllk_tod(mod8$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod8$estimate)


# AIC and BIC prefer the covariate model
# AIC and BIC prefer the model with tod


theta.star = mod$estimate
theta.star = mod_test$estimate
# theta.star = mod2$estimate
theta.star = mod2.4$estimate

mllk_na_sn(mod2$estimate, X = data2, N = 3)
# 16859.65

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




# Checking for global optimum (parallel) ----------------------------------

llks = rep(NA, 100)
r_thetalist = list()
for (k in 1:100){
  r_theta0 = theta0.3 = c(runif(6, 0, 1),
                          c(9, 16) + runif(2,-8,8), runif(1,0,4),
                          c(1.5, 4) + runif(2, -3, 3), runif(1,0,4),
                          runif(1,4,12), runif(1,0,3), runif(1,0,2),
                          c(55, 40) + runif(2,-30, 30), runif(1,0,6),
                          runif(3,0,0.4),
                          c(1, -0.9, 0) + runif(3, -2.5, 2.5),
                          runif(3,0,3))
  r_theta.star0 = c(log(theta0[1:18]),
                    qlogis(theta0[19:21]),
                    theta0[22:24],
                    log(theta0[25:27]))
  r_thetalist[[k]] = r_theta.star0
}

library(parallel)
estimation = function(r_theta.star0){
  mod = tryCatch(nlm(mllk, r_theta.star0, X = data2, N = 3, iterlim = 250),
                 error = function(e) e)
  return(mod)
}

# mods = mclapply(X = r_thetalist, FUN = estimation, mc.cores = 8)

for (k in 1:100){
  llks[k] = -mods[[k]]$minimum
}



 # Plotting ----------------------------------------------------------------

color = c("deepskyblue", "orange", "forestgreen", "tomato3")

par(mfrow = c(1,1))
# total
hist(data2$step, prob = T, breaks = 100, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 1000)

curve(
  delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]),
  add = T, lty = "dashed", lwd = 2, n = 1000
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

states = viterbi_na_sn(theta.star, X = data2, N = 3)

par(mfrow = c(3,1))
plot(data2$step[3001:3200], type = "h", col = color[states[3001:3200]], ylab = "Step length")
plot(data2$angle[3001:3200], type = "h", col = color[states[3001:3200]], ylab = "Turning angle")
plot(data2$height.fd[3001:3200], type = "h", col = color[states][3001:3200], ylab = "Height(fd)")



# Plotting HMM with skew norm ---------------------------------------------
N = 3
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


par(mfrow = c(1,1))
# total
hist(data2$step, prob = T, breaks = 100, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3], n = 1000)

curve(
  delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]),
  add = T, lty = "dashed", lwd = 2, n = 1000
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

hist(data2$height.fd, prob = T, breaks = 200, xlab = "Height.fd", xlim = c(-10,10))

curve(delta[1]*dsn(x, xi[1], omega[1], alpha = al[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dsn(x, xi[2], omega[2], alpha = al[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dsn(x, xi[3], omega[3], alpha = al[3]), add = T, lwd = 2, col = color[3],n =1000)

curve(
  delta[1]*dnorm(x, mu[1], sigma[1])+
    delta[2]*dnorm(x, mu[2], sigma[2])+
    delta[3]*dnorm(x, mu[3], sigma[3]),
  add = T, lty = "dashed", lwd = 2, n=1000
)


# Plotting 4 State HMM with skew norm -------------------------------------
N = 4
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


color = c("blue", "orange", "forestgreen", "deepskyblue")

par(mfrow = c(1,1))
# total
hist(data4$step, prob = T, breaks = 200, xlab = "Step length", main = "Resting")
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

hist(data4$angle, prob = T, breaks = 50, xlab = "Angle", xlim = c(0,1))
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

hist(data4$height.fd, prob = T, breaks = 200, xlab = "Height.fd", xlim = c(-5,10))
curve(delta[1]*dsn(x, xi[1], omega[1], al[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dsn(x, xi[2], omega[2], al[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dsn(x, xi[3], omega[3], al[3]), add = T, lwd = 2, col = color[3])
curve(delta[4]*dsn(x, xi[4], omega[4], al[4]), add = T, lwd = 2, col = color[4])


curve(
  delta[1]*dsn(x, xi[1], omega[1], al[1])+
    delta[2]*dsn(x, xi[2], omega[2], al[2])+
    delta[3]*dsn(x, xi[3], omega[3], al[3])+
    delta[4]*dsn(x, xi[4], omega[4], al[4])
  ,
  add = T, lty = "dashed", lwd = 2, n=1000
)

color = c("deepskyblue", "orange", "springgreen4", "dodgerblue3")

data4 = data2[1:5000,]
par(mfrow = c(3,1))
plot(data4$step[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Step lenght")
legend(x = 4700, y = 30, legend=c("State 1", "State 2", "State 3", "State 4"),
        col=color, lty = 1, cex=1)
plot(data4$angle[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Turning angle")
plot(data4$height.fd[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Height (fd.)")


tempseq = seq(min(data4$temp), max(data4$temp), length.out = 500)

# Getting hypothetical stationary distribution and transition prob --------

delta = solve_gamma_na_sn_cov(theta.star, tempseq, N = 4)
transprobs = get_transprobs(theta.star, tempseq)

# Plotting hypothetical delta ----------------------------------------------

par(mfrow = c(1,1))
plot(tempseq, delta[,1], type = "l", lwd = 2, col = color[1], ylim = c(0,0.8), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "temperature")
lines(tempseq, delta[,2], lwd = 2, col = color[2])
lines(tempseq, delta[,3], lwd = 2, col = color[3])
lines(tempseq, delta[,4], lwd = 2, col = color[4])
legend(15, 0.8, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1, box.lwd = 0)


tempseq2 = seq(15, 50, length.out = 3)

N = 4
mu.g = exp(theta.star[2*(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[2*(N-1)*N+N+1:N]) # sds of gamma distributions

# beta distribution: Turning angle/ pi
alpha = exp(theta.star[2*(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[2*(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions

# normal distribution: Height first difference
xi = theta.star[2*(N-1)*N+4*N+1:N] # means of normal distributions
omega = exp(theta.star[2*(N-1)*N+5*N+1:N]) # sds of normal distributions
al = theta.star[2*(N-1)*N+6*N+1:N]


par(mfrow = c(3,3))
par(mar = c(4, 4, 1.5, 1.5))
for (i in 1:length(tempseq2)){
  delta = solve_gamma_na_sn_cov(theta.star, tempseq2[i], N = 4)
  curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), col = color[1], lwd = 1, xlim = c(0,30), ylim = c(0,0.35), ylab = "density", xlab = "step length", n = 300, main = paste("Temperature =", round(tempseq2[i], 0),"°"))
  curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), col = color[2], add = T, lwd = 1, n = 300)
  curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), col = color[3], add = T, lwd = 1, n = 300)
  curve(delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), col = color[4], add = T, lwd = 1, n = 300)
  curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
          delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+
          delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3])+
          delta[4]*dgamma(x, shape = mu.g[4]^2/sigma.g[4]^2, scale = sigma.g[4]^2/mu.g[4]), lty = "dashed", lwd = 2, add = T, n = 300)
  
  curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), col = color[1], lwd = 1, xlim = c(0,1), ylim = c(0,6.5), ylab = "density", xlab = "turning angle", n = 300, main = paste("Temperature =", round(tempseq2[i], 0),"°"))
  curve(delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), col = color[2], lwd = 1, add = T, n = 300)
  curve(delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), col = color[3], lwd = 1, add = T, n = 300)
  curve(delta[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), col = color[4], lwd = 1, add = T, n = 300)
  curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
          delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
          delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3])+
          delta[4]*dbeta(x, shape1 = alpha[4], shape2 = beta[4]), lty = "dashed", lwd = 2, add = T, n = 300)
  
  curve(delta[1]*dsn(x, xi = xi[1], omega = omega[1], alpha = al[1]), col = color[1], lwd = 1, xlim = c(-10,10), ylim = c(0,0.25), ylab = "density", xlab = "height.fd", n = 300, main = paste("Temperature =", round(tempseq2[i], 0),"°"))
  curve(delta[2]*dsn(x, xi = xi[2], omega = omega[2], alpha = al[2]), col = color[2], lwd = 1, add = T, n = 300)
  curve(delta[3]*dsn(x, xi = xi[3], omega = omega[3], alpha = al[3]), col = color[3], lwd = 1, add = T, n = 300)
  curve(delta[4]*dsn(x, xi = xi[4], omega = omega[4], alpha = al[4]), col = color[4], lwd = 1, add = T, n = 300)
  curve(delta[1]*dsn(x, xi = xi[1], omega = omega[1], alpha = al[1])+
          delta[2]*dsn(x, xi = xi[2], omega = omega[2], alpha = al[2])+
          delta[3]*dsn(x, xi = xi[3], omega = omega[3], alpha = al[3])+
          delta[4]*dsn(x, xi = xi[4], omega = omega[4], alpha = al[4]), lty = "dashed", lwd = 2, add = T, n = 300)
}


# Plotting transition probabilities ---------------------------------------

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(tempseq, transprobs[,i], main = NULL, 
       ylab = colnames(transprobs)[i], xlab = "temperature", 
       type = "l", lwd = 2, col = color[st])
}

# Coordinate plot with decoded states -------------------------------------

par(mfrow = c(1,1))
plot(data4$x, data4$y, col = color[states])

library(simply3d)
simply_scatter(data4$x, data4$y, data4$height, colorvar = color[states])


# Look at landforms at different locations --------------------------------

par(mfrow = c(4,1))
barplot(prop.table(table(data4[which(states == 1),]$landform)), main = "Resting")
barplot(prop.table(table(data4[which(states == 2),]$landform)), main = "Soaring")
barplot(prop.table(table(data4[which(states == 3),]$landform)), main = "Gliding")
barplot(prop.table(table(data4[which(states == 4),]$landform)), main = "Resting")

data4$states = states
par(mfrow = c(1,1))
boxplot(data4$elevation ~ states)
boxplot(data4$temp ~ states)
# soaring warmer, gliding colder
# generally more temperature variation when flying

View(head(ider,200))
