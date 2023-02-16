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

cl = makeCluster(8); setDefaultCluster(cl=cl)
t1 = Sys.time()
mod_tod_bothways = optimParallel(par = theta.star0.tod, fn = mllk_tod, X = data7, N = 4, hessian = T, control = list(trace = 5, maxit = 10000))
Sys.time()-t1
stopCluster(cl)

theta.star = mod_tod_bothways$par