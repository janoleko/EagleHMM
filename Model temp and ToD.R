theta0.tt = c(Gamma[2:4,1], Gamma[c(1,3:4),2], Gamma[c(1:2,4),3], Gamma[1:3,4], rep(0,36),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.7, 10, 1, 3, # alphas
               2, 55, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5, # al
              0, 0, 0) # delta

theta.star0.tt = c(theta0.tt[1:48],
                    log(theta0.tt[49:64]),
                    theta0.tt[65:68],
                    log(theta0.tt[69:72]),
                    theta0.tt[73:79])

t1 = Sys.time()
mod9 = nlm(f = mllk_tt, p = theta.star0.tt, X = data4, N = 4, print.level = 2, iterlim = 2000, steptol = 1e-15)
Sys.time()-t1
theta.star = mod9$estimate
states = viterbi_tt(mod9$estimate, X = data2[1:5000,], N = 4)

# AIC
2*mllk_tt(mod9$estimate, X = data2[1:5000,], N = 4) + 2*length(mod9$estimate)
# BIC
2*mllk_tt(mod9$estimate, X = data2[1:5000,], N = 4) + log(nrow(data2[1:5000,]))*length(mod9$estimate)
# BIC better than model only with ToD

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
