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
states = viterbi_tod(mod8$estimate, X = data2[1:5000,], N = 4)2*mllk_tod(mod8$estimate, X = data2[1:5000,], N = 4) + 2*length(mod8$estimate)

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
plot(todseq, delta[,1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
lines(todseq, delta[,2], lwd = 2, col = color[2])
lines(todseq, delta[,3], lwd = 2, col = color[3])
lines(todseq, delta[,4], lwd = 2, col = color[4])
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
