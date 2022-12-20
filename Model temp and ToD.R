# theta0.tt = c(Gamma[2:4,1], Gamma[c(1,3:4),2], Gamma[c(1:2,4),3], Gamma[1:3,4], rep(0,36),
#                1, 8, 15, 2, # mu.gamma
#                1, 4, 10, 2, # sigma.gamma
#                0.7, 10, 1, 3, # alphas
#                2, 55, 40, 40, # betas
#                0, 0.2, -0.2, 0.4, # xi
#                0.05, 0.5, 0.5, 0.5, # omega
#                0, 5, -5, 5, # al
#               0, 0, 0) # delta
# 
# theta.star0.tt = c(theta0.tt[1:48],
#                     log(theta0.tt[49:64]),
#                     theta0.tt[65:68],
#                     log(theta0.tt[69:72]),
#                     theta0.tt[73:79])

# t1 = Sys.time()
# mod9 = nlm(f = mllk_tt, p = theta.star0.tt, X = data4, N = 4, print.level = 2, iterlim = 2000, steptol = 1e-15)
# Sys.time()-t1

mod9 = readRDS("mod9.rds")
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



