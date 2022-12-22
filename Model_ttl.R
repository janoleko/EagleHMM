theta0.ttl = c(rep(0, 108),
             1, 8, 15, 2, # mu.gamma
             1, 4, 10, 2, # sigma.gamma
             0.7, 10, 1, 3, # alphas
             2, 55, 40, 40, # betas
             0, 0.2, -0.2, 0.4, # xi
             0.05, 0.5, 0.5, 0.5, # omega
             0, 5, -5, 5, # al
             0, 0, 0) # delta

theta.star0.ttl = c(theta0.ttl[1:108],
                  log(theta0.ttl[109:124]),
                  theta0.ttl[125:128],
                  log(theta0.ttl[129:132]),
                  theta0.ttl[133:139])

data4 = data2[1:5000,]
data5 = dummy_cols(data4, select_columns = "landform.type", remove_first_dummy = TRUE)
colnames(data5)[19:23] = c("lower_slope", "mountain_divide", "peak_ridge", "upper_slope", "valley")

t1 = Sys.time()
mod11 = nlm(f = mllk_ttl, p = theta.star0.ttl, X = data5, N = 4, print.level = 2, iterlim = 2000, steptol = 1e-15)
Sys.time()-t1
theta.star = mod11$estimate
states = viterbi_ttl(mod11$estimate, X = data5, N = 4)

2*mllk_ttl(mod11$estimate, data5, N = 4) + 2*length(mod11$estimate)
2*mllk_ttl(mod11$estimate, data5, N = 4) + log(nrow(data5))*length(mod11$estimate)


# Results -----------------------------------------------------------------

color = c("deepskyblue", "orange", "springgreen4", "dodgerblue3")

# Scatterplot

par(mfrow = c(3,1))
plot(data5$step[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Step length")
legend(x = 4500, y = 30, legend=c("State 1", "State 2", "State 3", "State 4"),
       col=color, lty = 1, cex=1)
plot(data5$angle[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Turning angle")
plot(data5$height.fd[1:5000], pch = 20, col = color[states[1:5000]], ylab = "Height (fd.)")


# Hypothetical stationary for different lanforms

todseq = seq(min(data5$time), max(data5$time), length.out = 500)
tempseq = seq(min(data5$temp), max(data5$temp), length.out = 500)

h_delta1 = solve_gamma_ttl1(mod11$estimate, todseq, mean(data5$temp), 4)
h_delta2 = solve_gamma_ttl2(mod11$estimate, tempseq, mean(data5$time), 4)

# temperature is constant at its mean
par(mfrow = c(3,2))
landform = c("Cliff", "Lower slope", "Mountain/ Divide", "Peak/ Ridge", "Upper slope", "Valley")
for (i in 1:6){
  plot(todseq, h_delta1[,(i-1)*4+1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = paste("Landform:", landform[i]), ylab = "Stationary state probabilities", xlab = "time of day")
  lines(todseq, h_delta1[,(i-1)*4+2], lwd = 2, col = color[2])
  lines(todseq, h_delta1[,(i-1)*4+3], lwd = 2, col = color[3])
  lines(todseq, h_delta1[,(i-1)*4+4], lwd = 2, col = color[4])
  #legend(8.5, 1, legend=c("State 1", "State 2", "State 3", "State 4"),
        # col=color, lty = 1, cex=1, box.lwd = 0)
}

# time of day is constant at its mean
par(mfrow = c(3,2))
landform = c("Cliff", "Lower slope", "Mountain/ Divide", "Peak/ Ridge", "Upper slope", "Valley")
for (i in 1:6){
  plot(tempseq, h_delta2[,(i-1)*4+1], type = "l", lwd = 2, col = color[1], ylim = c(0,1), main = paste("Landform:", landform[i]), ylab = "Stationary state probabilities", xlab = "temperature")
  lines(tempseq, h_delta2[,(i-1)*4+2], lwd = 2, col = color[2])
  lines(tempseq, h_delta2[,(i-1)*4+3], lwd = 2, col = color[3])
  lines(tempseq, h_delta2[,(i-1)*4+4], lwd = 2, col = color[4])
  #legend(8.5, 1, legend=c("State 1", "State 2", "State 3", "State 4"),
  # col=color, lty = 1, cex=1, box.lwd = 0)
}


# Transition probabilities

transprobs_ttl1 = get_transprobs_ttl1(mod11$estimate, todseq, mean(data5$temp))

# Plot for time of day
par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
color2 = c("cadetblue", "cadetblue1", "chartreuse", "chartreuse4", "chocolate", "chocolate4")
for (i in 1:16){
  plot(todseq, transprobs_ttl1[((1:length(todseq))-1)*6+1,i], main = NULL, 
       ylab = colnames(transprobs_tt1)[i], xlab = "time of day", 
       type = "l", lwd = 2, col = color2[1], ylim = c(0,1))
  for (j in 2:6){
    lines(todseq, transprobs_ttl1[((1:length(todseq))-1)*6+j,i], lwd = 2, col = color2[j])
  }
}




