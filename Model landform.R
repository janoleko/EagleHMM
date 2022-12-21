theta0.l = c(rep(0, 72),
               1, 8, 15, 2, # mu.gamma
               1, 4, 10, 2, # sigma.gamma
               0.7, 10, 1, 3, # alphas
               2, 55, 40, 40, # betas
               0, 0.2, -0.2, 0.4, # xi
               0.05, 0.5, 0.5, 0.5, # omega
               0, 5, -5, 5, # al
              0, 0, 0) # delta

theta.star0.l = c(theta0.l[1:72],
                    log(theta0.l[73:88]),
                    theta0.l[89:92],
                    log(theta0.l[93:96]),
                    theta0.l[97:103])

data4 = data2[1:5000,]
data5 = dummy_cols(data4, select_columns = "landform.type", remove_first_dummy = TRUE)
colnames(data5)[19:23] = c("lower_slope", "mountain_divide", "peak_ridge", "upper_slope", "valley")

t1 = Sys.time()
mod10 = nlm(f = mllk_landform, p = theta.star0.l, X = data5, N = 4, print.level = 2, iterlim = 2000, steptol = 1e-15)
Sys.time()-t1
theta.star = mod10$estimate
states = viterbi_landform(mod10$estimate, data5, N = 4)

2*mllk_landform(mod10$estimate, X = data5, N = 4) + 2*length(mod10$estimate)
2*mllk_landform(mod10$estimate, X = data5, N = 4) + log(5000)*length(mod10$estimate)



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

stat = solve_gamma_landform(mod10$estimate, N = 4)

library(graphics)
par(mfrow = c(1,1))

plot(stat$landform, rep(2,6), pch = 16, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "landform type")
points(stat$landform, stat$delta[,1], type = "l", lty = "dashed", col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,1], pch = 16, col = color[1], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,2], type = "l", lty = "dashed", col = color[2], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,2], pch = 16, col = color[2], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,3], type = "l", lty = "dashed", col = color[3], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,3], pch = 16, col = color[3], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,4], type = "l", lty = "dashed", col = color[4], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")
points(stat$landform, stat$delta[,4], pch = 16, col = color[4], ylim = c(0,1), main = "Hypothetical stationary distribution", ylab = "Stationary state probabilities", xlab = "time of day")


# Transition probabilities

transprobs = get_transprobs_landform(mod10$estimate)

par(mfrow = c(4,4))
par(mar = c(4.5, 4.5, 1.5, 2))
for (i in 1:16){
  if(i %in% 1:4){st = 1}
  if(i %in% 5:8){st = 2}
  if(i %in% 9:12){st = 3}
  if(i %in% 13:16){st = 4}
  plot(transprobs$landform, rep(2,6), main = NULL, 
       ylab = colnames(transprobs$transprobs)[i], xlab = "landform type", 
       type = "l", lwd = 2, col = color[st], ylim = c(0,1))
  points(transprobs$landform, transprobs$transprobs[,i], type = "l", lty = "dashed", col = color[st],)
  points(transprobs$landform, transprobs$transprobs[,i], pch = 16, col = color[st],)
}
