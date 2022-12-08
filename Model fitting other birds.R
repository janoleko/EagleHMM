
# Fitting HMMs for the other birds ----------------------------------------

library(moveHMM)
library(tidyverse)

# Load in the data --------------------------------------------------------
setwd("/Users/jan-ole/R/HMM Project")
baidrag = read.csv("Baidrag_annotated.csv")
colnames(baidrag)

baidrag$height.fd = c(NA, diff(baidrag$height.above.msl))

N = nrow(baidrag)/30

measurement2 = rep(NA, nrow(ider))
for (i in 1:N){
  measurement2[(i-1)*30+1:30] = rep(i, 30)
}

d2 = baidrag[,5:6]
d2$ID = "ID"
colnames(d2) = c("x", "y", "ID")
d2 = prepData(d2)


# deleting wrong step lengths and angles
for (i in 1:N){
  d2$step[30*i] = NA
}
for (i in 1:(N-1)){
  d2$angle[30*i+1] = NA
}

data_ = as.data.frame(cbind(baidrag,d2, measurement2))

# aggregation with mean
data3 = data_ %>% 
  dplyr::select(step, angle, height = height.above.msl, height.fd, measurement2, x, y) %>% 
  group_by(measurement2) %>% 
  summarise(step = mean(step, na.rm = T)*1000, # converting to m/s
            count_angles = sum(!is.na(angle)), # find intervals where there are barely any values
            angle = abs(mean(angle, na.rm = T))/pi,
            height.fd = mean(diff(height), na.rm = T),
            x = mean(x, na.rm = T),
            y = mean(y, na.rm = T)) %>% 
  ungroup()
# seasonal dummy?

data3$angle[which(data3$count_angles <= 2)] = NA # set interval with barely any values to NA
data3$height.fd[which(data3$height.fd > 20)] = NA
data3$height.fd[which(data3$height.fd < -110)] = NA
data3$angle[which(data3$angle == 1)] = runif(length(data3$angle[which(data3$angle == 1)]), 0.95, 1)
data3$step[which(data3$step > 100)] = NA


# Histogramms (EDA) -------------------------------------------------------

hist(data3$step, prob = T, breaks = 1000)
hist(data2$step, prob = T, breaks = 100)
hist(data3$angle, prob = T, breaks = 100)
hist(data3$height.fd, prob = T, breaks = 100, xlim = c(-6,6))

par(mfrow = c(3,1))
plot(data3$step[1:500], type = "h")
plot(data3$angle[1:500], type = "h")
plot(data3$height.fd[1:500], type = "h")


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
