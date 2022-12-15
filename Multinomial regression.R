# Multinomial regression --------------------------------------------------
library(nnet)

data4$states = states
colnames(data4)
data5 = data4 %>% select(states, l_ = landform.type, temp, elevation)
data5$l_ = as.factor(data5$l_)
data5$states = relevel(as.factor(data5$states), ref = "2")
data5$elevation.fd = c(NA, diff(data5$elevation))

m = multinom(states ~ l_ + temp + elevation, data = data5)
summary(m)

