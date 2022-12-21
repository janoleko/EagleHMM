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


library(fastDummies)
help(package = "fastDummies")
datadummy = dummy_cols(data4, select_columns = "landform.type", remove_first_dummy = TRUE)
colnames(datadummy)[20:24] = c("lower_slope", "mountain_divide", "peak_ridge", "upper_slope", "valley")


