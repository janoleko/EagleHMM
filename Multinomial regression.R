# Multinomial regression --------------------------------------------------
library(fastDummies)
library(nnet)

colnames(data4)
data5 = data4 %>% select(states, l = landform.type)
dataf = dummy_cols(data5, select_columns = 'l')
colnames(dataf) = c("states", "landform", 
                    "cliff", "lower_slope",
                    "mountain_divide", "peak_ridge",
                    "upper_slope", "valley")

m = multinom(states ~ cliff + lower_slope + mountain_divide + peak_ridge +
               upper_slope + valley,
             data = dataf)
summary(m)