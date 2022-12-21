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

2*mllk_landform(mod10$estimate, X = data5, N = 4) + 2*length(mod10$estimate)
2*mllk_landform(mod10$estimate, X = data5, N = 4) + log(5000)*length(mod10$estimate)
