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