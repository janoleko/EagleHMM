
# Packages ----------------------------------------------------------------

library(tidyverse)
library(moveHMM)
library(CircStats)
library(simply3d)

# Load in the data --------------------------------------------------------

setwd("/Users/jan-ole/R/HMM Project")
ider = read.csv("Ider_annotated.csv")
colnames(ider)
ider$solar.time[1:1000]

ider$height.fd = c(NA, diff(ider$height.above.msl))

N = nrow(ider)/30

measurement = rep(1,30)
for (i in 2:N){
  measurement = c(measurement, rep(i, 30))
  if(i%%1000 == 0){
    print(i)
  }
}

d = ider[,5:6]
d$ID = "measurement"
colnames(d) = c("x", "y", "ID")
d = prepData(d)


# deleting wrong step lenghts and angles
for (i in 1:(nrow(d)/30)){
  d$step[30*i] = NA
}

for (i in 1:(nrow(d)/30-1)){
  d$angle[30*i+1] = NA
}



# Indizes -----------------------------------------------------------------

circle_ind = c(451:480, 481:510, 511:540, 541:570,  751:780 , 811:900, 961:990,
               1111:1140, 1951:1980 , 2191:2250, 2431:2460, 2521:2550, 5341:5370,
               5461:5520, 5551:5730, 5851:5940,6151:6180, 6871:6930, 7051:7140)
dir_ind = c(601:630, 781:810, 901:960 , 991:1110, 1141:1170, 2551:2610, 
            3211:3240, 3451:3540, 3601:3630, 3811:3900)
chill_ind = c(571:600 , 630:751 ,1171:1950, 1981:2190, 2251:2430, 2461:2520, 
              2611:3210, 3241:3450, 3541:3600, 3631:3810, 3901:5340, 5401:5460,
              5521:5550, 5731:5850, 5941:6150, 6181:6870, 6931:6870, 6931:7050,
              7141:7500)


# Indices Leon ------------------------------------------------------------

circle_ind = unique(c(circle_ind, 451:570, 601:630, 751:780, 811:900, 961:990, 1111:1140, 1951:1980, 
               2191:2250, 2431:2460,2521:2550,  5461:5520, 5851:5940, 6151:6180, 6871:6930, 7051:7140))

dir_ind = unique(c(dir_ind, 781:810, 901:960, 991:1020, 1081:1110, 1141:1170, 2551:2610, 
            3451:3540, 3811:3840, 5671:5730, 7201:7230))

chill_ind = unique(c(chill_ind, 571:600, 631:750, 1171:1950 , 1981:2190 , 2251:2430, 2461:2520, 
              2611:3210, 3241:3450, 3541:3600, 3631:3810, 3841:4260, 4291:5160,
              5221:5340, 5401:5460, 5521:5640, 5731:5850, 5941:6150,6181:6870,6961:6990, 7021:7050,
              7141:7200, 7231:7500))

fly_ind_l = c(1021:1080, 3601:3630, 4261:4290, 5161:5220, 5341:5370, 5641:5670, 6931:6960, 6991:7020)



# Indices Ole 2 -----------------------------------------------------------
# different part of TS

circle_ind2 = c(circle_ind, 75211:75270, 75361:75390, 75421:75480, 75511:75570, 75601:75720, 75781:75810,
                75871:75930, 75961:76020)
gliding_ind = c(dir_ind, 75121:75150, 75271:75360, 75391:75420, 75481:75510, 75571:75600, 75721:75780,
                75811:75840, 75931:75960, 76021:76050) 
rest_ind = c(chill_ind, 75001:75120, 76051:76080, 76111:76380)
fly_ind = c(fly_ind_l, 75181:75210, 75841:75870, 76081:76110)


label = rep(NA, nrow(ider))
label[circle_ind2] = 1
label[gliding_ind] = 2
label[rest_ind] = 3
label[fly_ind] = 2


data = as.data.frame(cbind(ider,d, measurement, label))

# aggregation with mean
data_aggr_mean = data %>% 
  dplyr::select(step, angle, height = height.above.msl, height.fd, measurement, label, x, y) %>% 
  group_by(measurement) %>% 
  summarise(step = mean(step, na.rm = T)*1000, # converting to m/s
            angle = abs(mean(angle, na.rm = T))/pi,
            height.fd = mean(diff(height), na.rm = T),
            x = mean(x, na.rm = T),
            y = mean(y, na.rm = T),
            label = mean(label))

color = c("cornflowerblue", "orange", "mediumspringgreen", "cyan1")

plot(data_aggr_mean$step[16:200], type = "h", col = color[data_aggr_mean$label[16:250]])
plot(data_aggr_mean$angle[16:200], type = "h", col = color[data_aggr_mean$label[16:250]])
plot(data_aggr_mean$height.fd[16:200], type = "h", col = color[data_aggr_mean$label[16:250]])
data_aggr_mean$height.fd[which(data_aggr_mean$height.fd > 20)] = NA


hist(data_aggr_mean$step, prob = T, breaks = 200, ylim = c(0,0.1))
hist(data_aggr_mean$angle, prob = T, breaks = 100)
hist(data_aggr_mean$height.fd, prob = T, breaks = 500, xlim = c(-3,3))




par(mfrow = c(1,1))
plot(data_aggr_mean$x, data_aggr_mean$y, lwd = 1)
ind = 2800
points(data_aggr_mean$x[ind], data_aggr_mean$y[ind], col = "blue", pch = 19)

# find out travelling indices
start = 2500
end = start + 2000
for (i in start:end){
  j = max(1, i-1000)
  plot(data_aggr_mean$x[j:i], data_aggr_mean$y[j:i], pch = 20, xlim = c(75,105), ylim = c(22,48))
  Sys.sleep(.1)
}
# 2500 is interesting 
2500*30


# aggregation with median (looks more promising)
data_aggr_median = data %>% 
  dplyr::select(step, angle, height = height.above.msl, height.fd, measurement, label) %>% 
  group_by(measurement) %>% 
  summarise(step = median(step, na.rm = T)*1000, # converting to m/s
            angle = abs(median(angle, na.rm = T))/pi,
            height.fd = median(height.fd, na.rm = T),
            label = mean(label))


# length of the time series in aggregated form:

7500/30 # 250

# we start at time point 16

# setting colors
color = c("cornflowerblue", "orange", "forestgreen", "cyan1")


plot(data_aggr_mean$var_angle[16:150], type = "h", col = color[data_aggr_median$label[16:150]])

plot(data_aggr_median$step[16:250], type = "h", col = color[data_aggr_median$label[16:250]], ylab = "median step length")

plot(data_aggr_mean$angle_m[16:250], type = "h", col = color[data_aggr_median$label[16:250]], ylab = "abs median angle / pi")
data_aggr_median$height.fd[which(abs(data_aggr_median$height.fd) > 30)] = NA # m/s > als 100 ist nicht möglich?
data_aggr_mean$height.fd[which(data_aggr_mean$height.fd > 20)] = NA # m/s > als 100 ist nicht möglich?

plot(data_aggr_mean$height.fd[16:250], type = "h", col = color[data_aggr_median$label[16:250]], ylab = "median height first diff.")
# abline(v = 93)

# Histograms for circles (rising)
hist(data_aggr_median$step[16:206][which(data_aggr_median$label[16:250] == 1)], prob = T, breaks = 10, xlab = "Step length circle")
curve(dgamma(x, shape = 12^2/5^2, scale = 5^2/12), add = T, col = color[1], lwd = 2)
sum(data_aggr_median$step[16:206][which(data_aggr_median$label[16:250] == 1)] == 0)/27
# point mass 0.03

hist(data_aggr_median$angle[16:206][which(data_aggr_median$label[16:250] == 1)], prob = T, breaks = 10, xlim = c(0,1))
curve(dbeta(x, shape1 = 2.5, shape2 = 20), add = T, col = color[1], lwd = 2)
sum(data_aggr_median$angle[16:206][which(data_aggr_median$label[16:250] == 1)] == 0)/27

hist(data_aggr_mean$height.fd[16:206][which(data_aggr_median$label[16:250] == 1)], prob = T, breaks = 30)
curve(dnorm(x, 1, 1), add = T, col = color[1], lwd = 2)


# Histograms for directed flight (down)
hist(data_aggr_median$step[16:206][which(data_aggr_median$label[16:250] == 2)], prob = T, breaks = 10)
curve(dgamma(x, shape = 25^2/10^2, scale = 10^2/25), add = T, col = color[2], lwd = 2)
sum(data_aggr_median$step[16:206][which(data_aggr_median$label[16:250] == 2)] == 0)/18
# 0.01

hist(data_aggr_median$angle[16:206][which(data_aggr_median$label[16:250] == 2)], prob = T, breaks = 20, xlim = c(0,1))
curve(dbeta(x, shape1 = .1, shape2 = 10), add = T, col = color[2], lwd = 2)
sum(data_aggr_median$angle[16:206][which(data_aggr_median$label[16:250] == 2)] == 0)/18

hist(data_aggr_median$height.fd[16:206][which(data_aggr_median$label[16:250] == 2)], prob = T, breaks = 20)
curve(dnorm(x, -2, 1), add = T, col = color[2], lwd = 2)


# Histograms for chilling 
hist(data_aggr_median$step[16:206][which(data_aggr_median$label[16:250] == 3)], prob = T, breaks = 10)
curve(dgamma(x, shape = 2^2/2^2, scale = 2^2/2), add = T, col = color[3], lwd = 2)
sum(na.omit(data_aggr_median$step[16:206][which(data_aggr_median$label[16:250] == 3)]) == 0)/160
# 0.65

hist(data_aggr_median$angle[16:206][which(data_aggr_median$label[16:250] == 3)], prob = T, breaks = 20, xlim = c(0,1))
curve(dbeta(x, shape1 = 1.2, shape2 = 1.2), add = T, col = color[3], lwd = 2)
sum(na.omit(data_aggr_median$angle[16:206][which(data_aggr_median$label[16:250] == 3)] == 0))/160

hist(data_aggr_median$height.fd[16:206][which(data_aggr_median$label[16:250] == 3)], prob = T, breaks = 20, xlim = c(-2,2))
curve(dnorm(x, 0, .05), add = T, col = color[3], lwd = 2)




# Looking at the intervals seperately -------------------------------------


par(mfrow = c(2,1))
plot(ider$location.long[451:480], ider$location.lat[451:480], pch = 20)
plot(ider$height.above.msl[451:480], pch = 20)
plot(ider$location.long[481:510], ider$location.lat[481:510], pch = 20)
plot(ider$location.long[511:540], ider$location.lat[511:540], pch = 20)
plot(ider$location.long[541:570], ider$location.lat[541:570], pch = 20)
plot(ider$location.long[571:600], ider$location.lat[571:600], pch = 20)
plot(ider$location.long[601:630], ider$location.lat[601:630], pch = 20) # 
# plot(ider$location.long[631:660], ider$location.lat[631:660], pch = 20) 
# plot(ider$location.long[661:690], ider$location.lat[661:690], pch = 20)
# plot(ider$location.long[691:720], ider$location.lat[691:720], pch = 20)
# plot(ider$location.long[721:750], ider$location.lat[721:750], pch = 20)
plot(ider$location.long[751:780], ider$location.lat[751:780], pch = 20)
plot(ider$location.long[781:810], ider$location.lat[781:810], pch = 20) #
plot(ider$location.long[811:840], ider$location.lat[811:840], pch = 20)
plot(ider$location.long[841:870], ider$location.lat[841:870], pch = 20)
plot(ider$location.long[871:900], ider$location.lat[871:900], pch = 20)
plot(ider$location.long[901:930], ider$location.lat[901:930], pch = 20) # 
plot(ider$location.long[931:960], ider$location.lat[931:960], pch = 20) # 
plot(ider$location.long[961:990], ider$location.lat[961:990], pch = 20)
plot(ider$location.long[991:1020], ider$location.lat[991:1020], pch = 20)
plot(ider$location.long[1021:1050], ider$location.lat[1021:1050], pch = 20) # 
# plot(ider$location.long[1051:1080], ider$location.lat[541:570], pch = 20) 
plot(ider$location.long[1081:1110], ider$location.lat[1081:1110], pch = 20) #
plot(ider$location.long[1111:1140], ider$location.lat[1111:1140], pch = 20)
plot(ider$location.long[1141:1170], ider$location.lat[1141:1170], pch = 20) #
# plot(ider$location.long[1171:1200], ider$location.lat[1171:1200], pch = 20)
# plot(ider$location.long[1201:1230], ider$location.lat[1201:1230], pch = 20)
# plot(ider$location.long[1231:1260], ider$location.lat[1231:1260], pch = 20)
# plot(ider$location.long[1261:1290], ider$location.lat[1261:1290], pch = 20)
# plot(ider$location.long[1291:1320], ider$location.lat[1291:1320], pch = 20)
# plot(ider$location.long[1321:1350], ider$location.lat[1321:1350], pch = 20)
# plot(ider$location.long[1351:1380], ider$location.lat[1351:1380], pch = 20)
# plot(ider$location.long[1381:1410], ider$location.lat[1381:1410], pch = 20)
# plot(ider$location.long[1411:1440], ider$location.lat[1411:1440], pch = 20)
# plot(ider$location.long[1441:1470], ider$location.lat[1441:1470], pch = 20)
# plot(ider$location.long[1471:1500], ider$location.lat[1471:1500], pch = 20)
# plot(ider$location.long[1501:1530], ider$location.lat[1501:1530], pch = 20)
# plot(ider$location.long[1531:1560], ider$location.lat[1531:1560], pch = 20)
# plot(ider$location.long[1561:1590], ider$location.lat[1561:1590], pch = 20)
# plot(ider$location.long[1591:1620], ider$location.lat[1591:1620], pch = 20)
# plot(ider$location.long[1621:1650], ider$location.lat[1621:1650], pch = 20)
# plot(ider$location.long[1651:1680], ider$location.lat[1651:1680], pch = 20)
# plot(ider$location.long[1681:1710], ider$location.lat[1681:1710], pch = 20)
# plot(ider$location.long[1711:1740], ider$location.lat[1711:1740], pch = 20)
# plot(ider$location.long[1741:1770], ider$location.lat[1741:1770], pch = 20)
# plot(ider$location.long[1771:1800], ider$location.lat[1771:1800], pch = 20)
# plot(ider$location.long[1801:1830], ider$location.lat[1801:1830], pch = 20)
# plot(ider$location.long[1831:1860], ider$location.lat[1831:1860], pch = 20)
# plot(ider$location.long[1861:1890], ider$location.lat[1861:1890], pch = 20)
# plot(ider$location.long[1891:1920], ider$location.lat[1891:1920], pch = 20)
# plot(ider$location.long[1921:1950], ider$location.lat[1921:1950], pch = 20)
plot(ider$location.long[1951:1980], ider$location.lat[1951:1980], pch = 20)
# plot(ider$location.long[1981:2010], ider$location.lat[1981:2010], pch = 20)
# plot(ider$location.long[2011:2040], ider$location.lat[2011:2040], pch = 20)
# plot(ider$location.long[2041:2070], ider$location.lat[2041:2070], pch = 20)
# plot(ider$location.long[2071:2100], ider$location.lat[2071:2100], pch = 20)
# plot(ider$location.long[2101:2130], ider$location.lat[2101:2130], pch = 20)
# plot(ider$location.long[2131:2160], ider$location.lat[2131:2160], pch = 20)
# plot(ider$location.long[2161:2190], ider$location.lat[2161:2190], pch = 20)
plot(ider$location.long[2191:2220], ider$location.lat[2191:2220], pch = 20)
plot(ider$location.long[2221:2250], ider$location.lat[2221:2250], pch = 20)
# plot(ider$location.long[2251:2280], ider$location.lat[2251:2280], pch = 20)
# plot(ider$location.long[2281:2310], ider$location.lat[2281:2310], pch = 20)
# plot(ider$location.long[2311:2340], ider$location.lat[2311:2340], pch = 20)
# plot(ider$location.long[2341:2370], ider$location.lat[2341:2370], pch = 20)
# plot(ider$location.long[2371:2400], ider$location.lat[2371:2400], pch = 20)
# plot(ider$location.long[2401:2430], ider$location.lat[2401:2430], pch = 20)
plot(ider$location.long[2431:2460], ider$location.lat[2431:2460], pch = 20)
# plot(ider$location.long[2461:2490], ider$location.lat[2461:2490], pch = 20)
# plot(ider$location.long[2491:2520], ider$location.lat[2491:2520], pch = 20)
plot(ider$location.long[2521:2550], ider$location.lat[2521:2550], pch = 20)
# plot(ider$location.long[2551:2580], ider$location.lat[2551:2580], pch = 20)
plot(ider$location.long[2581:2610], ider$location.lat[2581:2610], pch = 20) # 
# plot(ider$location.long[2611:2640], ider$location.lat[2611:2640], pch = 20)
# plot(ider$location.long[2641:2670], ider$location.lat[2641:2670], pch = 20)
# plot(ider$location.long[2671:2700], ider$location.lat[2671:2700], pch = 20)
# plot(ider$location.long[2701:2730], ider$location.lat[2701:2730], pch = 20)
# plot(ider$location.long[2731:2760], ider$location.lat[2731:2760], pch = 20)
plot(ider$location.long[2761:2790], ider$location.lat[2761:2790], pch = 20)
plot(ider$height.above.msl[2761:2790], pch = 20)

plot(ider$location.long[2791:2820], ider$location.lat[2791:2820], pch = 20)
plot(ider$height.above.msl[2791:2820], pch = 20)
# plot(ider$location.long[2821:2850], ider$location.lat[2821:2850], pch = 20)
# plot(ider$location.long[2851:2880], ider$location.lat[2851:2880], pch = 20)
# plot(ider$location.long[2881:2910], ider$location.lat[2881:2910], pch = 20)
# plot(ider$location.long[2911:2940], ider$location.lat[2911:2940], pch = 20)
# plot(ider$location.long[2941:2970], ider$location.lat[2941:2970], pch = 20)
# plot(ider$location.long[2971:3000], ider$location.lat[2971:3000], pch = 20)
# plot(ider$location.long[3001:3030], ider$location.lat[3001:3030], pch = 20)
# plot(ider$location.long[3031:3060], ider$location.lat[3031:3060], pch = 20)
# plot(ider$location.long[3061:3090], ider$location.lat[3061:3090], pch = 20)
# plot(ider$location.long[3091:3120], ider$location.lat[3091:3120], pch = 20)
# plot(ider$location.long[3121:3150], ider$location.lat[3121:3150], pch = 20)
# plot(ider$location.long[3151:3180], ider$location.lat[3151:3180], pch = 20)
# plot(ider$location.long[3181:3210], ider$location.lat[3181:3210], pch = 20)
plot(ider$location.long[3211:3240], ider$location.lat[3211:3240], pch = 20) # 
plot(ider$height.above.msl[3211:3240])
# plot(ider$location.long[3241:3270], ider$location.lat[3241:3270], pch = 20)
# plot(ider$location.long[3271:3300], ider$location.lat[3271:3300], pch = 20)
# plot(ider$location.long[3301:3330], ider$location.lat[3301:3330], pch = 20)
# plot(ider$location.long[3331:3360], ider$location.lat[3331:3360], pch = 20)
# plot(ider$location.long[3361:3390], ider$location.lat[3361:3390], pch = 20)
# plot(ider$location.long[3391:3420], ider$location.lat[3391:3420], pch = 20)
# plot(ider$location.long[3421:3450], ider$location.lat[3421:3450], pch = 20)
plot(ider$location.long[3451:3480], ider$location.lat[3451:3480], pch = 20) #
plot(ider$location.long[3481:3510], ider$location.lat[3481:3510], pch = 20) # 
plot(ider$location.long[3511:3540], ider$location.lat[3511:3540], pch = 20) # 
# plot(ider$location.long[3541:3570], ider$location.lat[3541:3570], pch = 20) 
# plot(ider$location.long[3541:3570], ider$location.lat[3541:3570], pch = 20)
# plot(ider$location.long[3571:3600], ider$location.lat[3571:3600], pch = 20)
plot(ider$location.long[3601:3630], ider$location.lat[3601:3630], pch = 20) # hochfliegen mit Flügelschlag? Siehe 3D Plot
# plot(ider$location.long[3631:3660], ider$location.lat[3631:3660], pch = 20)
# plot(ider$location.long[3661:3690], ider$location.lat[3661:3690], pch = 20)
# plot(ider$location.long[3691:3710], ider$location.lat[3691:3710], pch = 20)
# plot(ider$location.long[3711:3740], ider$location.lat[3711:3740], pch = 20)
# plot(ider$location.long[3741:3770], ider$location.lat[3741:3770], pch = 20)
# plot(ider$location.long[3771:3800], ider$location.lat[3771:3800], pch = 20)
plot(ider$location.long[3801:3830], ider$location.lat[3801:3830], pch = 20) # Fliegen
plot(ider$location.long[3831:3860], ider$location.lat[3831:3860], pch = 20) # Fliegen
plot(ider$location.long[3861:3890], ider$location.lat[3861:3890], pch = 20) # Fliegen





# EDA for different part of time series -----------------------------------

# start at 600.000
par(mfrow = c(2,1))

plot(ider$location.long[600001:600030], ider$location.lat[600001:600030], pch = 20)
plot(ider$height.above.msl[600001:600030], pch = 20) #3

plot(ider$location.long[600031:600060], ider$location.lat[600031:600060], pch = 20)
plot(ider$height.above.msl[600031:600060], pch = 20) #3

plot(ider$location.long[600061:600090], ider$location.lat[600061:600090], pch = 20)
plot(ider$height.above.msl[600061:600090], pch = 20) #3

plot(ider$location.long[600091:600120], ider$location.lat[600091:600120], pch = 20)
plot(ider$height.above.msl[600091:600120], pch = 20) #3

plot(ider$location.long[600121:600150], ider$location.lat[600121:600150], pch = 20)
plot(ider$height.above.msl[600121:600150], pch = 20) #3

plot(ider$location.long[600151:600180], ider$location.lat[600151:600180], pch = 20)
plot(ider$height.above.msl[600151:600180], pch = 20) #3

plot(ider$location.long[600181:600210], ider$location.lat[600181:600210], pch = 20)
plot(ider$height.above.msl[600181:600210], pch = 20) #3

plot(ider$location.long[600211:600240], ider$location.lat[600211:600240], pch = 20)
plot(ider$height.above.msl[600211:600240], pch = 20) #3

plot(ider$location.long[600241:600270], ider$location.lat[600241:600270], pch = 20)
plot(ider$height.above.msl[600241:600270], pch = 20) #3

plot(ider$location.long[600271:600300], ider$location.lat[600271:600300], pch = 20)
plot(ider$height.above.msl[600271:600300], pch = 20) #3

plot(ider$location.long[600301:600330], ider$location.lat[600301:600330], pch = 20)
plot(ider$height.above.msl[600301:600330], pch = 20) #3

plot(ider$location.long[600331:600360], ider$location.lat[600331:600360], pch = 20)
plot(ider$height.above.msl[600331:600360], pch = 20) #3

plot(ider$location.long[600361:600390], ider$location.lat[600361:600390], pch = 20)
plot(ider$height.above.msl[600361:600390], pch = 20) #only NAs

plot(ider$location.long[600391:600420], ider$location.lat[600391:600420], pch = 20)
plot(ider$height.above.msl[600391:600420], pch = 20) #3

plot(ider$location.long[600421:600450], ider$location.lat[600421:600450], pch = 20)
plot(ider$height.above.msl[600421:600450], pch = 20) #3

plot(ider$location.long[600451:600480], ider$location.lat[600451:600480], pch = 20)
plot(ider$height.above.msl[600451:600480], pch = 20) #3

plot(ider$location.long[600481:600510], ider$location.lat[600481:600510], pch = 20)
plot(ider$height.above.msl[600481:600510], pch = 20) #3

plot(ider$location.long[600511:600540], ider$location.lat[600511:600540], pch = 20)
plot(ider$height.above.msl[600511:600540], pch = 20) #3

plot(ider$location.long[600541:600570], ider$location.lat[600541:600570], pch = 20)
plot(ider$height.above.msl[600541:600570], pch = 20) #3

plot(ider$location.long[600571:600600], ider$location.lat[600571:600600], pch = 20)
plot(ider$height.above.msl[600571:600600], pch = 20) #3

plot(ider$location.long[600601:600630], ider$location.lat[600601:600630], pch = 20)
plot(ider$height.above.msl[600601:600630], pch = 20) #3

plot(ider$location.long[600631:600660], ider$location.lat[600631:600660], pch = 20)
plot(ider$height.above.msl[600631:600660], pch = 20) #3

plot(ider$location.long[600661:600690], ider$location.lat[600661:600690], pch = 20)
plot(ider$height.above.msl[600661:600690], pch = 20) #3

plot(ider$location.long[600691:600720], ider$location.lat[600691:600720], pch = 20)
plot(ider$height.above.msl[600691:600720], pch = 20) #3

plot(ider$location.long[600721:600750], ider$location.lat[600721:600750], pch = 20)
plot(ider$height.above.msl[600721:600750], pch = 20) #3

plot(ider$location.long[600751:600780], ider$location.lat[600751:600780], pch = 20)
plot(ider$height.above.msl[600751:600780], pch = 20) #3

plot(ider$location.long[600781:600810], ider$location.lat[600781:600810], pch = 20)
plot(ider$height.above.msl[600782:600810], pch = 20) #2 gliding

plot(ider$location.long[600811:600840], ider$location.lat[600811:600840], pch = 20)
plot(ider$height.above.msl[600811:600840], pch = 20) #3

plot(ider$location.long[600841:600870], ider$location.lat[600841:600870], pch = 20)
plot(ider$height.above.msl[600841:600870], pch = 20) #3

plot(ider$location.long[600871:600900], ider$location.lat[600871:600900], pch = 20)
plot(ider$height.above.msl[600871:600900], pch = 20) #3

plot(ider$location.long[600901:600930], ider$location.lat[600901:600930], pch = 20)
plot(ider$height.above.msl[600901:600930], pch = 20) #1 soaring

simply_scatter(ider$location.long[600901:600930], ider$location.lat[600901:600930], 
               ider$height.above.msl[600901:600930], colorvar = c(1, rep(0,29)))

plot(ider$location.long[600931:600960], ider$location.lat[600931:600960], pch = 20)
plot(ider$height.above.msl[600931:600960], pch = 20) # flapping?

plot(ider$location.long[600961:600990], ider$location.lat[600961:600990], pch = 20)
plot(ider$height.above.msl[600961:600990], pch = 20) #3

plot(ider$location.long[600991:601020], ider$location.lat[600991:601020], pch = 20)
plot(ider$height.above.msl[600991:601020], pch = 20) #3

plot(ider$location.long[601021:601050], ider$location.lat[601021:601050], pch = 20)
plot(ider$height.above.msl[601021:601050], pch = 20) #3

plot(ider$location.long[601051:601080], ider$location.lat[601051:601080], pch = 20)
plot(ider$height.above.msl[601051:601080], pch = 20) #?

plot(ider$location.long[601081:601110], ider$location.lat[601081:601110], pch = 20)
plot(ider$height.above.msl[601081:601110], pch = 20) #3

plot(ider$location.long[601111:601140], ider$location.lat[601111:601140], pch = 20)
plot(ider$height.above.msl[601111:601140], pch = 20) #3

plot(ider$location.long[601141:601170], ider$location.lat[601141:601170], pch = 20)
plot(ider$height.above.msl[601141:601170], pch = 20) #3

plot(ider$location.long[601171:601200], ider$location.lat[601171:601200], pch = 20)
plot(ider$height.above.msl[601171:601200], pch = 20) #3

plot(ider$location.long[601201:601230], ider$location.lat[601201:601230], pch = 20)
plot(ider$height.above.msl[601201:601230], pch = 20) # flapping flight again --> rising directed movement

plot(ider$location.long[601231:601260], ider$location.lat[601231:601260], pch = 20)
plot(ider$height.above.msl[601231:601260], pch = 20) #3

plot(ider$location.long[601261:601290], ider$location.lat[601261:601290], pch = 20)
plot(ider$height.above.msl[601261:601290], pch = 20) #3

plot(ider$location.long[601291:601320], ider$location.lat[601291:601320], pch = 20)
plot(ider$height.above.msl[601291:601320], pch = 20) #3

plot(ider$location.long[601321:601350], ider$location.lat[601321:601350], pch = 20)
plot(ider$height.above.msl[601321:601350], pch = 20) #3

plot(ider$location.long[601351:601380], ider$location.lat[601351:601380], pch = 20)
plot(ider$height.above.msl[601351:601380], pch = 20) #3

plot(ider$location.long[601381:601410], ider$location.lat[601381:601410], pch = 20)
plot(ider$height.above.msl[601381:601410], pch = 20) #3

plot(ider$location.long[601411:601440], ider$location.lat[601411:601440], pch = 20)
plot(ider$height.above.msl[601411:601440], pch = 20) #3

plot(ider$location.long[601441:601470], ider$location.lat[601441:601470], pch = 20)
plot(ider$height.above.msl[601441:601470], pch = 20) #3

plot(ider$location.long[601471:601500], ider$location.lat[601471:601500], pch = 20)
plot(ider$height.above.msl[601472:601500], pch = 20) #?

plot(ider$location.long[601501:601530], ider$location.lat[601501:601530], pch = 20)
plot(ider$height.above.msl[601501:601530], pch = 20) #3

plot(ider$location.long[601531:601560], ider$location.lat[601531:601560], pch = 20)
plot(ider$height.above.msl[601531:601560], pch = 20) #3

plot(ider$location.long[601561:601590], ider$location.lat[601561:601590], pch = 20)
plot(ider$height.above.msl[601561:601590], pch = 20) #3

plot(ider$location.long[601591:601620], ider$location.lat[601591:601620], pch = 20)
plot(ider$height.above.msl[601591:601620], pch = 20) # flapping flight

simply_scatter(ider$location.long[601591:601620], ider$location.lat[601591:601620], 
               ider$height.above.msl[601591:601620], colorvar = c(1, rep(0,29)))

plot(ider$location.long[601621:601650], ider$location.lat[601621:601650], pch = 20)
plot(ider$height.above.msl[601621:601650], pch = 20) #3

plot(ider$location.long[601651:601680], ider$location.lat[601651:601680], pch = 20)
plot(ider$height.above.msl[601651:601680], pch = 20) #3

plot(ider$location.long[601681:601710], ider$location.lat[601681:601710], pch = 20)
plot(ider$height.above.msl[601681:601710], pch = 20) # NAs

plot(ider$location.long[601711:601740], ider$location.lat[601711:601740], pch = 20)
plot(ider$height.above.msl[601711:601740], pch = 20) #3

plot(ider$location.long[601741:601770], ider$location.lat[601741:601770], pch = 20)
plot(ider$height.above.msl[601741:601770], pch = 20) #NAs

plot(ider$location.long[601771:601800], ider$location.lat[601771:601800], pch = 20)
plot(ider$height.above.msl[601771:601800], pch = 20) #3

plot(ider$location.long[601801:601830], ider$location.lat[601801:601830], pch = 20)
plot(ider$height.above.msl[601801:601830], pch = 20) #3

plot(ider$location.long[601831:601860], ider$location.lat[601831:601860], pch = 20)
plot(ider$height.above.msl[601831:601860], pch = 20) #3

plot(ider$location.long[601861:601890], ider$location.lat[601861:601890], pch = 20)
plot(ider$height.above.msl[601861:601890], pch = 20) #3

plot(ider$location.long[601891:601920], ider$location.lat[601891:601920], pch = 20)
plot(ider$height.above.msl[601891:601920], pch = 20) #?

plot(ider$location.long[601921:601950], ider$location.lat[601921:601950], pch = 20)
plot(ider$height.above.msl[601921:601950], pch = 20) #3

plot(ider$location.long[601951:601980], ider$location.lat[601951:601980], pch = 20)
plot(ider$height.above.msl[601951:601980], pch = 20) #3

plot(ider$location.long[601981:602010], ider$location.lat[601981:602010], pch = 20)
plot(ider$height.above.msl[601981:602010], pch = 20) #3

plot(ider$location.long[602011:602040], ider$location.lat[602011:602040], pch = 20)
plot(ider$height.above.msl[602011:602040], pch = 20) #3

plot(ider$location.long[602041:602070], ider$location.lat[602041:602070], pch = 20)
plot(ider$height.above.msl[602041:602070], pch = 20) #3

plot(ider$location.long[602071:602100], ider$location.lat[602071:602100], pch = 20)
plot(ider$height.above.msl[602071:602100], pch = 20) #3

plot(ider$location.long[602101:602130], ider$location.lat[602101:602130], pch = 20)
plot(ider$height.above.msl[602101:602130], pch = 20) #3

plot(ider$location.long[602131:602160], ider$location.lat[602131:602160], pch = 20)
plot(ider$height.above.msl[602131:602160], pch = 20) #3

plot(ider$location.long[602161:602190], ider$location.lat[602161:602190], pch = 20)
plot(ider$height.above.msl[602161:602190], pch = 20) #3

plot(ider$location.long[602191:602220], ider$location.lat[602191:602220], pch = 20)
plot(ider$height.fd[602191:602220], pch = 20) #?

plot(ider$location.long[602221:602250], ider$location.lat[602221:602250], pch = 20)
plot(ider$height.above.msl[602221:602250], pch = 20) #3

plot(ider$location.long[602251:602280], ider$location.lat[602251:602280], pch = 20)
plot(ider$height.above.msl[602251:602280], pch = 20) # Attack?

simply_scatter(ider$location.long[602251:602280], ider$location.lat[602251:602280], 
               ider$height.above.msl[602251:602280], colorvar = c(1, rep(2,29)))

plot(ider$location.long[602281:602310], ider$location.lat[602281:602310], pch = 20)
plot(ider$height.above.msl[602281:602310], pch = 20) #?

simply_scatter(ider$location.long[602281:602310], ider$location.lat[602281:602310], 
               ider$height.above.msl[602281:602310], colorvar = c(1, rep(2,29)))

plot(ider$location.long[602311:602340], ider$location.lat[602311:602340], pch = 20)
plot(ider$height.above.msl[602311:602340], pch = 20) #1

simply_scatter(ider$location.long[602311:602340], ider$location.lat[602311:602340], 
               ider$height.above.msl[602311:602340], colorvar = c(1, rep(2,29)))

plot(ider$location.long[602341:602370], ider$location.lat[602341:602370], pch = 20)
plot(ider$height.above.msl[602341:602370], pch = 20) #1

plot(ider$location.long[602371:602400], ider$location.lat[602371:602400], pch = 20)
plot(ider$height.above.msl[602371:602400], pch = 20) # random flight?

plot(ider$location.long[602401:602430], ider$location.lat[602401:602430], pch = 20)
plot(ider$height.above.msl[602401:602430], pch = 20) # 1

simply_scatter(ider$location.long[602401:602430], ider$location.lat[602401:602430], 
               ider$height.above.msl[602401:602430], colorvar = c(1, rep(2,29)))

plot(ider$location.long[602431:602460], ider$location.lat[602431:602460], pch = 20)
plot(ider$height.above.msl[602431:602460], pch = 20) # flying

plot(ider$location.long[602461:602490], ider$location.lat[602461:602490], pch = 20)
plot(ider$height.above.msl[602461:602490], pch = 20) #3

plot(ider$location.long[602491:602520], ider$location.lat[602491:602520], pch = 20)
plot(ider$height.above.msl[602491:602520], pch = 20) #3

plot(ider$location.long[602521:602550], ider$location.lat[602521:602550], pch = 20)
plot(ider$height.above.msl[602521:602550], pch = 20) #3

plot(ider$location.long[602551:602580], ider$location.lat[602551:602580], pch = 20)
plot(ider$height.above.msl[602551:602580], pch = 20) #3

plot(ider$location.long[602581:602610], ider$location.lat[602581:602610], pch = 20)
plot(ider$height.above.msl[602581:602610], pch = 20) #3

plot(ider$location.long[602611:602640], ider$location.lat[602611:602640], pch = 20)
plot(ider$height.above.msl[602611:602640], pch = 20) #3

plot(ider$location.long[602641:602670], ider$location.lat[602641:602670], pch = 20)
plot(ider$height.above.msl[602641:602670], pch = 20) #flying?

plot(ider$location.long[602671:602700], ider$location.lat[602671:602700], pch = 20)
plot(ider$height.above.msl[602671:602700], pch = 20) #3

plot(ider$location.long[602701:602730], ider$location.lat[602701:602730], pch = 20)
plot(ider$height.above.msl[602701:602730], pch = 20) #3

plot(ider$location.long[602731:602760], ider$location.lat[602731:602760], pch = 20)
plot(ider$height.above.msl[602731:602760], pch = 20) #3

plot(ider$location.long[602761:602790], ider$location.lat[602761:602790], pch = 20)
plot(ider$height.above.msl[602761:602790], pch = 20) #3

plot(ider$location.long[602791:602820], ider$location.lat[602791:602820], pch = 20)
plot(ider$height.above.msl[602791:602820], pch = 20) #3

plot(ider$location.long[602821:602850], ider$location.lat[602821:602850], pch = 20)
plot(ider$height.above.msl[602821:602850], pch = 20) #flying

plot(ider$location.long[602851:602880], ider$location.lat[602851:602880], pch = 20)
plot(ider$height.above.msl[602851:602880], pch = 20) #3

plot(ider$location.long[602881:602910], ider$location.lat[602881:602910], pch = 20)
plot(ider$height.above.msl[602881:602910], pch = 20) #3

plot(ider$location.long[602911:602940], ider$location.lat[602911:602940], pch = 20)
plot(ider$height.above.msl[602911:602940], pch = 20) #3

plot(ider$location.long[602941:602970], ider$location.lat[602941:602970], pch = 20)
plot(ider$height.above.msl[602941:602970], pch = 20) # flying

plot(ider$location.long[602971:603000], ider$location.lat[602971:603000], pch = 20)
plot(ider$height.above.msl[602971:603000], pch = 20) #3

plot(ider$location.long[603001:603030], ider$location.lat[603001:603030], pch = 20)
plot(ider$height.above.msl[603001:603030], pch = 20) #3

plot(ider$location.long[603031:603060], ider$location.lat[603031:603060], pch = 20)
plot(ider$height.above.msl[603031:603060], pch = 20) #3

plot(ider$location.long[603061:603090], ider$location.lat[603061:603090], pch = 20)
plot(ider$height.above.msl[603061:603090], pch = 20) #3

plot(ider$location.long[603091:603120], ider$location.lat[603091:603120], pch = 20)
plot(ider$height.above.msl[603091:603120], pch = 20) #3

plot(ider$location.long[603121:603150], ider$location.lat[603121:603150], pch = 20)
plot(ider$height.above.msl[603121:603150], pch = 20) # just NAs

plot(ider$location.long[603151:603180], ider$location.lat[603151:603180], pch = 20)
plot(ider$height.above.msl[603151:603180], pch = 20) #3

plot(ider$location.long[603181:603210], ider$location.lat[603181:603210], pch = 20)
plot(ider$height.above.msl[603181:603210], pch = 20) # flying

plot(ider$location.long[603211:603240], ider$location.lat[603211:603240], pch = 20)
plot(ider$height.above.msl[603211:603240], pch = 20) # flying

plot(ider$location.long[603241:603270], ider$location.lat[603241:603270], pch = 20)
plot(ider$height.above.msl[603241:603270], pch = 20) # flappy rising

simply_scatter(ider$location.long[603241:603270], ider$location.lat[603241:603270], 
               ider$height.above.msl[603241:603270], colorvar = c(1, rep(2,29)))

plot(ider$location.long[603271:603300], ider$location.lat[603271:603300], pch = 20)
plot(ider$height.above.msl[603271:603300], pch = 20) # 3

plot(ider$location.long[603301:603330], ider$location.lat[603301:603330], pch = 20)
plot(ider$height.above.msl[603301:603330], pch = 20) # 3

plot(ider$location.long[603331:603360], ider$location.lat[603331:603360], pch = 20)
plot(ider$height.above.msl[603331:603360], pch = 20) # 3

plot(ider$location.long[603361:603390], ider$location.lat[603361:603390], pch = 20)
plot(ider$height.above.msl[603361:603390], pch = 20) # 3

plot(ider$location.long[603391:603420], ider$location.lat[603391:603420], pch = 20)
plot(ider$height.above.msl[603391:603420], pch = 20) # 3

plot(ider$location.long[603421:603450], ider$location.lat[603421:603450], pch = 20)
plot(ider$height.above.msl[603421:603450], pch = 20) # 3

plot(ider$location.long[603451:603480], ider$location.lat[603451:603480], pch = 20)
plot(ider$height.above.msl[603451:603480], pch = 20) # 3

plot(ider$location.long[603481:603510], ider$location.lat[603481:603510], pch = 20)
plot(ider$height.above.msl[603481:603510], pch = 20) # 3

plot(ider$location.long[603511:603540], ider$location.lat[603511:603540], pch = 20)
plot(ider$height.above.msl[603511:603540], pch = 20) # 3

plot(ider$location.long[603541:603570], ider$location.lat[603541:603570], pch = 20)
plot(ider$height.above.msl[603541:603570], pch = 20) # 3


# New part of Ts ----------------------------------------------------------

75030/30
par(mfrow = c(2,1))

plot(ider$location.long[75001:75030], ider$location.lat[75001:75030], pch = 20)
plot(ider$height.above.msl[75001:75030], pch = 20) # 3
data$step[75001:75030]*1000

plot(ider$location.long[75031:75060], ider$location.lat[75031:75060], pch = 20)
plot(ider$height.above.msl[75031:75060], pch = 20) # 3

plot(ider$location.long[75061:75090], ider$location.lat[75061:75090], pch = 20)
plot(ider$height.above.msl[75061:75090], pch = 20) # 3

plot(ider$location.long[75091:75120], ider$location.lat[75091:75120], pch = 20)
plot(ider$height.above.msl[75091:75120], pch = 20) # 3
data$step[75091:75120]*1000

plot(ider$location.long[75121:75150], ider$location.lat[75121:75150], pch = 20)
plot(ider$height.above.msl[75121:75150], pch = 20) #2 gliding

plot(ider$location.long[75151:75180], ider$location.lat[75151:75180], pch = 20)
plot(ider$height.above.msl[75151:75180], pch = 20) # what is this?
data$step[75151:75180]*1000
simply_scatter(ider$location.long[75151:75180], ider$location.lat[75151:75180], 
               ider$height.above.msl[75151:75180], colorvar = c(1, rep(2,29)))

plot(ider$location.long[75181:75210], ider$location.lat[75181:75210], pch = 20)
plot(ider$height.above.msl[75181:75210], pch = 20) # flying?

plot(ider$location.long[75211:75240], ider$location.lat[75211:75240], pch = 20)
plot(ider$height.above.msl[75211:75240], pch = 20) # 1 soaring

plot(ider$location.long[75241:75270], ider$location.lat[75241:75270], pch = 20)
plot(ider$height.above.msl[75241:75270], pch = 20) # 1 soaring
simply_scatter(ider$location.long[75241:75270], ider$location.lat[75241:75270], 
               ider$height.above.msl[75241:75270], colorvar = c(1, rep(2,29)))

plot(ider$location.long[75271:75300], ider$location.lat[75271:75300], pch = 20)
plot(ider$height.above.msl[75271:75300], pch = 20) # 2 gliding (am ende anfang soaring)

plot(ider$location.long[75301:75330], ider$location.lat[75301:75330], pch = 20)
plot(ider$height.above.msl[75301:75330], pch = 20) # 2 gliding des todes

plot(ider$location.long[75331:75360], ider$location.lat[75331:75360], pch = 20)
plot(ider$height.above.msl[75331:75360], pch = 20) # 2 gliding

plot(ider$location.long[75361:75390], ider$location.lat[75361:75390], pch = 20)
plot(ider$height.above.msl[75361:75390], pch = 20) # 1 soaring

plot(ider$location.long[75391:75420], ider$location.lat[75391:75420], pch = 20)
plot(ider$height.above.msl[75391:75420], pch = 20) # 2 gliding
simply_scatter(ider$location.long[75391:75420], ider$location.lat[75391:75420], 
               ider$height.above.msl[75391:75420], colorvar = c(1, rep(2,29)))

plot(ider$location.long[75421:75450], ider$location.lat[75421:75450], pch = 20)
plot(ider$height.above.msl[75421:75450], pch = 20) # 1 soaring
simply_scatter(ider$location.long[75421:75450], ider$location.lat[75421:75450], 
               ider$height.above.msl[75421:75450], colorvar = c(1, rep(2,29)))

plot(ider$location.long[75451:75480], ider$location.lat[75451:75480], pch = 20)
plot(ider$height.above.msl[75451:75480], pch = 20) # 1 soaring

plot(ider$location.long[75481:75510], ider$location.lat[75481:75510], pch = 20)
plot(ider$height.above.msl[75481:75510], pch = 20) # 2 gliding

plot(ider$location.long[75511:75540], ider$location.lat[75511:75540], pch = 20)
plot(ider$height.above.msl[75511:75540], pch = 20) # 1 

plot(ider$location.long[75541:75570], ider$location.lat[75541:75570], pch = 20)
plot(ider$height.above.msl[75541:75570], pch = 20) # 1 
simply_scatter(ider$location.long[75541:75570], ider$location.lat[75541:75570], 
               ider$height.above.msl[75541:75570], colorvar = c(1, rep(2,29)))

plot(ider$location.long[75571:75600], ider$location.lat[75571:75600], pch = 20)
plot(ider$height.above.msl[75571:75600], pch = 20) # 2 gliding

plot(ider$location.long[75601:75630], ider$location.lat[75601:75630], pch = 20)
plot(ider$height.above.msl[75601:75630], pch = 20) # 1

plot(ider$location.long[75631:75660], ider$location.lat[75631:75660], pch = 20)
plot(ider$height.above.msl[75631:75660], pch = 20) # 1

plot(ider$location.long[75661:75690], ider$location.lat[75661:75690], pch = 20)
plot(ider$height.above.msl[75661:75690], pch = 20) # 1

plot(ider$location.long[75691:75720], ider$location.lat[75691:75720], pch = 20)
plot(ider$height.above.msl[75691:75720], pch = 20) # 1

plot(ider$location.long[75721:75750], ider$location.lat[75721:75750], pch = 20)
plot(ider$height.above.msl[75721:75750], pch = 20) # 2 gliding

plot(ider$location.long[75751:75780], ider$location.lat[75751:75780], pch = 20)
plot(ider$height.above.msl[75751:75780], pch = 20) # 2 gliding

plot(ider$location.long[75781:75810], ider$location.lat[75781:75810], pch = 20)
plot(ider$height.above.msl[75781:75810], pch = 20) # 1

plot(ider$location.long[75811:75840], ider$location.lat[75811:75840], pch = 20)
plot(ider$height.above.msl[75811:75840], pch = 20) # 2

plot(ider$location.long[75841:75870], ider$location.lat[75841:75870], pch = 20)
plot(ider$height.above.msl[75841:75870], pch = 20) # flying

plot(ider$location.long[75871:75900], ider$location.lat[75871:75900], pch = 20)
plot(ider$height.above.msl[75871:75900], pch = 20) #1

plot(ider$location.long[75901:75930], ider$location.lat[75901:75930], pch = 20)
plot(ider$height.above.msl[75901:75930], pch = 20) #1

plot(ider$location.long[75931:75960], ider$location.lat[75931:75960], pch = 20)
plot(ider$height.above.msl[75931:75960], pch = 20) # 2 gliding

plot(ider$location.long[75961:75990], ider$location.lat[75961:75990], pch = 20)
plot(ider$height.above.msl[75961:75990], pch = 20) # 1

plot(ider$location.long[75991:76020], ider$location.lat[75991:76020], pch = 20)
plot(ider$height.above.msl[75991:76020], pch = 20) # 1

plot(ider$location.long[76021:76050], ider$location.lat[76021:76050], pch = 20)
plot(ider$height.above.msl[76021:76050], pch = 20) # 2 oder flying

plot(ider$location.long[76051:76080], ider$location.lat[76051:76080], pch = 20)
plot(ider$height.above.msl[76051:76080], pch = 20) # 3
data$step[76051:76080]*1000

plot(ider$location.long[76081:76110], ider$location.lat[76081:76110], pch = 20)
plot(ider$height.above.msl[76081:76110], pch = 20) # flying

plot(ider$location.long[76111:76140], ider$location.lat[76111:76140], pch = 20)
plot(ider$height.above.msl[76111:76140], pch = 20) # 3
data$step[76111:76140]*1000

plot(ider$location.long[76141:76170], ider$location.lat[76141:76170], pch = 20)
plot(ider$height.above.msl[76141:76170], pch = 20) # 3

plot(ider$location.long[76171:76200], ider$location.lat[76171:76200], pch = 20)
plot(ider$height.above.msl[76171:76200], pch = 20) # 3
data$step[76171:76200]*1000

plot(ider$location.long[76201:76230], ider$location.lat[76201:76230], pch = 20)
plot(ider$height.above.msl[76201:76230], pch = 20) # 3

plot(ider$location.long[76231:76260], ider$location.lat[76231:76260], pch = 20)
plot(ider$height.above.msl[76231:76260], pch = 20) # 3

plot(ider$location.long[76261:76290], ider$location.lat[76261:76290], pch = 20)
plot(ider$height.above.msl[76261:76290], pch = 20) # 3

plot(ider$location.long[76291:76320], ider$location.lat[76291:76320], pch = 20)
plot(ider$height.above.msl[76291:76320], pch = 20) # 3

plot(ider$location.long[76321:76350], ider$location.lat[76321:76350], pch = 20)
plot(ider$height.above.msl[76321:76350], pch = 20) # 3

plot(ider$location.long[76351:76380], ider$location.lat[76351:76380], pch = 20)
plot(ider$height.above.msl[76351:76380], pch = 20) # 3

76380/30


circle_ind2 = c(circle_ind, 75211:75270, 75361:75390, 75421:75480, 75511:75570, 75601:75720, 75781:75810,
               75871:75930, 75961:76020)
gliding_ind = c(dir_ind, 75121:75150, 75271:75360, 75391:75420, 75481:75510, 75571:75600, 75721:75780,
                 75811:75840, 75931:75960, 76021:76050) 
rest_ind = c(chill_ind, 75001:75120, 76051:76080, 76111:76380)
fly_ind = c(75181:75210, 75841:75870, 76081:76110)


label = rep(NA, nrow(ider))
label[circle_ind2] = 1
label[gliding_ind] = 2
label[rest_ind] = 3
label[fly_ind] = 2

par(mfrow = c(3,1))
plot(data_aggr_mean$step[2501:2546], type = "h", lwd = 2, col = color[data_aggr_mean$label[2501:2546]])
plot(data_aggr_mean$angle[2501:2546], type = "h", lwd = 2, col = color[data_aggr_mean$label[2501:2546]])
plot(data_aggr_mean$height.fd[2501:2546], type = "h", lwd = 2, col = color[data_aggr_mean$label[2501:2546]])

# State 1: Soaring: Medium step length, angle close to 0.2 (small variance), height.fd > 0
# State 2: Gliding (down): Long step length, small angles, height.fd < 0
# State 3: Resting: very small step length, large angles with high variance, height.fd closely around 0 
# State 4: Flying (~same height): Large step length, very small angles, height.fd around zero

par(mfrow = c(4,1))

# Steplength
hist(data_aggr_mean$step[c(16:200,2501:2546)][which(data_aggr_mean$label[c(16:200,2501:2546)] == 1)], prob = T, breaks = 20, xlab = "Step length", xlim = c(0,30), main = "Soaring")
curve(dgamma(x, shape = 9^2/1.5^2, scale = 1.5^2/9), add = T)

hist(data_aggr_mean$step[c(16:200,2501:2546)][which(data_aggr_mean$label[c(16:200,2501:2546)] == 2)], prob = T, breaks = 20, xlab = "Step length", xlim = c(0,30), main = "Gliding")
curve(dgamma(x, shape = 16^2/6^2, scale = 6^2/16), add = T)

hist(data_aggr_mean$step[c(16:200,2501:2546)][which(data_aggr_mean$label[c(16:200,2501:2546)] == 3)], prob = T, breaks = 10, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(dgamma(x, shape = 0.3^2/0.2^2, scale = 0.2^2/0.3), add = T)

# hist(data_aggr_mean$step[c(2501:2546)][which(data_aggr_mean$label[c(2501:2546)] == 4)], prob = T, breaks = 30, xlab = "Step length", main = "Flying", xlim = c(0,20))
# curve(dgamma(x, shape = 12^2/5^2, scale = 5^2/12), add = T)


# Angle
hist(data_aggr_mean$angle[c(16:200, 2501:2546)][which(data_aggr_mean$label[c(16:200, 2501:2546)] == 1)], prob = T, breaks = 15, xlab = "Angle", xlim = c(0,1), main = "Soaring")
curve(dbeta(x, shape1 = 10, shape2 = 55), add = T)
data_aggr_mean$angle[c(2501:2546)][which(data_aggr_mean$label[c(2501:2546)] == 1)] == 0

hist(data_aggr_mean$angle[c(16:200, 2501:2546)][which(data_aggr_mean$label[c(16:200, 2501:2546)] == 2)], prob = T, breaks = 30, xlab = "Angle", xlim = c(0,1), main = "Gliding")
curve(dbeta(x, shape1 = 1, shape2 = 40), add = T)
data_aggr_mean$angle[c(2501:2546)][which(data_aggr_mean$label[c(2501:2546)] == 2)] == 0

hist(data_aggr_mean$angle[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 3)], prob = T, breaks = 60, xlab = "Angle", xlim = c(0,1), main = "Resting")
curve(dbeta(x, shape1 = .5, shape2 = 1.5), add = T)
data_aggr_mean$angle[c(2501:2546)][which(data_aggr_mean$label[c(2501:2546)] == 3)] == 0

# hist(data_aggr_mean$angle[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 4)], prob = T, breaks = 30, xlab = "Angle", xlim = c(0,1), main = "Flying")
# curve(dbeta(x, shape1 = 1, shape2 = 20), add = T)
# data_aggr_mean$angle[c(2501:2546)][which(data_aggr_mean$label[c(2501:2546)] == 4)] == 0



# Height fd
hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 1)], prob = T, breaks = 20, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Soaring")
curve(dnorm(x, 1, 0.7), add = T)

hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 2)], prob = T, breaks = 200, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Gliding")
curve(dnorm(x, -0.7, 0.7), add = T)

hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 3)], prob = T, breaks = 50, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Resting")
curve(dnorm(x, 0, 0.1), add = T)

# hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 4)], prob = T, breaks = 10, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Flying")
# curve(dnorm(x, 0, .3), add = T)

data_aggr_mean$step[which(data_aggr_mean$step == 0)] = 0.01


# Fitting a model without turning angle -----------------------------------


theta0 = c(rep(0.05, 12),
           8, 15, 0.3, 12, # mu.gamma
           0.7, 5, 0.2, 6, # sigma.gamma
           1, -1.3, 0, 0, # mu
           0.3, 0.5, 0.2, 0.3) 

theta.star0 = c(log(theta0[1:20]),
                theta0[21:24],
                log(theta0[25:28]))

data2 = data_aggr_mean
data2$height.fd[which(abs(data2$height.fd) > 20)] = NA

mllk2(theta.star0, data2, N = 4)

mod = nlm(f = mllk2, p = theta.star0, X = data2, N = 4, print.level = 2, iterlim = 10000)

# mod = optim(par = theta.star0, fn = mllk2, X = data2, N = 4, control = list(trace = 2, maxit = 20000))
# mod2 = optim(par = theta.star0, fn = mllk2, method = "BFGS", X = data2, N = 4, control = list(trace = 2, maxit = 20000))

theta.star = mod$estimate
N = 4
Gamma = diag(N)
Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
Gamma = Gamma/rowSums(Gamma)
delta = solve(t(diag(N)-Gamma+1),rep(1,N), tol = 1e-20) # stationary

mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions

# normal distribution: Height first difference
mu = theta.star[(N-1)*N+2*N+1:N] # means of normal distributions
sigma = exp(theta.star[(N-1)*N+3*N+1:N]) # sds of normal distributions




# Fitting a model with turning angle --------------------------------------

theta0 = c(rep(0.05, 6),
           9, 16, 0.3, # mu.gamma
           1.5, 6, 0.2, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 1.5, # betas
           0.01, 0.05, 0.05, # zero masses
           1, -0.5, 0, # mu
           1, 1, 0.1) # sigma

# Wir brauchen zero und one inflated beta Verteilung??

theta.star0 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24], 
                log(theta0[25:27]))

theta0.2 = c(rep(0.05, 6),
           9, 16, 0.3, # mu.gamma
           1.5, 6, 0.2, # sigma.gamma
           10, 1, 0.5, # alphas
           55, 40, 1.5, # betas
           0.01, 0.05, 0.05, # zero masses
           1, -0.9, 0, # mu
           0.5, 0.5, 0.1) # sigma

# Wir brauchen zero und one inflated beta Verteilung??

theta.star0.2 = c(log(theta0[1:18]),
                qlogis(theta0[19:21]),
                theta0[22:24], 
                log(theta0[25:27]))

theta0.3 = c(rep(0.05, 6),
             9, 16, 0.3, # mu.gamma
             1.5, 4, 0.2, # sigma.gamma
             10, 1, 0.5, # alphas
             55, 40, 1.5, # betas
             0.01, 0.05, 0.05, # zero masses
             1, -0.9, 0, # mu
             0.5, 0.5, 0.1) # sigma

# Wir brauchen zero und one inflated beta Verteilung??

theta.star0.3 = c(log(theta0[1:18]),
                  qlogis(theta0[19:21]),
                  theta0[22:24], 
                  log(theta0[25:27]))

data2 = data_aggr_mean
data2$height.fd[which(abs(data2$height.fd) > 20)] = NA
data2$angle[which(data2$angle == 1)] = .99

mllk(theta.star0, X = data2, N = 3)

t1 = Sys.time()
mod = nlm(f = mllk, p = theta.star0, X = data2, N = 3, print.level = 2, iterlim = 1000)
t2 = Sys.time()
t2-t1

mod2 = nlm(f = mllk, p = theta.star0.2, X = data2, N = 3, print.level = 2, iterlim = 1000)
mod3 = nlm(f = mllk, p = theta.star0.3, X = data2, N = 3, print.level = 2, iterlim = 1000)


# Checking for global optimum ---------------------------------------------

llk = rep(NA, 30)
mods = vector("list")
for (k in 1:30){
  r_theta0 = theta0.3 = c(runif(6, 0, 0.5),
                          c(9, 16) + runif(2,-5,5), runif(1,0,2),
                          c(1.5, 4) + runif(2, -1.5, 1.5), runif(1,0,2),
                          runif(1,5,10), runif(0,3), runif(0,2),
                          c(55, 40) + runif(2, -20, 20), runif(1,0,4),
                          runif(3,0,0.2),
                          c(1, -0.9, 0) + runif(3, -2,2),
                          runif(3,0,2))
  r_theta.star0 = c(log(theta0[1:18]),
                    qlogis(theta0[19:21]),
                    theta0[22:24], 
                    log(theta0[25:27]))
  mods[[k]] = nlm(mllk, r_theta.star0, X = data2, N = 3, iterlim = 300)
  llks[k] = -mods[[k]]$minimum
}

theta.star = mod$estimate
theta.star = mod2$estimate
theta.star = mod3$estimate

mllk(mod$estimate, data2, 3)
mllk(mod2$estimate, data2, 3)
mllk(mod3$estimate, data2, 3)

N = 3

Gamma = diag(N)
Gamma[!Gamma] = exp(theta.star[1:((N-1)*N)])
Gamma = Gamma/rowSums(Gamma)
delta = solve(t(diag(N)-Gamma+1),rep(1,N)) # stationary

# gamma distribution: Step length
mu.g = exp(theta.star[(N-1)*N+1:N]) # means of gamma distributions
sigma.g = exp(theta.star[(N-1)*N+N+1:N]) # sds of gamma distributions

# beta distribution: Turning angle/ pi
alpha = exp(theta.star[(N-1)*N+2*N+1:N]) # shape1 parameters of beta distributions
beta = exp(theta.star[(N-1)*N+3*N+1:N]) # shape2 parameters of beta distributions
p.b = plogis(theta.star[(N-1)*N+4*N+1:N])

# normal distribution: Height first difference
mu = theta.star[(N-1)*N+5*N+1:N] # means of normal distributions
sigma = exp(theta.star[(N-1)*N+6*N+1:N]) # sds of normal distributions



# Looking at results ------------------------------------------------------

# Steplength
hist(data_aggr_mean$step[c(16:200,2501:2546)][which(data_aggr_mean$label[c(16:200,2501:2546)] == 1)], prob = T, breaks = 20, xlab = "Step length", xlim = c(0,30), main = "Soaring")
curve(dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T)

hist(data_aggr_mean$step[c(16:200,2501:2546)][which(data_aggr_mean$label[c(16:200,2501:2546)] == 2)], prob = T, breaks = 20, xlab = "Step length", xlim = c(0,30), main = "Gliding")
curve(dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T)

hist(data_aggr_mean$step[c(16:200,2501:2546)][which(data_aggr_mean$label[c(16:200,2501:2546)] == 3)], prob = T, breaks = 10, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T)

par(mfrow = c(1,1))
# total
hist(data2$step, prob = T, breaks = 200, xlab = "Step length", xlim = c(0,30), main = "Resting")
curve(delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]), add = T, lwd = 2, col = color[3])

curve(
  delta[1]*dgamma(x, shape = mu.g[1]^2/sigma.g[1]^2, scale = sigma.g[1]^2/mu.g[1])+
    delta[2]*dgamma(x, shape = mu.g[2]^2/sigma.g[2]^2, scale = sigma.g[2]^2/mu.g[2])+ 
    delta[3]*dgamma(x, shape = mu.g[3]^2/sigma.g[3]^2, scale = sigma.g[3]^2/mu.g[3]),
  add = T, lty = "dashed", lwd = 2
)

# Angle
hist(data_aggr_mean$angle[c(16:200, 2501:2546)][which(data_aggr_mean$label[c(16:200, 2501:2546)] == 1)], prob = T, breaks = 15, xlab = "Angle", xlim = c(0,1), main = "Soaring")
curve(dbeta(x, shape1 = alpha[1], shape2 = beta[1]), add = T)

hist(data_aggr_mean$angle[c(16:200, 2501:2546)][which(data_aggr_mean$label[c(16:200, 2501:2546)] == 2)], prob = T, breaks = 30, xlab = "Angle", xlim = c(0,1), main = "Gliding")
curve(dbeta(x, shape1 = alpha[2], shape2 = beta[2]), add = T)

hist(data_aggr_mean$angle[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 3)], prob = T, breaks = 60, xlab = "Angle", xlim = c(0,1), main = "Resting")
curve(dbeta(x, shape1 = alpha[3], shape2 = beta[3]), add = T)

# total
hist(data2$angle, prob = T, breaks = 50, xlab = "Angle", xlim = c(0,1))
curve(delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]), add = T, lwd = 2, col = color[3])

curve(
  delta[1]*dbeta(x, shape1 = alpha[1], shape2 = beta[1])+
    delta[2]*dbeta(x, shape1 = alpha[2], shape2 = beta[2])+
    delta[3]*dbeta(x, shape1 = alpha[3], shape2 = beta[3]),
  add = T, lty = "dashed", lwd = 2
)

# Height fd
hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 1)], prob = T, breaks = 20, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Soaring")
curve(dnorm(x, mu[1], sigma[1]), add = T)

hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 2)], prob = T, breaks = 200, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Gliding")
curve(dnorm(x, mu[2], sigma[2]), add = T)

hist(data_aggr_mean$height.fd[c(16:250, 2501:2546)][which(data_aggr_mean$label[c(16:250, 2501:2546)] == 3)], prob = T, breaks = 50, xlab = "Height.fd", xlim = c(-2, 4.5), main = "Resting")
curve(dnorm(x, mu[3], sigma[3]), add = T)

# total
hist(data2$height.fd, prob = T, breaks = 100, xlab = "Height.fd")

lines(density(na.omit(data2$height.fd), bw = 0.03), lwd = 1)

curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1])
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2])
curve(delta[3]*dnorm(x, mu[3], sigma[3]), add = T, lwd = 2, col = color[3],n =1000)

curve(
  delta[1]*dnorm(x, mu[1], sigma[1])+
  delta[2]*dnorm(x, mu[2], sigma[2])+
  delta[3]*dnorm(x, mu[3], sigma[3]),
  add = T, lty = "dashed", lwd = 2, n=1000
)


# Getting probs -----------------------------------------------------------

states = viterbi(mod$estimate, data2, 3)
n = nrow(data2)
par(mfrow = c(1,1))

color = c("cornflowerblue", "orange", "forestgreen")

plot(data2$step[1:3000], pch = 20, col = color[states[1:3000]])
plot(data2$angle[1:3000], pch = 20, col = color[states[1:3000]])
plot(data2$height.fd[1:3000], pch = 20, col = color[states[1:3000]])












# LABELLING ---------------------------------------------------------------


2800*30
plot(ider$location.long[84001:84030], ider$location.lat[84001:84030], pch = 20)
plot(ider$height.above.msl[84001:84030], pch = 20) # 3

plot(ider$location.long[84031:84060], ider$location.lat[84031:84060], pch = 20)
plot(ider$height.above.msl[84031:84060], pch = 20) # 3

plot(ider$location.long[84061:84090], ider$location.lat[84061:84090], pch = 20)
plot(ider$height.above.msl[84061:84090], pch = 20) # 3

plot(ider$location.long[84091:84120], ider$location.lat[84091:84120], pch = 20)
plot(ider$height.above.msl[84091:84120], pch = 20) # 3

plot(ider$location.long[84121:84150], ider$location.lat[84121:84150], pch = 20)
plot(ider$height.above.msl[84121:84150], pch = 20) # 3

plot(ider$location.long[84151:84180], ider$location.lat[84151:84180], pch = 20)
plot(ider$height.above.msl[84151:84180], pch = 20) # 3

plot(ider$location.long[84181:84210], ider$location.lat[84181:84210], pch = 20)
plot(ider$height.above.msl[84181:84210], pch = 20) # 3

plot(ider$location.long[84211:84240], ider$location.lat[84211:84240], pch = 20)
plot(ider$height.above.msl[84211:84240], pch = 20) # 3

plot(ider$location.long[84241:84270], ider$location.lat[84241:84270], pch = 20)
plot(ider$height.above.msl[84241:84270], pch = 20) # 3

plot(ider$location.long[84271:84300], ider$location.lat[84271:84300], pch = 20)
plot(ider$height.above.msl[84271:84300], pch = 20) # 3

plot(ider$location.long[84301:84330], ider$location.lat[84301:84330], pch = 20)
plot(ider$height.above.msl[84301:84330], pch = 20) # 3

plot(ider$location.long[84331:84360], ider$location.lat[84331:84360], pch = 20)
plot(ider$height.above.msl[84331:84360], pch = 20) # 3

plot(ider$location.long[84361:84390], ider$location.lat[84361:84390], pch = 20)
plot(ider$height.above.msl[84361:84390], pch = 20) # 3
data$step[84361:84390]*1000

plot(ider$location.long[84391:84420], ider$location.lat[84391:84420], pch = 20)
plot(ider$height.above.msl[84391:84420], pch = 20) # 3

plot(ider$location.long[84421:84450], ider$location.lat[84421:84450], pch = 20)
plot(ider$height.above.msl[84421:84450], pch = 20) # flying / flapping flight?

plot(ider$location.long[84451:84480], ider$location.lat[84451:84480], pch = 20)
plot(ider$height.above.msl[84451:84480], pch = 20) #3

plot(ider$location.long[84481:84510], ider$location.lat[84481:84510], pch = 20)
plot(ider$height.above.msl[84481:84510], pch = 20) #3

plot(ider$location.long[84511:84540], ider$location.lat[84511:84540], pch = 20)
plot(ider$height.above.msl[84511:84540], pch = 20) #3

plot(ider$location.long[84541:84570], ider$location.lat[84541:84570], pch = 20)
plot(ider$height.above.msl[84541:84570], pch = 20) #3

plot(ider$location.long[84571:84600], ider$location.lat[84571:84600], pch = 20)
plot(ider$height.above.msl[84571:84600], pch = 20) #3

plot(ider$location.long[84601:84630], ider$location.lat[84601:84630], pch = 20)
plot(ider$height.above.msl[84601:84630], pch = 20) #3

plot(ider$location.long[84631:84660], ider$location.lat[84631:84660], pch = 20)
plot(ider$height.above.msl[84631:84660], pch = 20) #3

plot(ider$location.long[84661:84690], ider$location.lat[84661:84690], pch = 20)
plot(ider$height.above.msl[84661:84690], pch = 20) #3

plot(ider$location.long[84691:84720], ider$location.lat[84691:84720], pch = 20)
plot(ider$height.above.msl[84691:84720], pch = 20) #3






# Rest --------------------------------------------------------------------




simply_scatter(x = ider$location.long[3801:3830], y = ider$location.lat[3801:3830],
               z = ider$height.above.msl[3801:3830], colorvar = c(1, rep(0,29))) 
simply_scatter(x = ider$location.long[451:480], y = ider$location.lat[451:480], z = ider$height.above.msl[451:480])
simply_scatter(x = ider$location.long[511:540], y = ider$location.lat[511:540], z = ider$height.above.msl[511:540])
# upwards flight (soaring)
# maybe calculating turning angles in 3D is better

simply_scatter(x = ider$location.long[1111:1140], y = ider$location.lat[1111:1140], z = ider$height.above.msl[1111:1140], colorvar = c(1, rep(0, 29)))
# gliding down

simply_scatter(x = ider$location.long[2671:2700], y = ider$location.lat[2671:2700], z = ider$height.above.msl[2671:2700])
# basically no movement up or down 



hist(ider$height.fd[circle_ind], breaks = 1000, xlim = c(-20,20)) # Hauptsächlich hoch
hist(ider$height.fd[dir_ind], breaks = 1000, xlim = c(-20,20)) # Hauptsächlich runter
# auch hier Ausreißer raus
hist(ider$height.fd[chill_ind], breaks = 500, xlim = c(-20,20))




par(mfrow = c(3,1))
pacf(na.omit(d$angle[circle_ind]), main = "Pacf: Angle Circles")
pacf(na.omit(d$angle[dir_ind]), main = "Pacf: Angle Directed flight")
pacf(na.omit(d$angle[chill_ind]), main = "Pacf: Angle Chilling")
# looking very different

par(mfrow = c(1,1))
plot(d$angle[circle_ind], type = "l")
abline(v = seq(0,600,by = 30), col = "cornflowerblue")
acf(na.omit(d$angle[circle_ind]))

angle_test = d[ind2,]
angle_test$angle[which(angle_test$step == 0)] = 0

par(mfrow = c(3,1))
plot(d$angle[circle_ind], type = "h")
abline(v = seq(0,600,by = 30), col = "cornflowerblue")

plot(d$angle[dir_ind], type = "h")
abline(v = seq(0,600,by = 30), col = "cornflowerblue")

plot(d$angle[chill_ind], type = "h")
abline(v = seq(0,600,by = 30), col = "cornflowerblue")


hist(d$angle[circle_ind], prob = T, breaks = 100, xlim = c(-pi,pi), main = "Circles", xlab = "Angle")
hist(d$angle[dir_ind], prob = T, breaks = 100, xlim = c(-pi,pi), main = "Directed flight", xlab = "Angle")
helper = d$angle[chill_ind]
helper[which(d$step[chill_ind] == 0)] = 0
hist(helper, prob = T, breaks = 100, xlim = c(-pi,pi), main = "Chilling", xlab = "Angle")
curve(dvm(x, 0, 5), add = T, n = 1000, lwd = 2)

par(mfrow = c(1,1))
step = d$step*1000
hist(step[circle_ind], prob = T, breaks = 10000, xlim = c(0, 40), main = "Circles", xlab = "Step")
hist(step[dir_ind], prob = T, breaks = 10000, xlim = c(0,40), main = "Directed flight", xlab = "Step")
hist(step[chill_ind], prob = T, breaks = 10000, xlim = c(0,40), main = "Chilling", xlab = "Step")
curve(dgamma(x, shape = mu[3]^2/sigma[3]^2, scale = sigma[3]^2/mu[3]), add = T, n = 500, lwd = 2)

plot(ider$height.above.msl[circle_ind], main = "Cirlces: Rising up")
plot(ider$height.above.msl[dir_ind], main = "Directed: gliding down")
plot(ider$height.above.msl[chill_ind], main = "Chilling: Quiet random")


plot(ider$external.temperature[circle_ind])
plot(ider$external.temperature[dir_ind])
plot(ider$external.temperature[chill_ind])
# not very enlightening


plot(ider$landform[circle_ind])
plot(ider$landform[dir_ind])
plot(ider$landform[chill_ind])

# also very different!!




plot(d$x[451:480], d$y[451:480], pch = 20)
plot(d$angle[451:480], pch = 20)
plot(d$step[451:480]*1000, pch = 20) # Ausreißer stören
plot(d$step[451:479]*1000, pch = 20) # auch Periodisch
acf(na.omit(d$step[451:479]*1000))
dd = d[1:3000,]
dd$step = dd$step*1000
max(na.omit(dd$step))
dd$step[which(dd$step>50)] = NA # Ausreißer raus

for (i in 1:(99)){
  dd$step[30*i] = dd$angle[30*i] = NA
}
