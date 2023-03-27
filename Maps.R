
# Plot decoded states on map ----------------------------------------------
library(ggmap)
library(tidyverse)
library(ggplot2)
library(gridExtra)

color = c("deepskyblue", "orange", "springgreen3", "dodgerblue3")


# overall
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 4, lineend = "round")+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 1.5, alpha = 0.9, lineend = "round", linejoin = "mitre", colour = color[states])

# resting
data8 = data7[which(states == 1 | states == 4),]
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_point(data = data8, mapping = aes(x = x, y = y), colour = "deepskyblue4", size = 2) +
  geom_point(data = data8, mapping = aes(x = x, y = y), colour = "deepskyblue", size = 1)

# state 4
data9 = data7[which(states == 4),]
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_point(data = data9, mapping = aes(x = x, y = y), colour = "dodgerblue4", size = 2.5) +
  geom_point(data = data9, mapping = aes(x = x, y = y), colour = "dodgerblue", size = 1.5)


## Hinweg
data10 = data7 %>% filter(day <= 110)
first = ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_path(data = data10, aes(x = x, y = y), linewidth = 3, lineend = "round")+
  geom_path(data = data10, aes(x = x, y = y), linewidth = 1, alpha = 0.85, lineend = "round", linejoin = "mitre", colour = color[data10$states])+
  ggtitle("Autumn migration")

## Rückweg
data11 = data7 %>% filter(day >= 187)
second = ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_path(data = data11, aes(x = x, y = y), linewidth = 3, lineend = "round")+
  geom_path(data = data11, aes(x = x, y = y), linewidth = 1, alpha = 0.85, lineend = "round", linejoin = "mitre", colour = color[data11$states])+
  ggtitle("Spring migration")

grid.arrange(first, second, ncol=2)


## Map mit mTPI Werten

data7$mTPI2 = sign(data7$mTPI)*abs(data7$mTPI)^0.25

rbPal <- grDevices::colorRampPalette(c("yellow", "red"))
colors <- rbPal(nrow(data7))[as.numeric(cut(data7$mTPI2, breaks = nrow(data7)))]

ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 4, alpha = 0.7, lineend = "round", linejoin = "mitre", colour = colors)+
  geom_point(data = data8, mapping = aes(x = x, y = y), size = 1, colour = "black", alpha = 0.1, pch = 3)



rbPal <- grDevices::colorRampPalette(c("red", "yellow"))
colors <- rbPal(nrow(data8))[as.numeric(cut(data8$mTPI, breaks = nrow(data8)))]

ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_point(data = data8, mapping = aes(x = x, y = y), size = 1.5, colour = colors)


hist(data7$mTPI2, bor = "white")
hist(data7$mTPI)



library(sf)
library(elevatr)
library(sp)

prj_dd <- "EPSG:4326"

dev.off()

loc <- data.frame(x = runif(30, min = 70, max = 110), y = runif(30, min = 20, max = 50))
elevation_df <- get_elev_raster(loc, prj = prj_dd, z = 8)
raster::plot(elevation_df, maxpixels = 10000000, xlab = "longitude", ylab = "latitude", bty = "n")
points(data8$x, data8$y, pch = 20, col = "gray19")

