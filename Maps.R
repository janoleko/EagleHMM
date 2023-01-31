
# Plot decoded states on map ----------------------------------------------

library(ggmap)
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  geom_point(data = data7, mapping = aes(x = x, y = y, colour = color[states]))

data8 = data7[which(states == 1),]

ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  geom_point(data = data8, mapping = aes(x = x, y = y, colour = "skyblue"))

data9 = data7[which(states == 4),]

ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  geom_point(data = data9, mapping = aes(x = x, y = y))
