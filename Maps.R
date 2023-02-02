
# Plot decoded states on map ----------------------------------------------
library(ggmap)

# overall
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 4, lineend = "round")+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 1.5, alpha = 0.8, lineend = "round", linejoin = "mitre", colour = color[states])

# resting
data8 = data7[which(states == 1),]
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  geom_point(data = data8, mapping = aes(x = x, y = y), colour = "deepskyblue4", size = 2.5) +
  geom_point(data = data8, mapping = aes(x = x, y = y), colour = "deepskyblue", size = 1.5) +
  labs(title = "Resting locations")

# state 4
data9 = data7[which(states == 4),]
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  geom_point(data = data9, mapping = aes(x = x, y = y), colour = "dodgerblue4", size = 2.5) +
  geom_point(data = data9, mapping = aes(x = x, y = y), colour = "dodgerblue", size = 1.5) +
  labs(title = "State 4 locations")

