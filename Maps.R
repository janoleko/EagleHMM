
# Plot decoded states on map ----------------------------------------------
library(ggmap)

# overall
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 4, lineend = "round")+
  geom_path(data = data7, aes(x = x, y = y), linewidth = 1.5, alpha = 0.9, lineend = "round", linejoin = "mitre", colour = color[states])

# resting
data8 = data7[which(states == 1 | states == 4),]
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_point(data = data8, mapping = aes(x = x, y = y), colour = "deepskyblue4", size = 2.5) +
  geom_point(data = data8, mapping = aes(x = x, y = y), colour = "deepskyblue", size = 1.5)

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
  geom_path(data = data10, aes(x = x, y = y), linewidth = 4, lineend = "round")+
  geom_path(data = data10, aes(x = x, y = y), linewidth = 1.5, alpha = 0.9, lineend = "round", linejoin = "mitre", colour = color[data10$states])

## Rückweg
data11 = data7 %>% filter(day >= 187)
second = ggmap(get_map(location = c(75, 20, 110, 50), map_type = "toner", source = "stamen"))+
  annotate("rect", xmin = 75, xmax = 110, ymin = 20, ymax = 50, fill = "white", alpha = 0.4)+
  geom_path(data = data11, aes(x = x, y = y), linewidth = 4, lineend = "round")+
  geom_path(data = data11, aes(x = x, y = y), linewidth = 1.5, alpha = 0.9, lineend = "round", linejoin = "mitre", colour = color[data11$states])

grid.arrange(first, second, ncol=2)

