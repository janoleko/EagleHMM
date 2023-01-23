
# Plot decoded states on map ----------------------------------------------

library(ggmap)
ggmap(get_map(location = c(75, 20, 110, 50), map_type = "terrain", source = "stamen"))

