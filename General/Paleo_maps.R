### Make paleogeographic maps. 

install.packages("gplatesr")

coastline_gws_url <-
  "http://gws.gplates.org/reconstruct/coastlines/?time=0&model=GOLONKA"
polygons_gws_url <-
  "http://gws.gplates.org/reconstruct/static_polygons/?time=0&model=GOLONKA"

toarcian_coastlines <-
  rgdal::readOGR(coastline_gws_url)
toarcian_polygons <-
  rgdal::readOGR(polygons_gws_url)

library(broom)
library(ggplot2)
library(ggthemes)

toarcian_coastlines <-
  broom::tidy(toarcian_coastlines)
toarcian_polygons <-
  broom::tidy(toarcian_polygons)

toarcian_map <-
  ggplot() +
  geom_map(
    data = toarcian_polygons, map = toarcian_polygons,
    aes(x = long, y = lat, map_id = id),
    size = 0.15, fill = "#d8d8d8"
  ) +
  geom_map(
    data = toarcian_coastlines, map = toarcian_coastlines,
    aes(x = long, y = lat, map_id = id),
    size = 0.15, fill = NA, colour = "grey30"
  ) +
  geom_rect(
    data = data.frame(xmin = -180, xmax = 180, ymin = -90, ymax = 90),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = 1, fill = NA, size = 0.3
  ) +
  coord_map("mollweide") +
  ggthemes::theme_map()

toarcian_map +
  labs(
    title = " "
  )