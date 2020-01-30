library(tidyverse)
library(janitor)
map <- read_csv("map_vec")

btl <- read_delim("Soclim/data/data_btl_JULIA.csv", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
btl <- clean_names(btl)

names(btl) <- c('station', 'date', 'lon', 'lat', 'depth', 'temp', 'sal', 'ox_mml', 'chl_fluo', 'bbp', 'cp', 'fluo_bbp', 'fluo_cp', 'bbp_cp', 'poc_mmol_m3', 'pon_mmol_m3', 'station_bis')

ggplot(btl)+
  geom_polygon(aes( x = long, y = lat, group = group), data = map)+
  geom_point(aes(x = lon, y = lat, colour = chl_fluo))+
  coord_quickmap()+
  scale_color_viridis_c()+
  xlim(50,75)+
  ylim(-65, -25)+
  theme_bw()


cyto <- read_delim("Soclim/data/data_cyto.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
cyto <- clean_names(cyto)

micro <- read_delim("Soclim/data/data_micro.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
micro <- clean_names(micro)


