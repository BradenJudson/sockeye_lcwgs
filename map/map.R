setwd("~/sockeye_lcwgs/map")

# devtools::install_github("rspatial/geodata")

library(ggplot2); library(tidyverse); library(raster); library(sf)
library(bcmaps); library(ggspatial); library(sp); library(geodata)
library(ggrepel); library(rnaturalearth); library(cowplot)

# Custom "not in" operator.
"%ni%" <- Negate("%in%")

# Establish base map -----------------------------------------------------------

# Download American shape data and a subset of the Canadian data.
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))
prov <- sf::st_as_sf(geodata::gadm(country = "CAN", level = 1, path = "."))
bcn <- prov[prov$NAME_1 %in% c("Alberta", "Yukon", "Northwest Territories", "Nunavut"),]

# Use bcmaps package to get high resolution BC shape data.
# Necessary due to complexity of coastline.
bch <- st_transform(bc_bound_hres(), crs = 4326)

# Read in site information and reformat labels.
sites <- read.delim(file = "../data/sk2023_sequenced.txt") %>% 
  arrange(Latitude) %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population))),
         sitenum = as.numeric(rownames(.))) 

sites[sites$Population == "L_Thumb", "site"] <- "Thumb"
sites[sites$Population == "Big_Silver_Cr", "site"] <- "Big Silver"

(riv <- bcmaps::watercourses_5M())
rivers <- ne_load(scale = 10, type = "rivers_lake_centerlines", destdir = "maps", returnclass = "sf") 

# Above command misses the Canadian Okanagan, so I manually read that in.
okanagan <- st_read("ok_path.kml")      

(yv <- 55 - (55 - 44)*(1:ceiling((nrow(sites)/2)))/(nrow(sites))*2)

(pnw <- ggplot(data = sites) +
    geom_sf(data = USA, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = bcn, fill = "gray90", linewidth = 1/10) +
    geom_sf(data = bch, fill = "gray90", linewidth = 1/10) +
    ggspatial::annotation_north_arrow(location = "bl",
                                      pad_x = unit(3.5, "cm"),
                                      style = ggspatial::north_arrow_fancy_orienteering()) +
    ggspatial::annotation_scale(location = "bl", pad_x = unit(3.75, "cm"),
                                pad_y = unit(1/10, "cm"),
                                width_hint = 1/10) + 
    geom_sf(data = okanagan, color  = "skyblue", linewidth = 1/4) +
    geom_sf(data = riv,      colour = "skyblue", linewidth = 1/4) +
    geom_sf(data = rivers,   colour = "skyblue", linewidth = 1/4) +  
    # geom_sf(data = lakes,    colour = "skyblue", fill = "skyblue") +
    geom_point(data = sites, size = 2.5, stroke = 1/3,
               shape = 21, color = "black", fill = "white",
               aes(x = Longitude, y = Latitude)) +
    geom_text(data = sites, size = 1.3, fontface = "bold",
              aes(x = Longitude, y = Latitude, label = `sitenum`)) +
    scale_fill_manual(values = rep("white", nrow(sites)),
                      labels = paste(sites$sitenum, sites$site)) +
    guides(fill = guide_legend(override.aes = list(alpha = 0))) +
    coord_sf(xlim = c(-115, -155), ylim = c(44, 62)) +
    theme_minimal() +
    theme(legend.position = "right", panel.grid = element_blank(), 
          panel.background = element_rect(fill = alpha("skyblue", 1/10)),
          panel.border = element_rect(color = "black", fill = NA),
          plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_rect(aes(ymin = 43.5, ymax = 54.8, xmax = -145.5, xmin = -156.5),
              colour = "black", fill = "white", alpha = 1/6, size = 1/20) +
    geom_text(data = sites[sites$sitenum <= 25, ], aes(label = paste(sitenum, ". ", site), 
                                                       x = -156, y = yv, hjust = 0), size = 2, inherit.aes = F) +
    geom_text(data = sites[sites$sitenum > 25, ], aes(label = paste(sitenum, ". ", site), 
                                                      x = -151, y = yv[1:24], hjust = 0), size = 2, inherit.aes = F))

ggsave("plots/bc_map.tiff", dpi = 300, width = 6, height = 6)


# Inset -------------------------------------------------------------------


ca <- map_data("world", "Canada")
us <- map_data("world", "USA") 
me <- map_data("world", "Mexico")

(ins <- ggplot() +
    geom_polygon(data = us, aes(x = long, y = lat, group = group),
                 fill = "grey80", colour = "black", linewidth = 1/8) +
    geom_polygon(data = ca, aes(x = long, y = lat, group = group),
                 fill = "grey90", colour = "black", linewidth = 1/8) +
    geom_polygon(data = me, aes(x = long, y = lat, group = group),
                 fill = "grey70", colour = "black", linewidth = 1/8) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", 
                                      fill = NA, linewidth = 1/4),
          panel.background = element_rect(fill = "white"))  +
    annotate("rect", fill = NA, colour = "black",
             linewidth = 1/2,
             xmin = -155, xmax = -115,
             ymin = 44, ymax = 62) +
    # Important to maintain accurate proportions/orientations. 
    # Plot is cartesian otherwise and appears distorted.
    coord_map(ylim = c(72, 20),
              xlim = c(-57, -165)))



# Add inset --------------------------------------------------------------------


ggdraw(plot = pnw) +
  draw_plot({
    ins
  },
  x = 0.79,
  y = 0.53,
  width = 0.2,
  height = 0.5)

ggsave("plots/map_winset.tiff", dpi = 300, 
       width = 6, height = 6)
