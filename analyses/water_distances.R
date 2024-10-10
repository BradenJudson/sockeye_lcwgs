setwd("~/sockeye_lcwgs/analyses")

library(sf); library(gdistance); library(ggplot2); library(tidyverse)
library(raster); library(rmapshaper); library(igraph); library(vegan)

# Read in some large objects initially computed here:
# https://github.com/BradenJudson/chinook_hdGBS_lcWGS/blob/main/map/distance_matrix.R

# Read in waterway transition matrix.
tr1 <- readRDS("../data/waterway_transition_matrix.rds")

# Read in raster.
rast <- raster("../data/invraster.grd")
rr <- as.data.frame(rast$layer, xy = TRUE) %>% 
  filter(!is.na(layer))

# Population locations.
sites <- read.delim("../data/sk2023_sequenced.txt")
sites[sites$site_full == "Thumb (Lower)", c(10,9)] <- c(-154.4511, 57.5776)
sites[sites$site_full == "Klukshu", c(10,9)] <- c(-137.0018, 60.12)
sites[sites$site_full == "Kuthai", c(10,9)] <- c(-133.2336, 59.10)
sites[sites$site_full == "Kutlaku", c(10,9)] <- c(-134.200, 56.6048)
sites[sites$site_full == "Copper", c(10,9)] <- c(-131.7700, 53.16082)
sites[sites$site_full == "Mercer", c(10,9)] <- c(-132.9004, 53.5928)
sites[sites$site_full == "Bloomfield", c(10,9)] <- c(-128.71, 52.85661)
sites[sites$site_full == "Hetta", c(10,9)] <- c(-132.5943, 55.17258)
sites[sites$site_full == "McDonnell", c(10,9)] <- c(-127.9, 54.52)
sites[sites$site_full == "Sinta", c(10,9)] <- c(-125.74, 55.30229)
sites[sites$site_full == "Ashlulm", c(10,9)] <- c(-127.41, 51.68916)
sites[sites$site_full == "Canoe", c(10,9)] <- c(-127.0153, 51.4)
sites[sites$site_full == "Cayenne", c(10,9)] <- c(-119.415, 51.32071)
sites[sites$site_full == "Henderson", c(10,9)] <- c(-125.043, 49)
sites[sites$site_full == "Big Creek", c(10,9)] <- c(-123.2710, 46.5387)
sites[sites$site_full == "Azuklotz", c(10,9)] <- c(-126.7456, 56.0115)
sites[sites$site_full == "King Salmon", c(10,9)] <- c(-132.937, 58.7300)
sites[sites$site_full == "Koeye", c(10,9)] <- c(-127.8825,51.78039)

(q <- ggplot(NULL) + geom_raster(data = rr, aes(x=x,y=y), fill = "grey") +
  geom_point(data = sites, aes(x= Longitude,y=Latitude)) +
    scale_x_continuous(limits = c(-128, -127)) +
    scale_y_continuous(limits = c(51,52)))

sites <- sites[!sites$site_full %in% c("Redfish"), ]
# Convert sites to an object of class SpatialPointsDataFrame.
sitepoints <- SpatialPoints(coords = sites[,c("Longitude", "Latitude")],
                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))


spLine_list <- list()
# Nested for-loop that isolates unique and uni-directional combinations of populations.
for (j in 1:nrow(sitepoints@coords)) {
  
  for (k in j:nrow(sitepoints@coords)) {
    
    if(k == j) { next } # Skip self-comparisons (i.e., distance = 0).
    
    # Display loop progres over pairs of populations.
    pops <- sites$site_full
    print(paste0("Calculating distance from ", pops[j], " to ", pops[k]))
    
    # Shortest path between sites.
    SP <- shortestPath(tr1, sitepoints[j], sitepoints[k], output = "SpatialLines")
    
    # Fortify lines object for plotting purposes. 
    # # Also compute and add inter-population distances.
    SPDF <- fortify(SpatialLinesDataFrame(SP, data = data.frame(ID = 1))) %>%
      mutate(pop1 = sites[j, "site_full"], pop2 = sites[k, "site_full"],
             distance_m = round(geosphere::lengthLine(SP), 0))
    
    # Subset dataframe above and poplate list.
    spLine_list[[length(spLine_list)+1]] <- SPDF
    
  }
}

# Reformat list of spatial lines into a dataframe. 
# Make a 'label' for plotting purposes.
site_lines <- do.call("rbind", spLine_list) %>% 
  mutate(pair  = as.factor(paste(pop1, "-", pop2)),
         label = as.factor(paste0(pair, ": ", 
                           format(distance_m, format = "d", 
                                  big.mark = ","), " (m)")))

dists <- site_lines %>% group_by(pair) %>% 
  sample_n(1) %>% .[,c(7:10)] %>% arrange(pop1, pop2)

dist_mat <- as.data.frame(get.adjacency(graph.data.frame(dists), attr = "distance_m", sparse = F))
dist_mat[c(22:36,38:42),42] <- as.numeric(dist_mat[42, c(22:36,38:42)]); dist_mat[42, c(22:36,38:42)] <- 0
dist_mat[4:36, 37] <- as.numeric(dist_mat[37, 4:36]); dist_mat[37, 4:36] <- 0
tdistmat <- t(dist_mat)
write.csv(dist_mat, "../data/Oner_distances_mat.csv")


# PCNMs + visualization --------------------------------------------------------

dist_mat <- t(read.csv("../data/Oner_distances_mat.csv", row.names = 1))
pcnms <- as.data.frame(pcnm(dis = as.dist(dist_mat))$vector)

cd <- cbind(sites %>% arrange(site_full), 
            pcnms %>% arrange(rownames(.))) %>% 
  pivot_longer(cols = c("PCNM1", "PCNM2", "PCNM3", "PCNM4"))

ggplot(NULL) +
  geom_raster(data = as.data.frame(rast$layer, xy = TRUE) %>% 
                filter(!is.na(layer)),
              aes(x = x, y = y), fill = "gray80") +
  geom_point(data = cd, aes(x = Longitude, y = Latitude, fill = value),
             shape = 21, size = 2, colour = "black") +
  scale_fill_continuous(high = "navy", low = "lightskyblue1") +
  facet_wrap(~name, ncol = 2, nrow = 2) +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = "none") +
  labs(x = "Longitude", y = "Latitude") +
  scale_y_continuous(limits = c(45, 62), expand = c(0,0)) +
  scale_x_continuous(limits = c(-155, -116), expand = c(0,0))

ggsave("../plots/PCNMs.tiff", width = 14, height = 10, dpi = 300)

