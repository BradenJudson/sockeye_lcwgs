library(tidyverse); library(ggplot2)
library(terra); library(geodata)
library(ggrepel)

setwd("~/sockeye_lcwgs/bioclim")

sites <- read.delim("../data/sk2023_sequenced.txt") %>%
  filter(Population != "Redfish_L") %>% 
  select(c("Population", "Latitude", "Longitude")) 

c.list <- list()

# For each site, extract bioclimatic data.
for (i in 1:nrow(sites)) {
  
  # Isolate site-specific lat and long.
  lon <- sites[i, "Longitude"]; lat <- sites[i, "Latitude"]
  
  # Download data from WorldClim at 5 arcsmin.
  # Store downloaded *.tiff files locally.
  # First, download current (1970 - 2000) conditions.
  c.clim <- worldclim_tile(var  = "bio",
                           lon  = lon,
                           lat  = lat,
                           res  = 1/2,
                           path = "./wcdata",
                           download = FALSE)
  
  # Download projected bioclimatic variables.
  # Following Tigano et al. 2023 [Evo. App.]
  # Climate data at RCP2.6 (~best case).
  # Set download = TRUE for first time.
  f.clim26 <- cmip6_tile(ssp  = "126",
                         lon  = lon,
                         lat  = lat,
                         res  = 5,
                         path = "./fclim26",
                         time = "2041-2060",
                         var  = "bioc",
                         model = "UKESM1-0-LL",
                         download = FALSE)
  
  # Similar to above, but RCP8.5 (~worst case).
  f.clim85 <- cmip6_tile(ssp  = "585",
                         lon  = lon,
                         lat  = lat,
                         res  = 5,
                         path = "./fclim85",
                         time = "2041-2060",
                         var  = "bioc",
                         model = "UKESM1-0-LL",
                         download = FALSE)
  
  # Extract site coordinates and project.
  points <- vect(sites[,c(2:3)], crs = "EPSG:4326",
                 geom = c("Longitude", "Latitude"))
  
  # Set up a function to extract data from each of the 
  # three climate datasets isolated above.
  # Also renames columns consistently and adds two qualifiers.
  f <- function(clim, ssp, time) cbind(sites[i, 1],
                                 terra::extract(clim, y = points)[i, 2:20]) %>% 
    `colnames<-`(., c("Site", paste0("bio", seq(1, 19, 1)))) %>% 
    mutate(period = as.factor(time), ssp = as.factor(ssp)) 
  
  # For each site, bind together and specify time period and ssp.
  # NOTE: First (1970s-2000) data set is in a different order!!!
  # bio1, bio10, bio11,...bio9. Unlike the other two. 
  # Arrange separately otherwise there is a severe mismatch.
  c.list[[i]] <- rbind(
    f(f.clim26, time = "2041-2060", ssp = "26"),
    f(f.clim85, time = "2041-2060", ssp = "85"),
    cbind(sites[i, 1], terra::extract(c.clim, points)[i, 2:20]) %>% 
      `colnames<-`(., c("Site", paste0("bio", c(1, seq(10, 19, 1),
                                                seq(2, 9, 1))))) %>% 
      mutate(period = "1970-2000", ssp = NA),
    make.row.names = FALSE # Just use numbers.
  )
}

# Put the bioclim data into a single dataframe. 
biovar <- bind_rows(c.list); head(biovar)

# No missing data.
sum(is.na(biovar)) == nrow(biovar)/3

write.csv(biovar, "../data/sk_bioclim.csv", row.names = FALSE)



# PCA --------------------------------------------------------------------------

"%ni%" <- Negate("%in%") # Custom not-in operator.
# Retain only contemporary data and remove unnecessary columns.
cbio <- biovar[biovar$period == "1970-2000", colnames(biovar) %ni% c("period", "ssp")]

# Correlation matrix of bioclimatic variables. 
(bcor <- abs(cor(cbio[,-1]))) # Excludes site info.
bcor[!lower.tri(bcor)] <- 0; hist(bcor)

# Retain variables with less than 80% correlation.
(ev <- cbio[, !apply(bcor, 2, function(x) any(x > 0.80))])
dim(ev) # number of variables retained

# Scale and center (Z-transform) uncorrelated variables.
sval <- as.data.frame(scale(ev)) %>%  # Make numeric too.
  mutate_if(is.character, is.numeric)

# Make a population-specific dataframe with scores along various PC axes. 
(evpca <- prcomp(sval))
(pc1ve <- summary(evpca)$importance[2,1]) # PC1 % variation explained.
(pc2ve <- summary(evpca)$importance[2,2]) # PC2 % variation explained.

(loadings <- as.data.frame(evpca$rotation))

# Retain PC scores and attach to each site.
pop_climPC <- as.data.frame(prcomp(sval)$x) %>% 
  mutate(pop = tools::toTitleCase(tolower(gsub("\\_.*", "", cbio$Site))),
         Latitude = as.numeric(sites$Latitude))

# Lower Thumb gets incorrectly parsed to "l".
pop_climPC[pop_climPC == "l"] <- "Thumb"

# Plot PC1 and 2 scores.
ggplot(data = pop_climPC,
       aes(x = PC1, y = PC2,
           fill = Latitude)) +
  theme_bw() + 
  geom_segment(data = loadings, 
               # Expand 7-fold purely for visualization. 
               aes(x = 0, xend = PC1 * 5, 
                   y = 0, yend = PC2 * 5,
                   color = rownames(loadings)),
               arrow = arrow(length = unit(0.3, "cm")),
               inherit.aes = FALSE) +
  geom_point(size = 2, shape = 21) +  
  scale_fill_gradient(low = "#F0FFFF", high = "#1874CD") +
  labs(y = paste0("PC2 (", round(pc2ve*100, 1), "%)"),
       x = paste0("PC1 (", round(pc1ve*100, 1), "%)")) +
  geom_label_repel(aes(label = pop), 
                   max.overlaps = Inf,
                   min.segment.length = 0,
                   box.padding = 1/3,
                   size = 1.7,
                   segment.size = 1/5,
                   fill = "white") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        legend.margin = margin(t = 0, unit = "cm")) +
  guides(fill = guide_colourbar(barwidth = 25, label = F, reverse = T,
                                ticks = F, frame.colour = "black"),
         color = guide_legend(nrow = 1)) +
  coord_cartesian(ylim = c(-3.5, 2.1), clip = "off") +
  annotate("text", label = "North", x = -4.9, y = 3.18, hjust = 1) +
  annotate("text", label = "South", x = 2.0, y = 3.18, hjust = -1)





   

ggsave("../plots/envPCA.tiff", dpi = 300, width = 7, height = 7)






