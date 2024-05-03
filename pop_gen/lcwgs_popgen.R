setwd("~/sockeye_lcwgs/pop_gen")

library(ggplot2); library(tidyverse); library(ggrepel)

# Read in covriance matrix.
covmat <- read.table("../data/pca_full.cov")

# Read in file information.
files <- read.table("../data/bam_list_n379.txt", col.names = c("file")) %>% 
    mutate(fish_ID = gsub(".dedup.clip.rg.bam", "", 
                gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file)),
           pop = substr(x = fish_ID, start = 0, stop = 3)) 

# Slight adjustments to strings. 
files[files$fish_ID == "03_alignments_dedup_clip/NS.LH00487_0010.004.IDT_i7_223_58---IDT_i5_223_58.Was-2002-02-14.merged", "fish_ID"] <- "Was-2002-02-14"
files[files$pop == "03_", "pop"] <- "Was"

# Site information.
sites <- read.delim("../data/sk2023_sequenced.txt")[,c(4,5,8:9L)] %>% 
  mutate(pop = substr(x = Population, start = 0, stop = 3))
sites[sites$pop == "L_T", "pop"] <- "Thu"

# Individual sample information. 
samples <- read.csv("../data/LandGen-Sockeye-Ext-Data.csv")[,c(1,6)] %>% 
  group_by(fish_ID) %>% sample_n(1) %>% 
  mutate(site_full = gsub("Shushwap", "Shuswap", site_full))
# 
nrow(files[!files$fish_ID %in% samples$fish_ID, ]) == 0

# Compute and extract eigenvectors from covariance matrix.
mmepca <- eigen(covmat)
eigenvectors <- mmepca$vectors

# Put eigenvector values, and individual/population information together.
pcavectors <- as_tibble(cbind(files, data.frame(eigenvectors))) %>% 
  group_by(fish_ID) %>% sample_n(1) %>% # Randomly choose one of duplicate samples (n = 2).
  merge(., samples) %>% 
  select(c(1:5, ncol(.)))

# Get population-level average PC values for visualization.
pc_popav <- pcavectors %>% 
  group_by(site_full) %>% 
  summarise(X1A = mean(X1),
            X2A = mean(X2))

# Add pop-averages to main dataframe.
pcavectors2 <- merge(pcavectors, pc_popav, by = "site_full", all.x = T) %>% 
  merge(., sites, by = "site_full")

pca.eigenval.sum <- sum(mmepca$values)
varPC1 <-  (mmepca$values[1]/pca.eigenval.sum)*100
varPC2 <-  (mmepca$values[2]/pca.eigenval.sum)*100

ggplot(data = pcavectors2) +
  geom_segment(aes(x = X1A, y = X2A, xend = X1, yend = X2, 
                   colour = Latitude), linewidth = 3/4) +
  geom_point(aes(x = X1, y = X2, fill = Latitude),
             shape =21, size= 2) + theme_bw() +
  ggrepel::geom_label_repel(data = pc_popav, max.overlaps = Inf, 
                           min.segment.length = 0, box.padding = 3/4,
                            aes(x = X1A, y = X2A, label = site_full)) + 
  theme(legend.position = "top", legend.title = element_blank()) +
  labs(x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) +
  guides(fill = guide_colourbar(barwidth = 25, label = F, reverse = T, 
                ticks = F, frame.colour = "black"), colour = "none") +
  scale_y_continuous(breaks = seq(-1/2, 1/2, 0.05)) +
  coord_cartesian(ylim = c(-0.065, 0.08), clip = "off") +
  annotate("text", label = "N", x = -0.058, y = 0.097, size = 7) +
  annotate("text", label = "S", x = 0.053,  y = 0.097, size = 7)
  
ggsave("../plots/pca_6MSNPs_NS.tiff", dpi = 300, height = 8, width = 10)
