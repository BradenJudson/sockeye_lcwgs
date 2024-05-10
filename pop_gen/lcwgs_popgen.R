setwd("~/sockeye_lcwgs/pop_gen")

library(ggplot2); library(tidyverse); library(ggrepel)


# PCA: Full --------------------------------------------------------------------

 # Write a function for doing the PCA and organizing sample information. 
pca_func <- function(covmat, file_list) {
  
  # Read in covriance matrix.
  covmat <- read.table(covmat)
  
  # Read in file information.
  files <- read.table(file_list, col.names = c("file")) %>% 
    mutate(fish_ID = gsub(".dedup.clip.rg.bam", "", 
                          gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 
  
  # Site information.
  sites <- read.delim("../data/sk2023_sequenced.txt")[,c(4,5,8:9L)] %>% 
    mutate(pop = substr(x = Population, start = 0, stop = 3))
  sites[sites$pop == "L_T", "pop"] <- "Thu"
  
  # Individual sample information. 
  samples <- read.csv("../data/LandGen-Sockeye-Ext-Data.csv")[,c(1,6)] %>% 
    group_by(fish_ID) %>% sample_n(1) %>% 
    mutate(site_full = gsub("Shushwap", "Shuswap", site_full))
  
  # Compute and extract eigenvectors from covariance matrix.
  mmepca <- eigen(covmat)
  eigenvectors <- mmepca$vectors
  
  # Put eigenvector values, and individual/population information together.
  pcavectors <- as_tibble(cbind(files, data.frame(eigenvectors))) %>% 
    group_by(fish_ID) %>% 
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
  
  # Proportions of variance explained by each axis.
  pca.eigenval.sum <- sum(mmepca$values)
  varPC1 <-  (mmepca$values[1]/pca.eigenval.sum)*100
  varPC2 <-  (mmepca$values[2]/pca.eigenval.sum)*100
  
  # For automatic labeling of fill gradien bar.
  listlats <- list(North = max, South = min)
  listvals <- lapply(listlats, function(x) x(pcavectors2$Latitude))
  
  
  (p <- ggplot(data = pcavectors2, aes(fill = Latitude, colour = Latitude)) +
    geom_segment(aes(x = X1A, y = X2A, xend = X1, yend = X2), linewidth = 3/4) +
    geom_point(aes(x = X1, y = X2),color = "black", shape = 21, size = 2) + theme_bw() +
    scale_color_gradient(guide = 'none') +
    scale_fill_gradient(breaks = unlist(listvals)) +
    ggrepel::geom_label_repel(data = pc_popav, max.overlaps = Inf, 
                               min.segment.length = 0, box.padding = 3/4,
                               aes(x = X1A, y = X2A, label = site_full), 
                              inherit.aes = FALSE) +
    theme(legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 11)) +
    labs(x = paste0("PC1 (", round(varPC1, 1), "%)"),
          y = paste0("PC2 (", round(varPC2, 1), "%)")) +
    guides(fill = guide_colorbar(barwidth = 25, label = T, reverse = T, 
                                   ticks = F, frame.colour = "black", label.position = "top")))
    
}

(full_pca <- pca_func(covmat = "../data/pca_full.cov",
          file_list = "../data/bam_list_n377.txt"))

ggsave("../plots/pca_6MSNPs_NS.tiff", dpi = 300, height = 8, width = 10)


(full_maf1 <- pca_func(covmat = "../data/pca_full_maf1.cov",
                       file_list = "../data/bam_list_n377.txt"))

ggsave("../plots/pca_6MSNPs_maf1.tiff", dpi = 300, height = 8, width = 10)


# PCA: Fraser ------------------------------------------------------------------

(fra_pca <- pca_func(covmat = "../data/pca_fraser.cov",
                     file_list = "../data/bam_list_fraser.txt"))

ggsave("../plots/pca_fraser_maf5.tiff", dpi = 300, height = 8, width = 10)

(fra_pca <- pca_func(covmat = "../data/pca_fraser_maf1.cov",
                     file_list = "../data/bam_list_fraser.txt"))

ggsave("../plots/pca_fraser_maf1.tiff", dpi = 300, height = 8, width = 10)


# FST --------------------------------------------------------------------------


