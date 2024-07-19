setwd("~/sockeye_lcwgs/analyses")

library(vegan); library(tidyverse); library(adespatial)
library(sf)


# Site data --------------------------------------------------------------------

sites <- read.delim("../data/sk2023_sequenced.txt") %>% 
  mutate(site = tools::toTitleCase(tolower(gsub("\\_.*", "", Population)))) %>% 
  .[, c((ncol(.) -3):ncol(.))] %>% arrange(site) %>% filter(site != 'Redfish')

# Compute a linear distance matrix accounting for the curvature of the earth.
dm <- st_distance(st_as_sf(sites[,c(4,3,2)], coords = 2:3, crs = 4326))
hist(dm, main = NULL, xlab = "Distance between sites")

# Convert distance matrix into distance-based Moran's eigenvector maps.
db <- create.dbMEM.model(coord = sites[,c(2:3)], nsites = nrow(sites)) %>% 
  as.data.frame() %>% arrange(sites) %>% `rownames<-`(., c(sites$site))  


# Bioclimatic variables --------------------------------------------------------

bioclim <- read.csv("../data/sk_bioclim.csv") %>% 
  filter(is.na(ssp)) %>% .[,1:(ncol(.)-2)] %>% 
  mutate(Site = tools::toTitleCase(tolower(gsub("\\_.*", "", Site)))) %>% 
  arrange(Site) %>% column_to_rownames("Site") 

# Make sure ordering is consistent.
rownames(db) == rownames(bioclim)

png("../plots/predictor_cors.png", width = 1200, height = 1000, res = 100)
bioclim_sub <- as.matrix(cbind(subset(bioclim, select = c(bio2, bio4,
                                               bio8, bio10, bio15, bio16)),
                               sites[,c("Latitude", "Longitude")]))
gplots::heatmap.2(cor(bioclim_sub)^2, trace = 'none'); dev.off()


# RDA --------------------------------------------------------------------------

# Function for reading and formatting the output of vcftools --freq.
rf <- function(dir) paste0(dir, list.files(pattern = ".*.frq", 
                                           path = dir)) %>% set_names(.) %>% 
  map_df(read_table, .id = "FileName", 
         col_names = c("CHROM", "POS", "N_ALLELES", 
                       "N_CHR", "MajAll", "MinAll")) %>% 
  filter(MajAll != "{ALLELE:FREQ}") %>% 
  mutate(Pop  = gsub("\\_.*","", str_sub(FileName, start = 17)),
         MajAF = as.numeric(sub(".*:", "", MajAll)),
         MajNuc  = as.factor(substr(MajAll, start = 0, stop = 1)),
         MinAF = as.numeric(sub(".*:", "", MinAll)),
         MinNuc = as.factor(substr(MinAll, start = 0, stop = 1)),
         Gpos   = paste0(CHROM, "_", POS)) %>% 
  select(c("Gpos", "MajAF")) %>% 
  pivot_wider(values_from = MajAF, names_from = Gpos)

temp  <- paste0("../data/pop_frq/non_pruned/", list.files(path = "../data/pop_frq/non_pruned"))
lcafs <- as.data.frame(do.call(rbind, lapply(temp, rf))) %>% 
  `rownames<-`(., gsub("\\_.*","", str_sub(temp, start = 28)))

rownames(lcafs)

# Use below if starting fresh and don't want to wait for above to process.
# write.csv(lcafs, "../data/lcafs.csv", row.names = T)
# lcafs <- read.csv("../data/lcafs.csv")

nmod <- rda(lcafs ~ 1, db)
fmod <- rda(lcafs ~ ., db)
optdb <- ordiR2step(nmod, scope = formula(fmod), direction = "both")
summary(optdb)$call # Returns optimal dbMEMs to use.

(lcRDA <- vegan::rda(lcafs ~ bio2 + bio4 + bio8 + bio10 + bio15 + bio16 + Condition(dbMEM.1 + dbMEM.2 + dbMEM.8),
                     data = cbind(bioclim, db), scale = TRUE))
RsquareAdj(lcRDA) 
(eigvsumm <- round(summary(eigenvals(lcRDA, model = "constrained")), 3)) # % var explained by each axis.
screeplot(lcRDA) # Variation explained by each RDA axis.
round(vif.cca(lcRDA), 2) # Check model variance inflation factors.


# Define a function for extracting RDA elements and plotting the results.
# Will use this same function with the other data sources too.
rda.full <- \(model, dataset) {
  
  # Pull in the RDA model of interest.
  rdamodel <- get(x = as_label(eval(parse(text=enquo(model)))[[2]]),
                  envir = .GlobalEnv)
  
  # SNP loadings across first two RDA axes and the distributions of those loadings. 
  loadings <- scores(model, choices = c(1:2), display = 'species')
  png(paste0("../plots/", dataset, "_rda_loadings.png"), width = 1400, height = 600, res = 100)
  par(mfrow = c(1, 2)); hist(loadings[,1], main = "RDA1 Loadings", xlab = "RDA1")
  hist(loadings[,2], main = "RDA2 Loadings", xlab = "RDA2"); dev.off()
  print(paste0("Writing ", dataset, "_rda_loadings.png to ../plots/"), quote = F)
  
  # Identify and enumerate outlier loci along the first two RDA axes.
  cand1 <- outliers(loadings[,1], z = 3); print(paste0(length(cand1), " outliers on RDA1"), quote = F)
  cand2 <- outliers(loadings[,2], z = 3); print(paste0(length(cand2), " outliers on RDA2"), quote = F)
  
  # Isolate outlier loci into a pair of dataframes and adjust positional formatting.
  cand1DF <- cbind.data.frame(rep("RDA1", times=length(cand1)), names(cand1), unname(cand1))
  cand2DF <- cbind.data.frame(rep("RDA2", times=length(cand2)), names(cand2), unname(cand2))
  colnames(cand1DF) <- colnames(cand2DF) <- c("axis","snp","loading")
  all_cands <- rbind(cand1DF, cand2DF) %>% 
    separate(col = "snp", sep = 12,
             into = c("chromosome", "position"), remove = FALSE) %>% 
    mutate(chromosome = as.factor(gsub("_", "", chromosome)))
  
  # Only retain predictors (not conditioning variables).
  bio_vars_sub <- as.data.frame(bioclim_sub) %>% select(starts_with("bio"))
  
  # Make an empty matrix to populate in the following for loop.
  cand_mat <- matrix(nrow = nrow(all_cands), ncol = ncol(bio_vars_sub)) %>% 
    `colnames<-`(., c(colnames(bio_vars_sub)))
  
  # Estimate the correlation between the allele frequency of each outlier SNP and each bioclim variable.
  for (i in 1:nrow(cand_mat)) {
    nam <- all_cands[i, 2]; snp <- freqs[,nam]
    cand_mat[i,] <- apply(bio_vars_sub, 2, function(x) cor(x, snp)) } 
  
  # Add candidate SNP information to the correlation matrix computed above.
  cand_df_cor <- cbind.data.frame(all_cands, cand_mat) %>% rowwise() %>% 
    mutate(correlation = max(c_across(starts_with("bio")), na.rm = TRUE),
           predictor = names(.[6:(5+ncol(bio_vars_sub))])[which.max(c_across(6:(5+ncol(bio_vars_sub))))])  
  
  # Count how many outlier SNPs associate with each predictor variable.
  print(table(cand_df_cor$predictor))
  
  # Pivot correlation matrix to long form and merge SNP loadings along each RDA axis.
  dfwf <- cand_df_cor %>% 
    pivot_wider(names_from = "axis", values_from = "loading") %>% 
    merge(loadings, by.x = "snp", by.y = 0)
  
  # Join outlier SNPs with non-outlier SNPs for plotting purposes.
  RDAdf <- as.data.frame(loadings) %>% 
    merge(., dfwf[,c("snp", "predictor", "correlation")], 
          by.x = 0, by.y = "snp", all.x = TRUE) %>% 
    mutate(predictor = as.factor(case_when(
      !is.na(correlation) ~ predictor,
      is.na(correlation) ~ "non-outlier"))) %>% 
    separate(col = Row.names, c(12), into = c("chromosome", "position")) %>% 
    mutate(chromosome = as.factor(gsub("_", "", chromosome)),
           position = as.numeric(position))
  
  # Percent variation explained by each RDA axis.
  perc_var <- round(100*summary(rdamodel)$cont$importance[2, 1:2], 2)
  
  # Order site scores with respect to Latitude and join with RDA loadings.
  site_scores <- merge(sites, as.data.frame(scores(rdamodel, display = 'sites', scaling = 3)),
                       by.x = "Site", by.y = 0) %>% arrange(Latitude) %>% mutate(Latitude = as.factor(Latitude))
  
  # Script for scaling = 3 plot.
  (rda_scale3 <- ggplot() +
      geom_point(data = as.data.frame(scores(rdamodel, display = 'species', 
                                             choices = c(1,2), scaling = 3)),
                 aes(x = RDA1, y = RDA2), colour = 'grey70', shape = 21, fill = 'grey90') +
      geom_point(data = site_scores, aes(x = RDA1, y = RDA2, fill = factor(Latitude)), shape = 21, size = 2) +
      scale_fill_manual(values = c(viridis_pal(option = "D")(length(unique(site_scores$Site)))),
                        labels = levels(unique(site_scores$Site))) +
      geom_segment(data = scores(rdamodel, display = 'bp', choices = c(1,2), scaling = 3),
                   aes(xend = RDA1*25, yend = RDA2*25, x = 0, y = 0), colour = '#0868ac',
                   lineend = "round", arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(data = as.data.frame(scores(rdamodel, display = 'bp', choices = c(1,2), 
                                            scaling = 3)) %>% rownames_to_column("bio"), colour = '#0868ac',
                aes(x = 26.5*RDA1, y = 26.5*RDA2, label = bio)) +
      theme_bw() + guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
      theme(legend.position = "none", legend.title = element_blank()) +
      labs(x = paste0("RDA1 (", perc_var[1], "%)"),
           y = paste0("RDA2 (", perc_var[2], "%)")))
  ggsave(paste0("../plots/", dataset, "_rda_s3.tiff"), width = 8, height = 8)
  print(paste0("Writing ", dataset, "_rda_s3.tiff to ../plots/"), quote = F)
  
  # Script for scaling = 1 plot.
  (rda_scale1 <- ggplot() +
      geom_point(data = RDAdf[is.na(RDAdf$correlation),], 
                 aes(x = RDA1, y = RDA2), colour = 'grey70',
                 shape = 21, fill = "gray95") +
      geom_point(data = RDAdf[!is.na(RDAdf$correlation),],
                 aes(x = RDA1, y = RDA2, fill = predictor),  
                 shape = 21, size = 2, alpha = 3/5) +
      scale_fill_manual(values = c(viridis_pal(option = "H")(length(unique(RDAdf$predictor)))),
                        labels = levels(unique(RDAdf$predictor))) +
      labs(x = paste0("RDA1 (", perc_var[1], "%)"),
           y = paste0("RDA2 (", perc_var[2], "%)")) +
      geom_segment(data = scores(rdamodel, display = 'bp', choices = c(1,2), scaling = 1),
                   aes(xend = RDA1, yend = RDA2, x = 0, y = 0), colour = '#0868ac', size = 1,
                   lineend = "round", arrow = arrow(length = unit(0.1, "inches"))) +
      geom_text(data = as.data.frame(scores(rdamodel, display = 'bp', choices = c(1,2), 
                                            scaling = 1)) %>% rownames_to_column("bio"), colour = '#0868ac',
                aes(x = 1.06*RDA1, y = 1.06*RDA2, label = bio)) + theme_bw() +
      theme(legend.title = element_blank(), legend.position = c(0.08, 0.8), legend.background = element_blank())) 
  
  # Return base data frame and two ggplot objects as a list.
  return(list(RDAdf = RDAdf, rda_scale3 = rda_scale3, rda_scale1 = rda_scale1))
  
}

# Run above function on the RDA defined earlier.
# data argument is just for appending to output files/figures.
    # Use lcwgs or imputed_lcwgs, etc.
lcwgs_rda <- rda.full(model = lcRDA, data = "lcwgs")

# Plot both scaling = 1 and scaling = 3 plot types. 
(lcrda_plots <- cowplot::plot_grid(plotlist = list(
  lcwgs_rda$rda_scale3, lcwgs_rdaa$rda_scale1), 
  ncol = 2, labels = c("a)", "b)")))
ggsave("../plots/lcwgs_sox_rda.tiff", width = 16, height = 8)




