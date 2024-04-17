setwd("~/sockeye_lcwgs/genotypes/cov_reads")

library(tidyverse); library(cowplot); library(ggpmisc); library(ggExtra)

# Read n indvidual coverage data.
cov <- read.table("cov.txt", sep = "",header = F, 
                  col.names = c("file", "av", "eq", "coverage"))[,c(1,4)] 

# Read in individual read count data.
reads <- read.table("reads.txt", sep = "", header = F,
                    col.names = c("file", "reads")) 

# Combine everything into one object.
dat <- merge(cov, reads, by = "file") %>% 
  mutate(sample = gsub("^.*\\.","", str_sub(start = 0, end = nchar(file)-18, file)),
         abb    = as.factor(str_sub(start = 0, end = 3, sample))) %>% 
  relocate("sample", .before = "file")

# Plot relationship between reads and coverage. 
(lre <- ggplot(data = dat,
        aes(x = reads/1e6,
            y = coverage)) + 
    scale_x_continuous(breaks = seq(0, 90, 1)) +
    scale_y_continuous(breaks = seq(0,  4, 1/2)) +
    coord_cartesian(clip = 'on') +
    geom_segment(aes(y = mean(`coverage`), 
                     yend = mean(`coverage`), 
                     x = -Inf, xend = Inf),
                 linewidth = 1/2, linetype = 2) +
    geom_segment(aes(x = mean(`reads`)/1e6,
                     xend = mean(`reads`)/1e6,
                     y = -Inf, yend = Inf),
                 linewidth = 1/2, linetype = 2) +
    geom_smooth(method = "lm",
                colour = "red",
                alpha  = 1/6,
                linetype = 1) +
    geom_point(colour = "black",
               shape  = 21, 
               fill   = "gray80",
               size   = 2) +
    theme_bw() +
    theme(panel.grid = element_line(color = "gray95")) +
    labs(x = "Reads (millions)", 
         y = "Individual Coverage") +
    stat_poly_eq(use_label(c("R2", "p")),
                 label.x = "left",
                 label.y = "top",
                 small.p = TRUE))

# Add density plots of marginal distributions.
(mds <- ggMarginal(lre, type = "density", colour = "black", 
                   linewidth = 1, fill = "gray90", alpha = 0.7))

save_plot("../../plots/lcwgs_reads_coverage.tiff", mds, ncol = 2, 
          base_height = 6, base_asp = 1)

# Tabularized summary stats.
(sstat <- dat %>% 
    pivot_longer(cols = c("reads", "coverage")) %>% 
    group_by(name) %>% 
    summarise(mean_value = mean(value),
              med_value  = median(value),
              max_val    = max(value),
              min_val    = min(value),
              sd         = sd(value)))
