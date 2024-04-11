setwd("~/sockeye_lcwgs/genotypes/n_reads")

library(tidyverse); library(ggExtra); library(cowplot)

# Read in flagstat data and wrangle into a dataframe. 
# Only formatting relevant columns for the time being.
reads <-  as.data.frame(matrix(read.table("all_flagstat.txt", 
                               sep = "\t", col.names = "all")[,1],
                               ncol = 14, byrow = T))  %>% 
  `colnames<-`(., c("file", "reads", "secondary", "supp", "duplicates",
                    "mapped", "pairedseq", "read1" ,"read2", "p.pair",
                    "self.paired", "singleton", "diffchr", "diff5")) %>% 
  mutate(reads = as.numeric(gsub("\\D+", "", reads))/10,
         file  = gsub("^.*\\.","", str_sub(start = 0, 
                      end = nchar(file)-18, file)),
         abb   = as.factor(str_sub(start = 0, end = 3, file))) %>% 
  filter(file != "merged")
         
     
# Plot reads per samle per population.
(rbp <- ggplot(data = reads, aes(x = abb, y = reads/1e6)) +
  geom_hline(yintercept = mean(reads$reads/1e6),
             linetype = 2, colour = "blue") +
  geom_hline(yintercept = mean(reads$reads/1e6),
             linewidth = 1.5, colour = "blue",
             alpha = 1/15) +
  geom_boxplot(alpha  = 1/10) +
  geom_point(size  = 3/2, 
             shape = 21,
             color = "black",
             fill  = "grey90") +
  theme_bw() + labs(x = NULL, y = "Reads (millions)") +
  theme(axis.text.x = element_text(angle = 90,
                      vjust = 0.5, hjust = 1),
        axis.text   = element_text(size  = 14),
        axis.title  = element_text(size = 16)))

# Add marginal histogram of read counts.
(mH <- ggMarginal(rbp, type = "histogram", margins = "y",
                  binwidth = 1, fill = "gray90", size = 15))

# Use cowplot here to save multi-plot configuration more easily.
save_plot("../../plots/sockeye_indv_reads.tiff", mH, ncol = 2, base_height = 7)

