setwd("~/sockeye_lcwgs/genotypes/bam_pop_sep")

library(tidyverse)

# Writes a text file containing bam files per population per line.
# First chunk extracts sample identity.
files <- read.table("../../data/bam_list_n377.txt", col.names = c("file")) %>% 
  mutate(fish_ID = gsub(".dedup.clip.rg.bam", "", 
                        gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 

# Bring in sample information.
# Adds population names to individual samples.
# Requires some small amount of cleanup for current purposes.
samples <- read.csv("../../data/LandGen-Sockeye-Ext-Data.csv")[,c(1,6)] %>% 
  group_by(fish_ID) %>% sample_n(1) %>% 
  mutate(site_full = gsub("Shushwap", "Shuswap", site_full),
         site_full = gsub("\\(Lower)", "Lower", site_full))

# Converts to a list of dataframes.
# Each contains bam info per population.
# Select only file names at the end.
dat <- merge(files, samples, by = "fish_ID") %>% 
  mutate(siteName  = tolower(gsub(" ", "_", site_full)),
         direcfile = paste0("03_alignments_dedup_clip/", file)) %>% 
  group_by(site_full) %>% 
  split(., f = .$siteName) %>% 
  lapply(., function(x) x[,"direcfile"])

# Write each list of bam files per population as a separate and named txt file.
mapply(write.table, dat, file=paste0(tolower(names(dat)), "_bams.txt"),
       MoreArgs = list(quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\n"))

print(c(names(dat)), quote = F) # Check names.


# SCP to GPSC ------------------------------------------------------------------
# scp ../dir/*.txt ../gpsc/judsonb/sockeye_lcwGS/05_pop_bams/text_files
# dos2unix 05_pop_bams/*.txt