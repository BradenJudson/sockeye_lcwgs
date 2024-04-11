setwd("~/sockeye_lcwgs/genotypes/n_reads")

library(tidyverse)

reads <-  read.table("all_flagstat.txt",
                     sep = "\t")