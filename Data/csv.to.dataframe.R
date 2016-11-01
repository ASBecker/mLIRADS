# Data Readin for mLI-RADS Paper
# Author: Anton Becker (@ASBecker)

library("gdata")
library("plyr")

d.re <- read.csv("lirads.1-5.csv")
d.re[d.re == ""] <- NA
d.re <- d.re[rowSums(is.na(d.re)) != ncol(d.re), ]
