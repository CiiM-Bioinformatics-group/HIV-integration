suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
  library(reticulate)
  library(clusterProfiler)
  INT <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
})

#Which views are more affected by the different institute effect
