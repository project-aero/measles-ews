# plot-global-search.R
#  Script to plot the results from the initial global parameter search using
#  particle filter. The 500 parameter sets with the highest likelihood are
#  saved for subsequent analysis.
#
# Author:
#  Andrew Tredennick


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(GGally)
library(ggthemes)


# Load data ---------------------------------------------------------------

global_lls <- readRDS("initial-search-lls.RDS") %>%
  as_tibble() %>%
  arrange(-loglik) %>%
  mutate(
    highest = c(
      rep(TRUE, 500),
      rep(FALSE, n()-500)
    )
  )


# Plot all pairs ----------------------------------------------------------

p <- ggpairs(
  global_lls, 
  columns = c(1,3:5,12:16),
  mapping = ggplot2::aes(color = highest),
  upper = "blank",
  diag = "blank",
  lower = list(continuous = wrap("points", alpha = 0.5, size=0.5))
) 

for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_color_manual(values = ptol_pal()(2)) 
  }
}

p

