# plot-variance-ews.R
#  This script plots the distributions of variance across the 500 simulations
#  in the null and test intervals. We plot the variance for all cities and for
#  both re-emergence and elimination. The variance performs well in our study,
#  and is widely used in other studies, making it a good example.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)


# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(dplyr)


# Load EWS statistics -----------------------------------------------------

do_metric <- "variance"

emerge_ews <- read.csv("../results/ews-emergence.csv") %>%
  filter(metric == do_metric & susc_discount == 1e-04)  %>%
  mutate(
    metric = as.character(metric),
    metric = ifelse(metric == "variance", "Variance", metric),
    type = "Emerge."
  )

elimin_ews <- read.csv("../results/ews-elimination.csv") %>%
  filter(metric == do_metric & vacc_speed == 1.5e-05)  %>%
  mutate(
    metric = as.character(metric),
    metric = ifelse(metric == "variance", "Variance", metric),
    type = "Elimin."
  )

all_ews <- bind_rows(emerge_ews, elimin_ews)


# Plot the EWS distributions ----------------------------------------------

mycols <- c("#91CDF0", "#EF6677")

varplot <- ggplot(all_ews, aes(x = type, y = log(value), fill = half)) +
  geom_boxplot(outlier.size = 0.6, outlier.shape = 1, width = 0.5, lwd = 0.2) +
  facet_wrap(~city, scales = "free", nrow = 2) +
  scale_fill_manual(values = mycols, name = "Interval", labels = c("Null", "Test")) +
  labs(x = NULL, y = "log(Variance)") +
  theme_classic(base_size = 10) +
  theme(panel.spacing = unit(1, "lines"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(filename = "../figures/var-ews-example.pdf", plot = varplot, 
       width = 4, height = 4, units = "in")

# ggplot(filter(all_ews, type == "emergence"), aes(x = value)) +
#   geom_density(aes(fill = half), color = "white", alpha = 0.7) +
#   facet_wrap(~city, scales = "free") +
#   scale_fill_brewer(palette = "Set1", name = "Interval", labels = c("Null", "Test")) +
#   labs(x = "Variance", y = "Empirical density") +
#   theme_classic() +
#   theme(panel.spacing = unit(1, "lines"), strip.background = element_blank(),
#         strip.text = element_text(face = "bold"),
#         plot.title = element_text(face = "bold"),
#         legend.position = c(0.9, 0.3)) +
#   ggtitle("A. Variance EWS during emergence")
