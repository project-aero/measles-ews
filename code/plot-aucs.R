# plot-aucs.R
#  Script to plot the AUC results from emergence and endemic simulation
#  tests of EWS.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)



# Load libraries ----------------------------------------------------------

library(tidyverse)


# Load results ------------------------------------------------------------

endemic_aucs <- read.csv("../results/endemic-aucs.csv") %>%
  mutate(simulation = "Endemic")

emergence_aucs <- read.csv("../results/emergence-aucs.csv") %>%
  mutate(simulation = "Emergence")

auc_tbl <- bind_rows(
  endemic_aucs, 
  emergence_aucs
)


# Plot and save -----------------------------------------------------------

auc_plot <- ggplot(auc_tbl, aes(x = metric, y = abs(AUC-0.5), fill = AUC)) +
  geom_col(position = position_dodge()) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_fill_viridis_c(limits = c(0,1), direction = -1, option = "C") +
  facet_grid(simulation~city) +
  theme_minimal() +
  labs(x = NULL, y = "| AUC - 0.5 |")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.spacing = unit(1, "lines"))

ggsave(filename = "../figures/simulation-aucs.pdf", plot = auc_plot, width = 8.5, height = 5, units = "in")
