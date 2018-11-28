# plot-aucs.R
#  Script to plot the AUC results from emergence and endemic simulation
#  tests of EWS.
#
# Author:
#  Andrew Tredennick (atredenn@gmail.com)



# Load libraries ----------------------------------------------------------

library(tidyverse)
library(viridis)


# Load results ------------------------------------------------------------

emergence_aucs <- read.csv("../results/emergence-grid-aucs.csv")
elimination_aucs <- read.csv("../results/elimination-grid-aucs.csv")

star_tbl <- tibble(
  city = "Niamey",
  x = as.factor(rep(1e-04, 2)),
  y = c(5.7, 7.7)
)


# Make the plots ----------------------------------------------------------

emerge_plot <- ggplot() +
  geom_tile(data = emergence_aucs, aes(x = as.factor(susc_discount), y = metric, fill = abs(AUC-0.5))) +
  geom_text(data = star_tbl, aes(x = x, y = y, label = "*"), color = "white", size = 6) +
  scale_fill_viridis(limits = c(0,0.5), direction = -1, option = "C", name = "| AUC - 0.5 |") +
  facet_wrap(~city, nrow = 1) +
  labs(x = "Level of susceptible depletion", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.spacing = unit(1, "lines")) +
  ggtitle("Anticipating emergence")

eliminate_plot <- ggplot() +
  geom_tile(data = elimination_aucs, aes(x = as.factor(vacc_speed*10000), y = metric, fill = abs(AUC-0.5))) +
  scale_fill_viridis(limits = c(0,0.5), direction = -1, option = "C", name = "| AUC - 0.5 |") +
  facet_wrap(~city, nrow = 1) +
  labs(x = expression(paste("Rate to full vaccine coverage (", phantom()%*%phantom(), 10^4, ")")), y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.spacing = unit(1, "lines")) +
  ggtitle("Anticipating elimination")

auc_plot <- cowplot::plot_grid(emerge_plot, eliminate_plot, nrow = 2, align = "v", labels = "AUTO")

ggsave(filename = "../figures/simulation-grid-aucs.pdf", plot = auc_plot, width = 8.5, height = 5.5, units = "in")




# endemic_aucs <- read.csv("../results/endemic-aucs.csv") %>%
#   mutate(simulation = "Endemic")
# 
# emergence_aucs <- read.csv("../results/emergence-aucs.csv") %>%
#   mutate(simulation = "Emergence")
# 
# auc_tbl <- bind_rows(
#   endemic_aucs, 
#   emergence_aucs
# )
# 
# 
# # Plot and save -----------------------------------------------------------
# 
# auc_plot <- ggplot(auc_tbl, aes(x = metric, y = abs(AUC-0.5), fill = AUC)) +
#   geom_col(position = position_dodge()) +
#   scale_y_continuous(limits = c(0,0.5)) +
#   scale_fill_viridis_c(limits = c(0,1), direction = -1, option = "C") +
#   facet_grid(simulation~city) +
#   theme_minimal() +
#   labs(x = NULL, y = "| AUC - 0.5 |")+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(panel.spacing = unit(1, "lines"))
# 
# ggsave(filename = "../figures/simulation-aucs.pdf", plot = auc_plot, width = 8.5, height = 5, units = "in")
