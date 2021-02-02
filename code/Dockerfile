FROM rocker/geospatial:4.0.2
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN install2.r --skipinstalled --error cowplot \
 doParallel \
 ggthemes \
 ggrepel \
 spaero \
 viridis
RUN install2.r --repos "https://mran.microsoft.com/snapshot/2019-04-26" --error pomp
