# PlotPosterior.R
# Written by Tal Zinger
# 
# This script plots a posterior distribution given by FITS, 
# for all given positions.
#
# This script is intentionally written with "hard-coded" output files,
# with a novice user in mind. Users are encouranged to modify this file
# to suit their specfic needs.
#

# in order to install the packages - run the commands:
# install.packages("ggplot2")
# install.packages("reshape2")

library(ggplot2)
library(reshape2)

args.vec <- commandArgs(trailingOnly = T)

file.name <- args.vec[1]

posterior.df <- read.table( file.name, header = T, sep = "\t" )
posterior.df$distance <- NULL
posterior.long.df <- reshape2::melt( id.var="position", data=posterior.df )

per.parameter.plot <- ggplot(posterior.long.df) + geom_histogram( aes(x=value, fill=variable), bins = 100 ) + facet_wrap(~position, scales = "free_x")
per.position.plot <- ggplot(posterior.long.df) + geom_histogram( aes(x=value, fill=position), bins = 100 ) + facet_wrap(~variable, scales = "free_x")

ggsave( per.parameter.plot, filename = "posterior_per_parameter.png" )
ggsave( per.position.plot, filename = "posterior_per_position.png" )
