# PlotSimulation.R
# Written by Tal Zinger
#
# Receives a simulation output from FITS ans plots the trajectories of
# all simulated alleles, as well as mutant alleles only.
#
# This script is intentionally written with "hard-coded" output files,
# with a novice user in mind. Users are encouranged to modify this file
# to suit their specfic needs.

library(ggplot2)

args.vec <- commandArgs(trailingOnly = T)

file.name <- args.vec[1]

simulation.df <- read.table( file.name, header = T, sep = "\t" )
simulation.df$allele <- as.factor( simulation.df$allele )

first.generation <- simulation.df$gen[1]
first.gen.df <- simulation.df[ simulation.df$gen==first.generation,]
first.gen.df <- first.gen.df[ order(-first.gen.df$freq), ]
common.allele <- first.gen.df$allele[1]

all.plot <- ggplot(simulation.df) + geom_line( aes(x=gen, y=freq, colour=allele) )
mutants.plot <- ggplot(simulation.df[ simulation.df$allele != common.allele, ]) + geom_line( aes(x=gen, y=freq, colour=factor(allele)) )

ggsave( all.plot, filename = "simulation_all_allales.png" )
ggsave( mutants.plot, filename = "simulation_mutant_alleles.png" )
