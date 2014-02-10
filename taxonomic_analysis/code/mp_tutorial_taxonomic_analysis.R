# mp_tutorial_taxonomic_analysis.R
# This script summarizes all R-based analysis from the 
# Taxonomic Analysis presentation of the MetaPathways tutorial
# Tuesday, February 11, 2014
# Niels W. Hanson (nielsh@mail.ubc.ca)

# 1. Loading data
# set the working directory to where you have the HOT_mean_table.txt
setwd("~/Desktop/HOT")
HOT_data <- read.table("HOT_megan_table.txt", header=TRUE, sep="\t", row.names=1)



HOT_data.t <- t(HOT_data) # transpose HOT_data
HOT_data.t.dist <- dist(HOT_data.t) # calculate the euclidian distance between the rows
HOT_data.t.dist
