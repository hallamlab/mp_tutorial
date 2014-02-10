# mp_tutorial_taxonomic_analysis.R
# This script summarizes all R-based analysis from the 
# Taxonomic Analysis presentation of the MetaPathways tutorial
# Tuesday, February 11, 2014
# Niels W. Hanson (nielsh@mail.ubc.ca)

## 1. Loading Megan Table into R
# set the working directory to where you have the HOT_mean_table.txt

setwd("~/Desktop/HOT")
HOT_data <- read.table("HOT_megan_table.txt", header=TRUE, sep="\t", row.names=1)
# if you are having a hard time you can just download it from the host
HOT_data <- read.table("HOT_megan_table.txt", header=TRUE, sep="\t", row.names=1)

## 2. Distance Matrices
HOT_data.t <- t(HOT_data) # transpose HOT_data
HOT_data.t.dist <- dist(HOT_data.t) # calculate the euclidian distance between the rows
HOT_data.t.dist
# Note that dist uses the Euclidian distance by default

# Next we will hierarhcally cluster our distance matrix using hclust
HOT_data.euclid.fit <- hclust(HOT_data.t.dist)
plot(HOT_data.euclid.fit, main="Distance: Euclidian, Clustering:Complete")

# a popular distance with biologists is the 