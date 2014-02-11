# mp_tutorial_taxonomic_analysis.R
# This script summarizes all R-based analysis from the 
# Taxonomic Analysis presentation of the MetaPathways tutorial
# Tuesday, February 11, 2014
# Niels W. Hanson (nielsh@mail.ubc.ca)

## 1. Loading Megan Table into R
# set the working directory to where you have the HOT_mean_table.txt
# Download: https://raw.github.com/nielshanson/mp_tutorial/master/taxonomic_analysis/files/HOT_megan_table.txt
setwd("~/mp_tutorial/taxonomic_analysis/files/")
HOT_data <- read.table("HOT_megan_table.txt", header=TRUE, sep="\t", row.names=1)

## 2. Data Scaling
# alway important to plot your data to get a sense of its distribution
par(mfrow=c(1,2)) # multiple plots in one graph
hist(as.matrix(HOT_data))
hist(as.matrix(log(HOT_data + 1)))
HOT_data <- log(HOT_data + 1)
par(mfrow=c(1,2))

## 3. Distance Matrices and Clustering
HOT_data.t <- t(HOT_data) # transpose HOT_data
HOT_data.t.dist <- dist(HOT_data.t) # calculate the euclidian distance between the rows
HOT_data.t.dist
# Note that dist uses the Euclidian distance by default

# Next we will hierarhcally cluster our distance matrix using hclust
HOT_data.euclid.fit <- hclust(HOT_data.t.dist)
plot(HOT_data.euclid.fit, main="Distance: Euclidian, Clustering:Complete")

# a popular distance with biologists is the Bray-Curtis simmilarity
try(library("ecodist"), install.packages("ecodist")) # try load, otherwise install
library("ecodist") # try to load again
HOT_data.t.bcdist <- bcdist(HOT_data.t)
HOT_data.bcdist.ward.fit <- hclust(HOT_data.t.bcdist, method="ward")
plot(HOT_data.bcdist.ward.fit, main="Distance: Bray-Curtis, Clustering:Ward's")

# a good way to visualize your clustering and distance matrix is via a heatmap
# remember to use install.packages("package_name") if you can not load
library("gplots")
library("RColorBrewer") # for colours
euclid_dend <- as.dendrogram(HOT_data.euclid.fit) # get ordering for heatmap
my_colours <- brewer.pal(8,"GnBu") # my colours
heatmap.2(as.matrix(HOT_data.t.dist), margin=c(14,14), 
                                      Rowv=euclid_dend, 
                                      Colv=euclid_dend, 
                                      col=my_colours, 
                                      trace="none", 
                                      denscol="black")

# do the same for the Bray-Cutris clustering
bc_dend <- as.dendrogram(HOT_data.bcdist.ward.fit) # get ordering for heatmap
my_bc_colours <- brewer.pal(8,"BuPu") # my colours
heatmap.2(as.matrix(HOT_data.t.bcdist), margin=c(14,14), 
                                        Rowv=bc_dend, 
                                        Colv=bc_dend, 
                                        col=my_bc_colours, 
                                        trace="none", 
                                        denscol="black")

# finally the last think is to have some statistical confidence about our clustering
# the main way is via bootstrapped p-values in the pvclust package. Its a good idea to 
# resample about 1000 times

library(pvclust)
HOT_data.pv_fit <- pvclust(HOT_data, method.hclust="complete", method.dist="euclidian", n=1000) # in this case no transform is need
plot(HOT_data.pv_fit)

# unfortunately, pvcust does not have bcdist implemented but it was not that hard to add it
library("devtools") # used to source functions from the internet
source_url('http://raw.github.com/nielshanson/mp_tutorial/master/taxonomic_analysis/code/pvclust_bcdist.R')
# alternatively, you might want to download the above function and load using source()
# e.g. source('/where/I/downloaded/pvclust_bcdist.R')
HOT_data.bcdist.pv_fit <- pvclust(HOT_data, method.hclust="ward", method.dist="bray–curtis", n=1000)
plot(HOT_data.bcdist.pv_fit)

## 4. Visualizing taxonomic abundances

# one quick way is to revisit our heatmap.2 function

my_colours <- brewer.pal(8,"Blues") # another colour scheme
HOT_heat <- heatmap.2(as.matrix(HOT_data), margin=c(14,14), 
                                               col=my_colours, 
                                               Colv=bc_dend, 
                                               trace="none", 
                                               denscol="black")

# finally a rather impressive looking bubble plot can be made with 
# the visualizion package ggplot2
library(ggplot2)
library(reshape2) # transform our Data from wide to long format
HOT_data$taxa = rownames(HOT_data) # add taxa from rownames to Data Frame
HOT_data.m <- melt(HOT_data)
colnames(HOT_data.m)[2] = "sample" # rename variable column to sample

# in order to plot things in properly, the order of each variable has to be explicitly set
name_order <- HOT_data.bcdist.ward.fit$labels[HOT_data.bcdist.ward.fit$order] # get order of samples from clustering
HOT_data.m$sample <- factor(HOT_data.log.m$sample, levels=name_order) # set order of samples
HOT_data.m$taxa <- factor(HOT_data.log.m$taxa, levels=unique(HOT_data.log.m$taxa)) # set order of taxa

# cut bray-curtis clustering to get groups
bc_ward_groups <- cutree(HOT_data.bcdist.ward.fit, h=0.2) # slice dendrogram for groups (hight=0.2)
HOT_data.m$clust_group <- as.vector(bc_ward_groups[as.vector(HOT_data.m[,"sample"])])
HOT_data.m$clust_group <- as.factor(HOT_data.m$clust_group) # set group numbers as factors

# finally create the bubble plot
g <- ggplot(subset(HOT_data.m, value >0), aes(x=sample, y=taxa, color=clust_group))
g <- g + geom_point(aes(size=value)) # plot the points and scale them to value
g <- g + theme_bw() # use a white background
g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # rotate and centre labels
g # plot it, whew!
