## HOT_pathway_analysis.r

# set working directory to the HOT_pathway_analysis directory
setwd("~/mp_tutorial/downstream_analysis_r/code/")

# read in data
# the "wide", the "lookup" table, and the "metadata" table of 
# samples to experimental/environmental conditions
pathways_wide <- read.table("../files/HOT_Sanger_pwy.wide.txt", sep="\t", header=T, row.names=1)
hot_metadata <- read.table("../files/HOT_Sanger_ex_var.csv.txt", sep="\t", header=T)

# quick check to look at the head of our table to check naming conventions
head(pathways_wide)
head(hot_metadata)

# slicing and dicing
pathways_wide[1,] # first row
pathways_wide[,1] # first column

# To find which pathways had a more than one ORF we can use a logical operator
# to give a logical vector
pathways_wide[,1] > 1

# and these can be quite complex:
# "both column 1 and column 2 are true but column 3 is less than 10
complex_cond <- (((pathways_wide[,1] > 1) & (pathways_wide[,2] > 1) ) | pathways_wide[,3] > 0)

# Also to find out how many there are you can just sum the vector
sum( complex_cond )

# these are useful because we can now use these to select partiular parts of the matrix
# e.g. this only displays the pathways where column 1 was greater than 1
pathways_wide[pathways_wide[,1] > 1,]
# there is a nice convenience script if you want to use the header names
subset(pathways_wide, "X1_upper_euphotic" > 0)

## 1. E.g. Venn Diagrams
# load venn diagram functions
try(library("devtools"), install.packages("devtools")) # used to source functions from the internet
library("devtools")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram2.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")

rownames(pathways_wide)
colnames(pathways_wide)

# use our newly found skills to identify which pathways had a signal in each 
# sample
pwys_10m <- rownames(pathways_wide)[pathways_wide[,"X1_upper_euphotic"] > 0]
pwys_70m <- rownames(pathways_wide)[pathways_wide[,"X6_upper_euphotic"] > 0]
pwys_130m <- rownames(pathways_wide)[pathways_wide[,"X2_chlorophyllmax"] > 0]
pwys_200m <- rownames(pathways_wide)[pathways_wide[,"X3_below_euphotic"] > 0]
pwys_500m <- rownames(pathways_wide)[pathways_wide[,"X5_uppermesopelagic"] > 0]
pwys_770m <- rownames(pathways_wide)[pathways_wide[,"X7_omz"] > 0]
pwys_4000m <- rownames(pathways_wide)[pathways_wide[,"X4_deepabyss"] > 0]

# We can now use these as inputs to the our venn_diagram scripts
# 
venn_10m_and_4000m <- venn_diagram2(pwys_10m, pwys_4000m,
                                    "10m", "4000m")
venn_10m_70m_130m <- venn_diagram3(pwys_10m, pwys_70m, pwys_130m,
                                   "10m", "70m", "130m")
venn_500m_770m_4000m <- venn_diagram3(pwys_500m, pwys_770m, pwys_4000m,
                                      "500m", "770m", "4000m")
venn_10m_70m_130m_200m <- venn_diagram4(pwys_10m, pwys_70m, pwys_130m, pwys_200m,
                                        "10m", "70m", "130m", "200m")
venn_200m_500m_770m_4000m <- venn_diagram4(pwys_200m, pwys_500m, pwys_770m, pwys_4000m,
                                           "200m", "500m", "770m", "4000m")

# it can also be valuable to compare interesting pathway sets from the above venn_diagrams
# e.g. pathways common to 10m, 70m, and 130m, against pathways common to
# 500m, 770m, and 4000m
my_colors = c("#96B4D4", "#ED808B") # custom colors

compare_cores <- venn_diagram2(venn_10m_70m_130m$"10m_70m_130m", venn_500m_770m_4000m$"500m_770m_4000m",
              "Surface_Core", "Deep_Core", colors=my_colors)

# the euler option attempts to scale the relative sizes of the 
# cirles (not always possible for more complex diagrams)
venn_diagram2(venn_10m_70m_130m$"10m_70m_130m", venn_500m_770m_4000m$"500m_770m_4000m",
            "Surface_Core", "Deep_Core", colors=my_colors, euler=TRUE)

# Finally the out and file_name options allow you to write out a pdf file to the 
# current working directory
venn_diagram2(venn_10m_70m_130m$"10m_70m_130m", venn_500m_770m_4000m$"500m_770m_4000m",
            "Surface_Core", "Deep_Core", colors=my_colors, out = TRUE, file_name = "my_filename", euler=TRUE)

## 2. Hierarchcial Clustering, Heatmaps, Bubble Plots
# load some required packages, otherwise install them and try again
try( library("ecodist"), install.packages("ecodist") ) # Bray-Curtus dissimilarity
library("ecodist")
# Bootstrapped clustering
source_url('http://raw.github.com/nielshanson/mp_tutorial/master/taxonomic_analysis/code/pvclust_bcdist.R')

# Fit pvclust with bray-curtis distance and 10 bootstraps (in the interest of time)
hist(as.matrix(pathways_wide))
hist(as.matrix(log(pathways_wide+1)))
pathways_wide_sqrt <- sqrt(pathways_wide)
path_bcdist_sqrt.pvfit = pvclust(pathways_wide_sqrt, nboot=10, method.hclust="ward", method.dist="bray–curtis")
plot(path_bcdist_sqrt.pvfit)

# we note that our clustering via pathways is not the same as when we did it with taxonomy

# visualization libraries ggplot2
try( library("ggplot2"), install.packages("ggplot2") ) # Bubble plots
library("ggplot2")
try( library("RColorBrewer"), install.packages("RColorBrewer")) # Color palletes
library("RColorBrewer")


# now that we are looking at the pathways directory it will help to have the translated pathway names
# the MetaCyc hierarchy
meta_17_hier <- read.table("ref/meta_17_hierarchy.txt", sep="\t", header=F, row.names=1)
meta_17 <- read.table("ref/meta_17.txt", sep="\t", header=F, row.names=1)
meta_17_hier$V3 <- factor(meta_17_hier$V3, levels=unique(meta_17_hier$V3))
meta_17_hier$V4 <- factor(meta_17_hier$V4, levels=unique(meta_17_hier$V4))
meta_17_hier$V5 <- factor(meta_17_hier$V5, levels=unique(meta_17_hier$V5))

# another way to visualize a clustering is to plot a heatmap or bubble plot
# there are many options and below is (in my opinion) a good set:
# - sqrt transform: pathways tend to have lots of large and small values
# - cluster on samples with bcdist (columns)
# - Color Brewer Palette (http://colorbrewer2.org)
# - remove the trance line
# - outline of density black
# - increase the margins to fit better

# rownames(pathways_wide) <- as.vector(meta_17[order,1])
# quartz()
# heatmap.2(sqrt(as.matrix(pathways_wide[order,])), 
#           distfun = bcdist,
#           dendrogram = "col",
#           col = brewer.pal(8, "Blues"),
#           trace = "none",
#           density.info = "histogram",
#           denscol="black",
#           margins = c(10,10),
#           xlab="Samples",
#           ylab="Pathways")
 
# to see all the colorbrewer paletes
# quartz()
# display.brewer.all()
 
# we might bring in ggplot at this point as a visualizaion framework
# bubble plots which use bubble size rather than color have a better dynamic range 
# than heatmaps

# load pathways lookup table
pathways_lookup <- read.table("data/HOT_pwys_lookup.txt", sep="\t", header=T)

# small correction for pathways not found in hierarchy
missing <- setdiff(intersect(pathways_lookup$PWY_NAME, pathways_lookup$PWY_NAME), rownames(meta_17_hier))
pathways_lookup <- pathways_lookup[!(pathways_lookup$PWY_NAME %in% missing),]

pwy_order <- intersect(rownames(meta_17_hier[order(meta_17_hier[,"V3"], meta_17_hier[,"V4"]),]), unique(pathways_lookup$PWY_NAME))
pathways_lookup$PWY_NAME <- factor(pathways_lookup$PWY_NAME, levels = pwy_order)
pathways_lookup$SAMPLE <- factor(pathways_lookup$SAMPLE, levels = hot_metadata$sample)

pwy_level1 <- meta_17_hier[pathways_lookup$PWY_NAME,2]
pwy_level2 <-meta_17_hier[pathways_lookup$PWY_NAME,3]

pathways_lookup <- cbind(pathways_lookup, pwy_level1)
pathways_lookup <- cbind(pathways_lookup, pwy_level2)
colnames(pathways_lookup)

# rearrange according to clustering
pathways_lookup$SAMPLE <- factor(pathways_lookup$SAMPLE, levels = path_bcdist.pvfit$hclust$labels[path_bcdist.pvfit$hclust$order])
# cut clustering dedrogram to get groups
cluster_groups <- cutree(path_bcdist.pvfit$hclust, h=0.35)
pathways_lookup$cluster_group = cluster_groups[as.vector(pathways_lookup[,"SAMPLE"])]
pathways_lookup$cluster_group = as.factor(pathways_lookup$cluster_group)

# bubble plots
# decline and sum according to metacyc hierarchy
quartz()
g <- ggplot(subset(pathways_lookup, ORF_COUNT >0), aes(x=SAMPLE,y=pwy_level1)) + 
     geom_point(aes(size=sqrt(ORF_COUNT), color=cluster_group)) + theme_bw() +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways") 
g

quartz()
g <- ggplot(subset(pathways_lookup, ORF_COUNT >0), aes(x=SAMPLE,y=pwy_level2)) + 
     geom_point(aes(size=sqrt(ORF_COUNT), color=factor(cluster_group))) + theme_bw() +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways")
g

brew_colors = brewer.pal(8,"Blues")

# heatmap in ggplot
quartz()
p <- ggplot(subset(pathways_lookup, ORF_COUNT >0), aes(SAMPLE, pwy_level2)) + 
     geom_tile(aes(fill = sqrt(ORF_COUNT)), colour = "white") +
     scale_fill_gradient(low = brew_colors[1], high=brew_colors[8]) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways")
p

# rearrange according to depth
pathways_lookup$SAMPLE <- factor(pathways_lookup$SAMPLE, levels = hot_metadata$sample)

# so many pathways might be interesting to see what is common between 10m, 70m, 130m
pwys_10m_only <- subset(pathways_lookup, (PWY_NAME %in% compare_cores$Surface_Core_only))
quartz()
g <- ggplot(subset(pwys_10m_only, ORF_COUNT >0), aes(x=SAMPLE,y=PWY_COMMON_NAME)) + 
     geom_point(aes(size=sqrt(ORF_COUNT))) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways") 
g

# so many pathways might be interesting to see what is common between 10m, 70m, 130m
pwys_10m_only <- subset(pathways_lookup, (PWY_NAME %in% compare_cores$Deep_Core_only))
quartz()
g <- ggplot(subset(pwys_10m_only, ORF_COUNT >0), aes(x=SAMPLE,y=PWY_COMMON_NAME)) + 
     geom_point(aes(size=sqrt(ORF_COUNT))) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways") 
g

## 3. Dimentionality Reduction: PCA, NMDS
library(ade4)
library(vegan)
library(gclus)
library(ape)

# Princple Component Analysis (PCA)
# common to do a hellinger transform before PCA
pathways_wide.h <- t(decostand(t(pathways_wide), "hellinger"))
pathways_wide.h.pca <- rda(t(pathways_wide.h))
p <- length(pathways_wide.h.pca$CA$eig)
pathways_wide.h.pca.sc1 <- scores(pathways_wide.h.pca, display="wa", scaling=1, choices=c(1:p))
quartz("Pathways: PCA")
qplot(pathways_wide.h.pca.sc1[,1], pathways_wide.h.pca.sc1[,2], label=rownames(pathways_wide.h.pca.sc1), size=2, geom=c("point"), xlab="PC1 (25.1% Variance)", ylab="PC2 (19.5% Variance)", color=factor(cluster_groups)) + 
geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() + theme(legend.position="none") + xlim(-0.6, 0.6)

# alternatively one can use the regular plots command
quartz()
par(cex=1.0, cex.lab=1.0, cex.axis=1.0, las=1, bty="n", pch=16)
plot(pathways_wide.h.pca.sc1[,1], pathways_wide.h.pca.sc1[,2], col=cluster_groups, xlim=c(-0.7,0.4), ylim=c(-0.3,0.4), xlab="PC1", ylab="PC2")
# 


# Non-metric Multi-dimentional Scaling (NMDS)

pathways_wide.nmds <- metaMDS(t(pathways_wide), distance = "bray")

quartz("Pathways NMDS - Bray")
qplot(pathways_wide.nmds$points[,1], pathways_wide.nmds$points[,2], label=rownames(pathways_wide.nmds$points), size=2, geom=c("point"), 
      xlab="MDS1", ylab="MDS2", main=paste("NMDS/Bray - Stress =", round(pathways_wide.nmds$stress,3)), color=factor(cluster_groups)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() +theme(legend.position="none") + xlim(-0.5,1.0)

# alternatively one can always use the regular plots command
quartz()
par(cex=1.0, cex.lab=1.0, cex.axis=1.0, las=1, bty="n", pch=16)
plot(pathways_wide.nmds$points[,1], pathways_wide.nmds$points[,2], col=cluster_groups, xlab="MDS1", ylab="MDS2", xlim=c(-0.4,0.8), ylim=c(-0.0002,0.0002))