## mp_tutorial_downstream_analysis.R
# Niels Hanson (nielsh@mail.ubc.ca)
# Summary of Introduction to Downstream Analysis in R talk
# Hydrocarbon MetaPathways Workshop, Thursday February 13 2014

## 1. Loading Data into R 

# set working directory to where your data is
# e.g. the downstream_analysis_r/code/ directory
setwd("~/mp_tutorial/downstream_analysis_r/code/")

# read in data
# the "wide", the "lookup" table, and the "metadata" table of 
# samples to experimental/environmental conditions
pathways_wide <- read.table("../files/HOT_Sanger_pwy.wide.txt", sep="\t", header=T, row.names=1)
hot_metadata <- read.table("../files/HOT_Sanger_ex_var.csv.txt", sep="\t", header=T)

# quick check to look at the head of our table to check naming conventions
head(pathways_wide)
tail(pathways_wide)
dim(pathways_wide)
head(hot_metadata)
tail(hot_metadata)
dim(hot_metadata)

## 2. Slicing and Dicing Data Frames
pathways_wide[1,] # first row
pathways_wide[,1] # first column

# select by row and column names instead of numbers
pathways_wide["SUCSYN-PWY",]
pathways_wide[,"X1_upper_euphotic"]

# the dollarsign operator can be used to access columns
pathways_wide$"X1_upper_euphotic"

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

# there is a nice convenience function if you want to use the header names
subset(pathways_wide, "X1_upper_euphotic" > 0)

## 1. E.g. Venn Diagrams
# load venn diagram functions
try(library("devtools"), install.packages("devtools")) # used to source functions from the internet
library("devtools")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram2.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")

# can use rownames and colnames to get pathway names and column variable names
rownames(pathways_wide)
colnames(pathways_wide)

# use our newly found skills to identify which pathways had a signal in each sample
pwys_10m <- rownames(pathways_wide)[pathways_wide[,"X1_upper_euphotic"] > 0]
pwys_70m <- rownames(pathways_wide)[pathways_wide[,"X6_upper_euphotic"] > 0]
pwys_130m <- rownames(pathways_wide)[pathways_wide[,"X2_chlorophyllmax"] > 0]
pwys_200m <- rownames(pathways_wide)[pathways_wide[,"X3_below_euphotic"] > 0]
pwys_500m <- rownames(pathways_wide)[pathways_wide[,"X5_uppermesopelagic"] > 0]
pwys_770m <- rownames(pathways_wide)[pathways_wide[,"X7_omz"] > 0]
pwys_4000m <- rownames(pathways_wide)[pathways_wide[,"X4_deepabyss"] > 0]

# We can now use these as inputs to the our venn_diagram scripts
quartz()
venn_10m_and_4000m <- venn_diagram2(pwys_10m, pwys_4000m,
                                    "10m", "4000m")
quartz()
venn_10m_70m_130m <- venn_diagram3(pwys_10m, pwys_70m, pwys_130m,
                                   "10m", "70m", "130m")
quartz()
venn_500m_770m_4000m <- venn_diagram3(pwys_500m, pwys_770m, pwys_4000m,
                                      "500m", "770m", "4000m")
quartz()
venn_10m_70m_130m_200m <- venn_diagram4(pwys_10m, pwys_70m, pwys_130m, pwys_200m,
                                        "10m", "70m", "130m", "200m")
quartz()
venn_200m_500m_770m_4000m <- venn_diagram4(pwys_200m, pwys_500m, pwys_770m, pwys_4000m,
                                           "200m", "500m", "770m", "4000m")

# it can also be valuable to compare interesting pathway sets from the above venn_diagrams
# e.g. pathways common to 10m, 70m, and 130m, against pathways common to
# 500m, 770m, and 4000m
my_colors = c("#29A7A7", "#C000C0") # custom colors
# combine function 'c' turns arguments into vectors

# the euler option attempts to scale the relative sizes of the 
# cirles (not always possible for more complex diagrams)
quartz()
compare_cores <- venn_diagram2(venn_500m_770m_4000m$"500m_770m_4000m", venn_10m_70m_130m$"10m_70m_130m",
                               "Deep", "Surface", colors=my_colors, euler=TRUE)

## 2. Data distributions

# relativeize counts by their column sums
pathways_wide.rel <- scale(pathways_wide, center=FALSE, scale=colSums(pathways_wide))

# take a look at the two distributions
par(mfrow=c(1,2))
hist(as.matrix(pathways_wide), xlab="Raw ORFs", main="Raw Counts", col="light blue")
hist(as.matrix(pathways_wide.rel), xlab="Relative ORFs", main="Relative Counts", col="light green")
par(mfrow=c(1,1)) # reset graph to 1x1
# both relative and raw distributions looks pretty similar

# also good to try a few transforms
par(mfrow=c(2,3))
hist(pathways_wide.rel, xlab="Relative ORFs", main="hist", col="dark blue")
hist(log(pathways_wide.rel+1), xlab="Relative ORFs (Log)", main="hist (log)", col="dark green")
hist(sqrt(pathways_wide.rel), xlab="Relative ORFs (Square Root)", main="hist (sqrt)",col="tomato")
plot(density(pathways_wide.rel), xlab="Relative ORFs", main="density", col="dark blue")
plot(density(log(pathways_wide.rel+1)), xlab="Relative ORFs (Log)", main="density (log)", col="dark green")
plot(density(sqrt(pathways_wide.rel)), xlab="Relative ORFs (Square Root)", main="density (sqrt)", col="tomato")
par(mfrow=c(1,1))

# r has a few special names for colors
colors()

# we might also take this opportunity to highlight the pairs function to
# compare envinronmental paramter metadata
quartz()
pairs(hot_metadata[2:length(hot_metadata)])

# we can do the same using qqplot using the histgram and density geom patterns
try (library("ggplot2"), install.packages("ggplot2"))
library("ggplot2")

# since ggplot is a different visualization environment than the standard R graphics
# a separate function as been created to display plots together
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/multiplot.r")

# in this case ggplot will expect our data as a vector of numbers
pathways_wide.vector <- as.vector(as.matrix(pathways_wide))
pathways_wide.rel.vector <- as.vector(pathways_wide.rel)

# same as above but in qplot
p1 <- qplot(pathways_wide.vector, 
            geom="histogram",
            xlab="Raw ORFs", main="Raw Counts")
p2 <- qplot(pathways_wide.rel.vector, 
            geom="histogram",
            xlab="Relative ORFs",
            main="Relative Counts")
quartz()
multiplot(p1, p2, cols=2)

# same as above but in qplot
p1 <- qplot(pathways_wide.rel.vector, xlab="Relative ORFs", main="hist", geom="histogram")
p2 <- qplot(log(pathways_wide.rel.vector + 1), xlab="Relative ORFs (Log)", main="hist (log)", geom="histogram")
p3 <- qplot(sqrt(pathways_wide.rel.vector), xlab="Relative ORFs (Square Root)", main="hist (sqrt)", geom="histogram")
p4 <- qplot(pathways_wide.rel.vector, xlab="Relative ORFs", main="density", geom="density")
p5 <- qplot(log(pathways_wide.rel.vector + 1), xlab="Relative ORFs (Log)", main="density (log)", geom="density")
p6 <- qplot(sqrt(pathways_wide.rel.vector), xlab="Relative ORFs (Square Root)", main="density (sqrt)", geom="density")
quartz()
multiplot(p1, p4, p2, p5, p3, p6, cols=3)

# though this is deprecated for GGally, its ability to quickly work with dataframes
# is still valuable
quartz()
plotmatrix(log(hot_metadata[2:length(hot_metadata)]), colour="gray20") + geom_smooth(method="lm")

## 3. Hierarchical Clustering
source_url('http://raw.github.com/nielshanson/mp_tutorial/master/taxonomic_analysis/code/pvclust_bcdist.R')
HOT_data.bcdist.pv_fit <- pvclust(HOT_data, method.hclust="ward", method.dist="bray–curtis", n=1000)
plot(HOT_data.bcdist.pv_fit)




# time to switch from wide to long table formats using melt
try( library("reshape2"), install.packages("reshape2") )
library("reshape2")

# add pathways from rowsnames to matrix
pathways_wide$pwy = rownames(pathways_wide)
# go from wide to long table format
pathways_long <- melt(pathways_wide)
colnames(pathways_long)[2] = "samp"
pathways_long$samp <- sub("X","", pathways_long$samp)

# add long pathway names 
meta_17 <- read.table("../files/meta_17.txt", sep="\t", header=T, row.names=1)
pathways_long$pwy_long = meta_17[pathways_long$pwy,1]

# download metacyc hierarchy & level the factors
meta_17_hier <- read.table("../files/meta_17_hierarchy.txt", sep="\t", header=F, row.names=1)
meta_17_hier$V3 <- factor(meta_17_hier$V3, levels=unique(meta_17_hier$V3))
meta_17_hier$V4 <- factor(meta_17_hier$V4, levels=unique(meta_17_hier$V4))
meta_17_hier$V5 <- factor(meta_17_hier$V5, levels=unique(meta_17_hier$V5))

# add additional pathway levels 
pwy_level1 <- meta_17_hier[pathways_long$pwy,2]
pwy_level2 <- meta_17_hier[pathways_long$pwy,3]
pathways_long <- cbind(pathways_long, pwy_level1)
pathways_long <- cbind(pathways_long, pwy_level2)

# small correction for pathways not found in hierarchy
missing <- setdiff(unique(pathways_long$pwy), rownames(meta_17_hier))
pathways_long <- pathways_long[!(pathways_long$pwy %in% missing),]

# order pathways first by level V3 and then by V4
pwy_order <- intersect(rownames(meta_17_hier[order(meta_17_hier[,"V3"], meta_17_hier[,"V4"]),]), 
                       unique(pathways_long$pwy))

# factor pathways and then the samples 
pathways_long$pwy <- factor(pathways_long$pwy, levels = pwy_order)
pathways_long$samp <- factor(pathways_long$samp, levels = hot_metadata$sample)

# whew, you got the longtable

# load ggplot2


## Visualization in ggplot2
# load some required packages, otherwise install them and try again

library(ggplot2)

g <- ggplot(subset(pathways_long, value >0), aes(x=samp,y=pwy_level1)) +
     geom_point(aes(size=sqrt(value), color=samp)) + 
     theme_bw() + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways")
g

g <- ggplot(subset(pathways_long, value >0), aes(x=samp,y=pwy_level2)) +
     geom_point(aes(size=sqrt(value), color=samp)) + 
     theme_bw() + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     labs(x = "Samples", y = "Pathways")
g

g <- ggplot(pathways_long, aes(pwy_level1)) +
     geom_bar(aes(fill=samp)) + 
     theme_bw() + 
     coord_flip() + 
     facet_wrap(~ samp, ncol = 7)
g

g <- ggplot(pathways_long, aes(pwy_level2)) +
     geom_bar(aes(fill=samp)) + 
     theme_bw() + 
     coord_flip() + 
     facet_wrap(~ samp, ncol = 7)
g

g <- ggplot(pathways_long, aes(pwy_level2, fill=samp)) +
     geom_bar() + 
     theme_bw() + 
     coord_flip()
g

g <- ggplot(pathways_long, aes(value)) +
     geom_density( aes(fill=samp, alpha=0.6)) + 
     facet_wrap(~ samp, ncol=7) +
     theme_bw() + 
     theme(legend.position="none")
g

g <- ggplot(pathways_long, aes(log(value+1))) +
     geom_density( aes(fill=samp, alpha=0.6)) + 
     facet_wrap(~ samp, ncol=7) +
     theme_bw() + 
     theme(legend.position="none")
g

g <- ggplot(pathways_long, aes(value)) +
     geom_histogram( aes(fill=samp)) + 
     facet_wrap(~ samp, ncol=7) +
     theme_bw() + 
     theme(legend.position="none")
g  

g <- ggplot(pathways_long, aes(log(value + 1))) +
     geom_histogram( aes(fill=samp)) + 
     facet_wrap(~ samp, ncol=7) +
     theme_bw() + 
     theme(legend.position="none")
g
