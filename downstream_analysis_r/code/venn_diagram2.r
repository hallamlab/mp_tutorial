venn_diagram2 <- function(a, b, name_a = "A", name_b = "B", colors =  c("#e41a1c","#377eb8"), out = FALSE, file_name = "default.pdf", euler = FALSE){
	# Niels Hanson
	# Info: A function to create a venn digagram between three (3) overapping sets a, b, and c.
	# Requies: The VennDiagram R package, the program will try to download and install if you dont
	#          have it.
	# Inputs: a,b,c - vertical vectors or the first column of matrixes as sets of strings to be compared
	#         name_a, name_b, and name_c - names that you give to sets a, b, and c, respectively
	#         these will appear in the final image
	#         colors - a vector of colors defined in hex format ie. colors <- c("#B80830","#EACC33","#46E2D9")
	#         it might be a good idea that these mix well as the interlapping classes will be a mix of colors
	#         Colors must be specified in hex ie. "#123456", google "hex colors" for more info
	#         name_output - the name of the output image that will be put in R's current working directory
	#         you can set this yourself before running. ie. setwd("~/my-director/")
	
	# A: find all classes for set a, a_b, a_c, a_b_c
	# number in a, we will use this at the end for a sanity check
	
	# make sure that a,b,c are all sets, get rid of duplicates
	a <- intersect(a,a)
	b <- intersect(b,b)
	# length of a
	a_tot <- length(a)
	# length of b
	b_tot <- length(b)
	# number shared between a and b and c
	a_b <- intersect(a,b)
	a_b_tot <- length(a_b)
	a_b_tot
	a_only <- setdiff(a,b)
	a_only_tot <- length(a_only)
	b_only <- setdiff(b,a)
	b_only_tot <- length(b_only)

	# try to load VennDiagram else try to download it
	# More info: see Hanbo Chen and Paul C Boutros. BMC Bioinformatics 2011, 12:35 doi:10.1186/1471-2105-12-35
	try(library(VennDiagram), install.packages("VennDiagram")) 
	library(VennDiagram)
	

	# create the venn diagram and output into the current working directory
	temp <- list(
		name_a = c(1:length(a)),
		name_b = c((length(a_only)+1):((length(a_only))+length(b)))
		)
	names(temp) <- c(name_a,name_b)
	
	output <- venn.diagram(
		x = temp,
		filename = NULL,
		col = "black",
		fill = c(colors[1], colors[2]),
		alpha = 0.5,
		label.col = c("black"),
		cex = 2.5,
		fontfamily = "serif",
		fontface = "bold",
		cat.default.pos = "text",
		cat.col = c("black"),
		cat.cex = 2.5,
		cat.fontfamily = "serif",
		cat.dist = c(0.06, 0.06),
		cat.pos = c(-20, 14),
		euler.d = euler,
		scaled = euler
		);
		
		try(library("grid"), install.packages("grid")) 
    library(grid)
        if (out == TRUE) {
        # only produce ouptut files upon request
        png(filename = paste(file_name,'.png', sep=''))
        grid.draw(output)
        dev.off()       
        pdf(file = paste(file_name,'.pdf', sep=''), width=8, height=8)
        grid.draw(output)
        dev.off()
        } else {
          # just plot the figure
          grid.draw(output)
        }
		
        # pack all the sets up and return to the user
        out <- NULL;
        # a <- c(out, name_a, as.vector(a))
        # a_only <-  as.vector(a_only)
        # a_b <- as.vector(a_b)
        # b <- c(out,name_b,as.vector(b))
        # b_only <- c(out,paste(name_b,"only",sep="_"),as.vector(b_only))
        out <- list(a,a_only,a_b,b,b_only)
        names(out) <- c(name_a, 
                        paste(name_a,"only",sep="_"), 
                        paste(name_a,name_b,sep="_"),
                        name_b,
                        paste(name_b,"only",sep="_"))
  
        return(out)

}