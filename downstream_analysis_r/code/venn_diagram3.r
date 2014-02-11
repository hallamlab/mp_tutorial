venn_diagram3 <- function(a, b, c, name_a = "A", name_b = "B", name_c = "C", colors =  c("#B80830","#EACC33","#46E2D9"), out = FALSE, file_name = "default.pdf", euler=FALSE, special=TRUE){
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
	c <- intersect(c,c)
	# length of a
	a_tot <- length(a);
	a_tot
	# number shared between a and b and c
	a_b_c <- intersect(a,intersect(b,c))
	a_b_c_tot <- length(a_b_c)
	a_b_c_tot
	# number shared with b
	a_b <- setdiff(intersect(a,b),a_b_c)
	a_b_tot <- length(a_b)
	a_b_tot
	# number shared with c
	a_c <- setdiff(intersect(a,c),a_b_c)
	a_c_tot <- length(a_c)
	a_c_tot
	# number in a only
	a_only <- setdiff(a, union(union(a_c,a_b), a_b_c))
	a_only_tot <- length(a_only)
	a_only_tot
	# we should have defined all classificaitons for set 'a'
	# sanity check should resolve TRUE if everything adds up
	a_tot == (a_b_c_tot + a_b_tot + a_c_tot + a_only_tot)
	
	# B: do the same for intersecting classes for b
	b_tot <- length(b)
	b_tot
	# number shared between b and c
	b_c <- setdiff(intersect(b,c),a_b_c)
	b_c_tot <- length(b_c)
	b_c_tot
	# b only
	b_only <- setdiff(b, union(union(b_c,a_b), a_b_c))
	b_only_tot <- length(b_only)
	b_only_tot
	
	# sanity check for B
	b_tot == (a_b_tot + a_b_c_tot + b_c_tot + b_only_tot)
	
	# C: do the same for intersecting classes for c
	c_tot <- length(c)
	c_tot
	c_only <- setdiff(c,union(union(b_c,a_c), a_b_c))
	c_only_tot <- length(c_only)
	c_only_tot
	
	# sanity check
	c_tot == (a_c_tot + b_c_tot + a_b_c_tot + c_only_tot)
	# try to load VennDiagram else try to download it
	# More info: see Hanbo Chen and Paul C Boutros. BMC Bioinformatics 2011, 12:35 doi:10.1186/1471-2105-12-35
	try(library(VennDiagram), install.packages("VennDiagram")) 
	library(VennDiagram)
	

	# create the first list
	offset_a <- 1*(10^6)
	offset_a_b <- 2*(10^6)
	offset_a_c <- 3*(10^6)
	offset_b <- 4*(10^6)
	offset_b_c <- 5*(10^6)
	offset_c <- 6*(10^6)
	offset_a_b_c <- 7*(10^6)
	
	list_a_only <- NULL
	list_a_b <- NULL
	list_a_c <- NULL
	list_b_only <- NULL
	list_b_c <- NULL
	list_c_only <- NULL
	list_a_b_c <- NULL
	
	if(length(a_only) > 0){
	    list_a_only <- offset_a:(offset_a+length(a_only)-1);
	}
    if(length(a_c) > 0){
        list_a_c <- offset_a_c:(offset_a_c+length(a_c)-1);
    }	
    if(length(a_b) > 0){
        list_a_b <- offset_a_b:(offset_a_b+length(a_b)-1);
    }	
	if(length(b_only) > 0){
	    list_b_only <- offset_b:(offset_b+length(b_only)-1);
	}	
    if(length(b_c) > 0){
        list_b_c <- offset_b_c:(offset_b_c+length(b_c)-1);
    }	
	if(length(c_only) > 0){
	    list_c_only <- offset_c:(offset_c+length(c_only)-1);
    }
    if(length(a_b_c) > 0){
        list_a_b_c <- offset_a_b_c:(offset_a_b_c+length(a_b_c)-1)
    }
	
	
	temp <- list(
		name_a = c(list_a_only, 
		          list_a_b,
		          list_a_b_c,
                  list_a_c),
		name_b = c(list_b_only, 
		           list_a_b, 
		           list_a_b_c, 
		           list_b_c),
		name_c = c(list_c_only, 
		           list_a_c,
		           list_a_b_c,
		           list_b_c
		           )
		)
		
        # temp <- list(
        #           name_a = c(1:(length(a_only)), 
        #                     (length(a_only)):(length(a_only)-1+length(a_b)),
        #                     (length(a_only)+length(a_b)+1):(length(a_only)+length(a_b)+length(a_b_c)),
        #                       (length(a_only)+length(a_b)+length(a_b_c)+1):(length(a_only)+length(a_b)+length(a_b_c)+length(a_c))),
        #           name_b = c((a_tot+1):(a_tot+length(b_only)), 
        #                      (length(a_only)+1):(length(a_only)+length(a_b)), 
        #                      (length(a_only)+length(a_b)+1):(length(a_only)+length(a_b)+length(a_b_c)), 
        #                      (a_tot+length(b_only)+1):(a_tot+length(b_only)+length(b_c))),
        #           name_c = c((c_ind+1):(c_ind+1+length(c_only)), 
        #                      (length(a_only)+length(a_b)+1):(length(a_only)+length(a_b)+length(a_b_c)),
        #                      (length(a_only)+length(a_b)+length(a_b_c)+1):(length(a_only)+length(a_b)+length(a_b_c)+length(a_c)),
        #                      (a_tot+length(b_only)+1):(a_tot+length(b_only)+length(b_c)))
        #           )
	names(temp) <- c(name_a,name_b,name_c)
	
	output<-venn.diagram(
		x = temp,
		sp.cases = special,
		filename = NULL,
		col = "black",
		fill = c(colors[1], colors[2], colors[3]),
		alpha = 0.5,
		label.col = c("black", "white", "black", "white", "white", "white", "black"),
		cex = 2.5,
		fontfamily = "serif",
		fontface = "bold",
		cat.default.pos = "text",
		cat.col = c("black", "black", "black"),
		cat.cex = 2.5,
		cat.fontfamily = "serif",
		cat.dist = c(0.06, 0.06, 0.03),
		cat.pos = 0,
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
      pdf(file = file_name, width=8, height=8)
      grid.draw(output)
      dev.off()
    } else {
      # just plot the figure
      grid.draw(output)
    }
        
		
	out <- list(a,
              a_only,
              a_b,
              a_c,
	            b,
              b_only,
              b_c,
	            c,
              c_only,
	            a_b_c)
	
	names(out) <- c(name_a,
                  paste(name_a,"only",sep="_"),
	                paste(name_a,name_b,sep="_"),
	                paste(name_a,name_c,sep="_"),
	                name_b,
	                paste(name_b,"only",sep="_"),
	                paste(name_b,name_c,sep="_"),
	                name_c,
	                paste(name_c,"only",sep="_"),
	                paste(name_a,name_b,name_c,sep="_"))
	return(out)
}