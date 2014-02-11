venn_diagram4 <- function(a, b, c, d, name_a = "A", name_b = "B", name_c = "C", name_d = "D", colors =  c("#542BFD","#FCFF4F","#5CFF4A", "#F63C4B"), out = FALSE, file_name = "default.pdf", euler=FALSE, special=TRUE){
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
	d <- intersect(d,d)
	# length of a
	a_tot <- length(a)
	b_tot <- length(b)
	c_tot <- length(c)
	d_tot <- length(d)

	# number shared between a and b and c and d
	a_b_c_d <- intersect(intersect(a,b),intersect(c,d))
	a_b_c_d_tot <- length(a_b_c_d)
	a_b_c_d_tot
	
	# a, b, and c, check to make sure not in d
	a_b_c <- setdiff(intersect(intersect(a,b), c), d)
	a_b_c_tot <- length(a_b_c)
	
	# b,c, and d, checking to make sure not in a
	b_c_d <-setdiff(intersect(intersect(b,c), d), a)
	b_c_d_tot <- length(b_c_d)
	
	# a,c, and d, checking to make sure not in b
	a_c_d <- setdiff(intersect(intersect(a,c), d),b)
	a_c_d_tot <- length(a_c_d)
	
	# a,b, and d, checking to make sure not in c
	a_b_d <- setdiff(intersect(intersect(a,b),d),c)
	a_b_d_tot <- length(a_b_d)
	
	# all pair sets only
	a_b <- setdiff(intersect(a,b),union(c,d))
	a_b_tot <- length(a_b)
	a_c <- setdiff(intersect(a,c),union(b,d))
	a_c_tot <- length(a_c)
	a_d <- setdiff(intersect(a,d),union(b,c))
	a_d_tot <- length(a_d)
	b_c <- setdiff(intersect(b,c),union(a,d))
	b_c_tot <- length(b_c)
	b_d <- setdiff(intersect(b,d),union(a,c))
	b_d_tot <- length(b_d)
	c_d <- setdiff(intersect(c,d),union(a,b))
	c_d_tot <- length(c_d)
	
	# a, b, c, d only
	a_only <- setdiff(a,union(union(b,c),d))
	a_only_tot <- length(a_only)
	b_only <- setdiff(b,union(union(a,c),d))
	b_only_tot <- length(b_only)
	c_only <- setdiff(c,union(union(a,b),d))
	c_only_tot <- length(c_only)
	d_only <- setdiff(d,union(union(a,b),c))
	d_only_tot <- length(d_only)
	
	# sanity checks
	a_tot == (a_only_tot + a_b_tot + a_b_c_tot + a_b_c_d_tot + a_b_d_tot + a_d_tot + a_c_d_tot + a_c_tot)
	b_tot == (b_only_tot + b_c_tot + b_c_d_tot + b_d_tot + a_b_d_tot + a_b_c_d_tot + a_b_c_tot + a_b_tot)
	c_tot == (c_only_tot + c_d_tot + b_c_d_tot + a_b_c_d_tot + a_c_d_tot + a_c_tot + a_b_c_tot + b_c_tot)
	d_tot == (d_only_tot + b_d_tot + a_b_d_tot + a_d_tot + a_c_d_tot + a_b_c_d_tot + b_c_d_tot + c_d_tot)
	
	# try to load VennDiagram else try to download it
	# More info: see Hanbo Chen and Paul C Boutros. BMC Bioinformatics 2011, 12:35 doi:10.1186/1471-2105-12-35
	try(library(VennDiagram), install.packages("VennDiagram")) 
	library(VennDiagram)
	

	# create the first list
	offset_a <- 1*(10^6)
	offset_b <- 2*(10^6)
	offset_c <- 3*(10^6)
	offset_d <- 4*(10^6)
	offset_a_b <- 5*(10^6)
	offset_a_c <- 6*(10^6)
	offset_a_d <- 7*(10^6)
	offset_b_c <- 8*(10^6)
	offset_b_d <- 9*(10^6)
	offset_c_d <- 10*(10^6)
	offset_a_b_c <- 11*(10^6)
	offset_b_c_d <- 12*(10^6)
	offset_a_c_d <- 13*(10^6)
	offset_a_b_d <- 14*(10^6)
	offset_a_b_c_d <- 15*(10^6)
	
	list_a_only <- NULL
	list_b_only <- NULL
	list_c_only <- NULL
	list_d_only <- NULL
	list_a_b <- NULL
	list_a_c <- NULL
	list_a_d <- NULL
	list_b_c <- NULL
	list_b_d <- NULL
	list_c_d <- NULL
	list_a_b_c <- NULL
	list_b_c_d <- NULL
	list_a_c_d <- NULL
	list_a_b_d <- NULL
	list_a_b_c_d <- NULL
	
	if(length(a_only) > 0){
	    list_a_only <- offset_a:(offset_a+length(a_only)-1);
	}
	if(length(b_only) > 0){
	    list_b_only <- offset_b:(offset_b+length(b_only)-1);
	}
	if(length(c_only) > 0){
	    list_c_only <- offset_c:(offset_c+length(c_only)-1);
	}
	if(length(d_only) > 0){
	    list_d_only <- offset_d:(offset_d+length(d_only)-1);
	}			
    if(length(a_b) > 0){
        list_a_b <- offset_a_b:(offset_a_b+length(a_b)-1);
    }	
    if(length(a_c) > 0){
        list_a_c <- offset_a_c:(offset_a_c+length(a_c)-1);
    }
    if(length(a_d) > 0){
        list_a_d <- offset_a_d:(offset_a_d+length(a_d)-1);
    }
    if(length(b_c) > 0){
        list_b_c <- offset_b_c:(offset_b_c+length(b_c)-1);
    }
    if(length(b_d) > 0){
        list_b_d <- offset_b_d:(offset_b_d+length(b_d)-1);
    }	
    if(length(c_d) > 0){
        list_c_d <- offset_c_d:(offset_c_d+length(c_d)-1);
    }
    if(length(a_b_c) > 0){
        list_a_b_c <- offset_a_b_c:(offset_a_b_c+length(a_b_c)-1);
    }
    if(length(b_c_d) > 0){
        list_b_c_d <- offset_b_c_d:(offset_b_c_d+length(b_c_d)-1);
    }
    if(length(a_c_d) > 0){
        list_a_c_d <- offset_a_c_d:(offset_a_c_d+length(a_c_d)-1);
    }
    if(length(a_b_d) > 0){
        list_a_b_d <- offset_a_b_d:(offset_a_b_d+length(a_b_d)-1);
    }
	if(length(a_b_c_d) > 0){
		list_a_b_c_d <- offset_a_b_c_d:(offset_a_b_c_d+length(a_b_c_d)-1);
	}
	
	temp <- list(
		name_a = c(list_a_only, 
		          list_a_b,
		          list_a_b_c,
                  list_a_c,
				  list_a_c_d,
				  list_a_d,
                  list_a_b_c_d,				
				  list_a_b_d
                  ),
		name_d = c(list_d_only, 
		           list_c_d, 
		           list_b_c_d, 
		           list_b_d,
				   list_a_b_c_d,
				   list_a_b_d,
				   list_a_c_d,
				   list_a_d
		          ),
		name_b = c(list_a_b, 
		           list_a_b_c,
		           list_a_b_c_d,
		           list_a_b_d,
				   list_b_only,
				   list_b_c,
				   list_b_c_d,
				   list_b_d
		           ),
		name_c = c(list_c_only, 
		           list_b_c,
		           list_a_b_c,
		           list_a_c,
				   list_a_c_d,
				   list_a_b_c_d,
				   list_b_c_d,
				   list_c_d
		           )		
		)
		

	names(temp) <- c(name_a,name_d,name_b,name_c)
	
	output<-venn.diagram(
		x = temp,
		sp.cases = special,
		filename = NULL,
		col = "black",
		fill = c(colors[1], colors[4], colors[2], colors[3]),
		alpha = 0.5,
		label.col = c("black"),
		cex = 2.5,
		fontfamily = "serif",
		fontface = "bold",
		cat.default.pos = "text",
		cat.col = c("black"),
		cat.cex = 3.0,
		cat.fontfamily = "serif",
		cat.pos = 0,
		euler.d = euler,
		scaled = euler
		);

    try(library("grid"), install.packages("grid")) 
    library(grid)
    if (out == TRUE) {
       png(filename = paste(name_output,'.png', sep=''))
       grid.draw(output)
       dev.off()       
       pdf(file = paste(name_output,'.pdf', sep=''), width=8, height=8)
       grid.draw(output)
       dev.off()
    } else {
      # just plot the graph
      grid.draw(output)
    }
		
		# pack all the sets up and return to the user
		out <- NULL;
		out <- list(a,a_only,a_b,a_c,a_d,
                b,b_only,b_c,b_d,
                c,c_only,c_d,
                d,d_only,
                a_b_c,b_c_d,a_c_d,a_b_d,
                a_b_c_d);
  
    names(out) <- c(name_a,paste(name_a,"only",sep="_"),
                    paste(name_a,name_b,sep="_"),
                    paste(name_a,name_c,sep="_"),
                    paste(name_a,name_d,sep="_"),
                    name_b,
                    paste(name_b,"only",sep="_"),
                    paste(name_b,name_c,sep="_"),
                    paste(name_b,name_d,sep="_"),
                    name_c,
                    paste(name_c,"only",sep="_"),
                    paste(name_c,name_d,sep="_"),
                    name_d,
                    paste(name_d,"only",sep="_"),
                    paste(name_a,name_b,name_c,sep="_"),
                    paste(name_b,name_c,name_d,sep="_"),
                    paste(name_a,name_c,name_d,sep="_"),
                    paste(name_a,name_b,name_d,sep="_"),
                    paste(name_a,name_b,name_c,name_d,sep="_"))
		return(out)
}