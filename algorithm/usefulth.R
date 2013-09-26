#####useful functions
#useful packages
library(laser)
library(geiger)
library(ape)
library(apTreeshape)

RandomDna <- function(n, p = c(0.25,0.25,0.25,0.25)) {
    paste(sample(c("A","C","G","T"), n,
		prob=p, replace=TRUE))
}

strsplit(s,"")[[1]] #for some string 's'

#ask for input
x.f <- function(n){
	x <- list()
	for (i in 1:n){
		z <- readline(paste("Parameter ", i, ": ", sep=''))
		# just store the input in a list with the same key
		x[[i]] <- as.numeric(z)
	}
return(x)
}


# in ape package
# nj()       #run neighbor joining using this function
# dist.dna() #once in the bin use this function with jc69 model to get t
# as.DNAbin()#to make dnabin to write into distance matrix form