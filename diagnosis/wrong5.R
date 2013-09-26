#check the wrong trees with their distance matrix
#geneating procedule
dat <- simn(1000) 			#generating dna sequences of the simulating tree
d <- dmat(dat)/alpha		#get the distance matrix
tr <- algorithm(d)			#calculate the tree
all.equal.phylo(unroot(mytree), tr, use.edge.length = F)	#comparing with the true tree


#####data and functions are stored in the Rdata file.
# 5 wrong trees with their distance matrces
tr1; d1; plot(tr1)
tr2; d2; plot(tr2)
tr3; d3; plot(tr3)
tr4; d4; plot(tr4)
tr5; d5; plot(tr5)

###here are some stored references
#the true rooted tree with its true distance matrix
mytree; truemat; plot(mytree)
t1 = 1	#long edge length
t2 = 0.1#short edge length
t3 = 4	#edge length for the outgroup
alpha = 0.05

#a list of 1000 random generated distance matrix
matlist

#the average mean of 1000 list of random generated distance matrix
#(avgs <- meanch(matlist))
avgs

#this is the variance of that
#(vars <- varch(matlist))
vars

#some main functions i have used here
#dna sequence generater of the tree
simn  #usually i do simn(1000)

#distance matrix generater
dmat

#algorithm function
algorithm



