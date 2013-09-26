#check the wrong trees with their distance matrix
#geneating procedule
library(ape)
(mytree <- read.tree(text = "(((sp1:1, sp2:0.1):0.1, sp3:1):
		0.1,((((sp8:1, sp7:1):0.1, sp6:0.1):0.1, sp5:0.1):1, sp4:1):1);"))
source(file.choose()) #newalgorithm.R
tranm <- function(t, alpha){
	P <- diag(4)
	d <- 1/4 + 3/4 * exp(-4*alpha*t)
	x <- rep(1/4 - exp(-4*alpha*t)/4, 6)
	P <- d * P
	P[row(P) > col(P)] <- x
	P[row(P) < col(P)] <- x
	rownames(P) <- colnames(P) <- c('A','C','G','T')
return(P)
}


prob1 <- function(nucleotide){
	if(nucleotide == 'A'){
		p = P1['A', ]
	} else if(nucleotide == 'C'){
		p = P1['C', ]
	} else if(nucleotide == 'G'){
		p = P1['G', ]
	} else{
		p = P1['T', ]
		}
return(p)
}

#this one is for short edge length
prob2 <- function(nucleotide){
	if(nucleotide == 'A'){
		p = P2['A', ]
	} else if(nucleotide == 'C'){
		p = P2['C', ]
	} else if(nucleotide == 'G'){
		p = P2['G', ]
	} else{
		p = P2['T', ]
		}
return(p)
}

#this one is for the edge length of the outgroup
prob3 <- function(nucleotide){
	if(nucleotide == 'A'){
		p = P3['A', ]
	} else if(nucleotide == 'C'){
		p = P3['C', ]
	} else if(nucleotide == 'G'){
		p = P3['G', ]
	} else{
		p = P3['T', ]
		}
return(p)
}

generate <- function(p){
	nucleotide <- sample(c('A', 'C', 'G', 'T'), 1, prob=p)
return(nucleotide)
}

#function to simulate mytree
simulate <- function(){
	proot = c(0.25, 0.25, 0.25, 0.25)
	root  <- generate(proot)
	outgp <- generate(prob3(root))
	node1 <- generate(prob2(root))
	node2 <- generate(prob2(node1))
	tip1  <- generate(prob1(node2))
	tip2  <- generate(prob2(node2))
	tip3  <- generate(prob1(node1))
	node3 <- generate(prob1(root))
	tip4  <- generate(prob1(node3))
	node4 <- generate(prob1(node3))
	tip5  <- generate(prob2(node4))
	node5 <- generate(prob2(node4))
	tip6  <- generate(prob2(node5))
	node6 <- generate(prob2(node5))
	tip7  <- generate(prob1(node6))
	tip8  <- generate(prob1(node6))
	tips  <- c(outgp, tip1, tip2, tip3, tip4, tip5, tip6, tip7, tip8)
return(tips)
}

##function for n numbers of sites of the simulated tree
simn <- function(n){
	ls  <- list()
	out <- as.character()
	sp1 <- as.character()
	sp2 <- as.character()
	sp3 <- as.character()
	sp4 <- as.character()
	sp5 <- as.character()
	sp6 <- as.character()
	sp7 <- as.character()
	sp8 <- as.character()
	for(i in 1:n){
		ls[[i]] <- simulate()
		out[i] <- ls[[i]][1]
		sp1[i] <- ls[[i]][2]
		sp2[i] <- ls[[i]][3]
		sp3[i] <- ls[[i]][4]
		sp4[i] <- ls[[i]][5]
		sp5[i] <- ls[[i]][6]
		sp6[i] <- ls[[i]][7]
		sp7[i] <- ls[[i]][8]
		sp8[i] <- ls[[i]][9]
	}
	dataframe <- rbind(out, sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
return(dataframe)
}

t1 = 1	#long edge length
t2 = 0.1#short edge length
t3 = 4	#edge length for the outgroup
alpha = 0.05
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
dat <- simn(1000) 			#generating dna sequences of the simulating tree
d <- dmat(dat)/alpha		#get the distance matrix
tr <- algorithm(d)			#calculate the tree
all.equal.phylo(unroot(mytree), tr, use.edge.length = F)	#comparing with the true tree

#contour plot
loglik1 <- function(p1, p2, p3, vars){
	n1 = vars[1]
	n2 = vars[2]
	n3 = vars[3]
	n4 = vars[4]
	n5 = vars[5]
	logl =  n1 * log(1/4*((p1*p2*p3)+(1-3*p1)*p2*p3 + p1*(1-3*p2)*p3 + p1*p2*(1-3*p3))) +	#type ACT
			n2 * log(1/4*(1-3*p1)*(1-3*p2)*(1-3*p3) + 3/4*p1*p2*p3) +						#type AAA
			n3 * log(1/4*((1-3*p1)*p2*p3+p1*(1-3*p2)*(1-3*p3)) + 2/4*p1*p2*p3) +			#type ACC
			n4 * log(1/4*((1-3*p1)*(1-3*p2)*p3+p1*p2*(1-3*p3)) + 2/4*p1*p2*p3) +			#type AAC
			n5 * log(1/4*((1-3*p1)*p2*(1-3*p3)+p1*(1-3*p2)*p3) + 2/4*p1*p2*p3)				#type ACA
return(-logl)
}
optim(c(0.15, 0.15, 0.15),  fn = loglik, x = v, 
	lower = 0.0001, upper = 0.249999, method = 'L-BFGS-B')

v <- ext(dat[1,], dat[6,], dat[7,])
v2 <- ext(dat[1,], dat[2,], dat[3,])
v3 <- ext(dat[1,], dat[7,], dat[8,])
x <- y <- seq(0.0001, .24999, by = 0.001)
z <- outer(x, y, FUN = loglik1, p1 = 4*0.05, vars=v)
z2 <- outer(x, y, FUN = loglik1, p1 = 4*0.05, vars=v2)
z3 <- outer(x, y, FUN = loglik1, p1 = 4*0.05, vars=v3)
image(x, y, z)#, xlim = c(0, 0.05), ylim = c(0, 0.05))
contour(x, y, z, nlevels = 50, add = T)

image(x, y, z2, xlim = c(0, 0.05), ylim = c(0, 0.05))
contour(x, y, z2, nlevels = 50, add = T)
image(x, y, z3, xlim = c(0, 0.05), ylim = c(0, 0.05))
contour(x, y, z3, nlevels = 50, add = T)