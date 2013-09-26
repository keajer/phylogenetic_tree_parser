###algorithm.R
library(ape)

#mytree with an additional outgroup
#(mytree <- read.tree(text = "(outgp:1,((sp1:1, sp2:0.01):0.01, sp3:1):
#		0.01,((((sp8:1, sp7:1):0.01, sp6:0.01):0.01, sp5:0.01):1, sp4:1):1);"))
#plot(mytree)

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
#two transition matrix with t1 = 1; t2 = 0.01; alpha = 0.15;
#t1 = 1
#t2 = 0.01
#alpha = 0.15
#P1 <- tranm(t1, alpha)
#P2 <- tranm(t2, alpha)

#functions to get the probablity of a node
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

#function to generate a nucleotide with probablity from TransMat
generate <- function(p){
	nucleotide <- sample(c('A', 'C', 'G', 'T'), 1, prob=p)
return(nucleotide)
}

#function to simulate mytree
simulate1 <- function(){
	proot = c(0.25, 0.25, 0.25, 0.25)
	root  <- generate(proot)
	outgp <- generate(prob2(root))
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
simn1 <- function(n){
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
		ls[[i]] <- simulate1()
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

#dat <- simn(1000) #create 1000 length sequence for a outgroup and 8 sp's
#will use this data as my test data.

#this function is to extract the info from 3 dna seq.
ext <- function(seq1, seq2, seq3){
	tp1 <- tp2 <- tp3 <- tp4 <- tp5 <- 0
	for(i in 1:length(seq1)){
		if(seq1[i] != seq2[i] & seq2[i] != seq3[i] & seq1[i] != seq3[i]){
			tp1 = tp1 + 1	#type ACT
		}
		else if(seq1[i] == seq2[i] & seq1[i] == seq3[i]){
			tp2 = tp2 + 1	#type AAA
		}
		else if(seq1[i] != seq2[i] & seq2[i] == seq3[i]){
			tp3 = tp3 + 1	#type ACC
		}
		else if(seq1[i] == seq2[i] & seq2[i] != seq3[i]){
			tp4 = tp4 + 1	#type AAC
		}
		else{
			tp5 = tp5 + 1	#type ACA
		}
	}
return(c(tp1, tp2, tp3, tp4, tp5))
}

#fetch the functions for optimization
source("C:\\Users\\Kan\\Desktop\\Thesis\\newalgorithm.R")

###testing
#set the initial values to be initval <- c(0.1,.15,.2)
#initval <- c(0.15,.15,.15)#all p's has constrain 0 <= pi < 0.25, i = 1,2,3
#options(scipen = 20)	#suppress "scientific" notation
#options(scipen = NULL)	#back to default setting
#op1 <- optim(initval,  fn = loglik, gr = loglikgr, x = ext(dat[1, ], dat[2, ], dat[8, ]), #
#		lower = 0.00001, upper = 0.249999, method = 'L-BFGS-B')
#op1$par[2] + op1$par[3]
#op2 <- optim(initval,  fn = loglik, gr = loglikgr, x = ext(dat[1, ], dat[4, ], dat[9, ]), 
#		lower = 0.00001, upper = 0.249999, method = 'L-BFGS-B')
#op2$par[2] + op2$par[3]


#start to write a function getting the initial D distrance matrix
#the function is good for all dna sequences with first row is the outgroup. 
dmat <- function(dna){
	D <- matrix(0, nrow = nrow(dna), ncol = nrow(dna), 
			dimnames = list(rownames(dna), rownames(dna)))
	D <- D[-1, -1]	#here deduct the first column and row for the outgroup
	ls <- list()	#list to store the optimal solutions
	for(k in 1:(nrow(D) - 1)){
		for(i in 1:(ncol(D) - k)){
			ls[[i + k - 1]] <- optim(c(0.15,.15,.15),  fn = loglik, gr = loglikgr, 
							x = ext(dna[1,], dna[1+k,], dna[i+1+k,]), 
							lower = 0.00001, upper = 0.249999, method = 'L-BFGS-B')
			D[k, i + k]	 <- ls[[i + k - 1]]$par[2]
			D[i + k, k]	 <- ls[[i + k - 1]]$par[3]
		}
	}
	D <- -1/4*log(1-4*D)
return(D)
}
#the function returns Tij's which has been transformed by CJ69 model
#such that: T = alpha*t in the model
#D <- dmat(dat)/0.15
#calculate the distance matrix with test data


#this function finds indecies(leaf node names) for a matrix
matdef <- function(d, mini){
	a <- list()
	m <- 1
	for(i in 1:nrow(d)){
		for(j in 1:ncol(d)){
			if(d[i, j] == mini){
				a[[m]] <- c(rownames(d)[i],rownames(d)[j]) 
				m <- m + 1
			}
		}
	}
	a <- sample(a, 1)
return(unlist(a))
}

#this function is to match nodes to its dim# in a special way
#i.e. if node name is from original distance matrix (D) return same dim#
#if node is a internal node, return a backward order + the # of col of D
matchnode <- function(mat1, mat2, node){
	names1 <- colnames(mat1)
	num	<- which(node == names1)
	span <- (ncol(mat1) - ncol(mat2)):1
	if(num <= ncol(mat2)){
		num <- num
	}else{
		num <- ncol(mat2) + span[num - ncol(mat2)]
	}
return(num)
}

#now to write a function for the algorithm
algorithm <- function(D){
	C <- D
	L <- D
	Dim <- nrow(L)
	iter <- 1
	ls <- list()
	while(Dim > 2){ #stop at dim = 2
		#each iteration construct a d-criterion
		d <- matrix(0, ncol = ncol(L), nrow = nrow(L), dimnames = list(rownames(L), rownames(L)))
		for(i in 1:nrow(L)){
			for(j in 1:ncol(L)){
				if(i != j){
					d[i, j] <- L[i, j]/min(L[i, L[i, ] != 0]) + L[j, i]/min(L[j, L[j, ] != 0])
				}
			}
		}
		#with d-criterion, to calculate the minimal pair
		mini <- min(d[d != 0])
		pair <- matdef(d, mini)
		ls[[iter]] <- pair
		
		#now get the new node for C within each iteration
		#C expends one dim each iteration
		C0 <- matrix(0, nrow = (ncol(C)+1), ncol = (ncol(C)+1), 
				dimnames = list(c(colnames(C), paste("n",iter, sep = "")), c(colnames(C), paste("n",iter, sep = ""))))
		C0[1:nrow(C), 1:ncol(C)] <- C
		C <- C0
		iter <- iter + 1
		for(i in rownames(L)){
			if(i != pair[1] & i != pair[2]){
				C[i, ncol(C)] <- 1/2*(C[i, pair[1]] + C[i, pair[2]])
			}
		}
		for(j in colnames(L)){
			if(j != pair[1] & j != pair[2]){
				C[nrow(C), j] <- 1/2*(C[pair[2], j] + C[pair[1], j] - C[pair[1], pair[2]] - C[pair[2], pair[1]])
			}
		}
			
		#allocate L each iteration from C, deduct two dim and and one dim for new node
		#first give row and col names for L
		name <- as.character()
		for(i in colnames(L)){
			if(i != pair[1] & i != pair[2]){
			name[i] <- i  
			}
		}
		name[length(name) + 1] <- colnames(C)[ncol(C)]
		#setup the matrix
		L1 <- matrix(0, nrow = (ncol(L)-1), ncol = (ncol(L)-1), dimnames = list(name, name))
		for(i in rownames(L1)){
			for(j in colnames(L1)){
				L1[i, j] <- C[i, j]
			}
		}
		L <- L1

	Dim <- nrow(L)
	}
	
	#manage the output
	ls[[length(ls)+1]] <- colnames(L)
	ls1 <- unlist(ls)

	m <- matrix(0, ncol = 3, nrow = ncol(D)*2-2)
	for(i in 1:(ncol(D)*2-2)){
		m[i, 1] <- ncol(D)+((ncol(D)*2-2)-i)%/%2
		m[i, 2] <- matchdimn(C, D, ls1[i])
		m[i, 3] <- rexp(1, 5)	#assign a random exp(1) number 
	}
	if(m[nrow(m),2] == m[nrow(m)-1, 1]){
		m[nrow(m)-1,c(1,2)] <- c(m[nrow(m)-1,2],m[nrow(m),2])
		m <- m[-nrow(m),]
	}else{
		m[nrow(m)-1,c(1,2)] <- c(m[nrow(m),2],m[nrow(m)-1,2])
		m <- m[-nrow(m),]
	}
	
	labels <- colnames(D)
	obj <- list(edge = cbind(m[,1], m[,2]), edge.length = m[,3],
				tip.label = labels, Nnode = ncol(C)-ncol(D))
	class(obj) <- "phylo"
	obj <- reorder(obj)

return(obj)
}

#tr1 <- algorithm(D)	#make a tree list

#the actual mytree, i.e. drop the artificial outgroup
#mytree <- drop.tip(mytree, 'outgp')
#all.equal.phylo(unroot(mytree), tr1, use.edge.length = F)

