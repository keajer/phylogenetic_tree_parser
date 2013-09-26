#####functions for the new algorithm

#first to use MLE for optimization
###BFGS
#log likelihood function
loglik <- function(pars, x){
	p1 = pars[1]
	p2 = pars[2]
	p3 = pars[3]
	n1 = x[1]
	n2 = x[2]
	n3 = x[3]
	n4 = x[4]
	n5 = x[5]
	logl =  n1 * log(1/4*((p1*p2*p3)+(1-3*p1)*p2*p3 + p1*(1-3*p2)*p3 + p1*p2*(1-3*p3))) +	#type ACT
			n2 * log(1/4*(1-3*p1)*(1-3*p2)*(1-3*p3) + 3/4*p1*p2*p3) +						#type AAA
			n3 * log(1/4*((1-3*p1)*p2*p3+p1*(1-3*p2)*(1-3*p3)) + 2/4*p1*p2*p3) +			#type ACC
			n4 * log(1/4*((1-3*p1)*(1-3*p2)*p3+p1*p2*(1-3*p3)) + 2/4*p1*p2*p3) +			#type AAC
			n5 * log(1/4*((1-3*p1)*p2*(1-3*p3)+p1*(1-3*p2)*p3) + 2/4*p1*p2*p3)				#type ACA
return(-logl)
}


#gradient of loglik
loglikgr <- function(pars, x){
	p1 = pars[1]
	p2 = pars[2]
	p3 = pars[3]
	n1 = x[1]
	n2 = x[2]
	n3 = x[3]
	n4 = x[4]
	n5 = x[5]
	
	dldp1 = (n1*((p2*p3)/2 + (p2*(3*p3 - 1))/4 + (p3*(3*p2 - 1))/4))/
			((p1*p2*(3*p3 - 1))/4 + (p1*p3*(3*p2 - 1))/4 + (p2*p3*(3*p1 - 1))/4 - (p1*p2*p3)/4) - 
			(n3*((p2*p3)/4 - ((3*p2 - 1)*(3*p3 - 1))/4))/
			((p1*(3*p2 - 1)*(3*p3 - 1))/4 - (p2*p3*(3*p1 - 1))/4 + (p1*p2*p3)/2) + 
			(n4*((p2*p3)/2 - (p2*(3*p3 - 1))/4 + (3*p3*(3*p2 - 1))/4))/
			((p3*(3*p1 - 1)*(3*p2 - 1))/4 - (p1*p2*(3*p3 - 1))/4 + (p1*p2*p3)/2) + 
			(n5*((p2*p3)/2 + (3*p2*(3*p3 - 1))/4 - (p3*(3*p2 - 1))/4))/
			((p2*(3*p1 - 1)*(3*p3 - 1))/4 - (p1*p3*(3*p2 - 1))/4 + (p1*p2*p3)/2) + 
			(n2*((3*p2*p3)/4 - (3*(3*p2 - 1)*(3*p3 - 1))/4))/
			((3*p1*p2*p3)/4 - ((3*p1)/4 - 1/4)*(3*p2 - 1)*(3*p3 - 1))
 
	dldp2 = (n1*((p1*p3)/2 + (p1*(3*p3 - 1))/4 + (p3*(3*p1 - 1))/4))/
			((p1*p2*(3*p3 - 1))/4 + (p1*p3*(3*p2 - 1))/4 + (p2*p3*(3*p1 - 1))/4 - (p1*p2*p3)/4) - 
			(n5*((p1*p3)/4 - ((3*p1 - 1)*(3*p3 - 1))/4))/
			((p2*(3*p1 - 1)*(3*p3 - 1))/4 - (p1*p3*(3*p2 - 1))/4 + (p1*p2*p3)/2) + 
			(n3*((p1*p3)/2 + (3*p1*(3*p3 - 1))/4 - (p3*(3*p1 - 1))/4))/
			((p1*(3*p2 - 1)*(3*p3 - 1))/4 - (p2*p3*(3*p1 - 1))/4 + (p1*p2*p3)/2) + 
			(n4*((p1*p3)/2 - (p1*(3*p3 - 1))/4 + (3*p3*(3*p1 - 1))/4))/
			((p3*(3*p1 - 1)*(3*p2 - 1))/4 - (p1*p2*(3*p3 - 1))/4 + (p1*p2*p3)/2) + 
			(n2*((3*p1*p3)/4 - 3*((3*p1)/4 - 1/4)*(3*p3 - 1)))/
			((3*p1*p2*p3)/4 - ((3*p1)/4 - 1/4)*(3*p2 - 1)*(3*p3 - 1))

	dldp3 = (n1*((p1*p2)/2 + (p1*(3*p2 - 1))/4 + (p2*(3*p1 - 1))/4))/
			((p1*p2*(3*p3 - 1))/4 + (p1*p3*(3*p2 - 1))/4 + (p2*p3*(3*p1 - 1))/4 - (p1*p2*p3)/4) - 
			(n4*((p1*p2)/4 - ((3*p1 - 1)*(3*p2 - 1))/4))/
			((p3*(3*p1 - 1)*(3*p2 - 1))/4 - (p1*p2*(3*p3 - 1))/4 + (p1*p2*p3)/2) + 
			(n3*((p1*p2)/2 + (3*p1*(3*p2 - 1))/4 - (p2*(3*p1 - 1))/4))/
			((p1*(3*p2 - 1)*(3*p3 - 1))/4 - (p2*p3*(3*p1 - 1))/4 + (p1*p2*p3)/2) + 
			(n5*((p1*p2)/2 - (p1*(3*p2 - 1))/4 + (3*p2*(3*p1 - 1))/4))/
			((p2*(3*p1 - 1)*(3*p3 - 1))/4 - (p1*p3*(3*p2 - 1))/4 + (p1*p2*p3)/2) + 
			(n2*((3*p1*p2)/4 - 3*((3*p1)/4 - 1/4)*(3*p2 - 1)))/
			((3*p1*p2*p3)/4 - ((3*p1)/4 - 1/4)*(3*p2 - 1)*(3*p3 - 1))
 
return(c(-dldp1, -dldp2, -dldp3))
}

###testing
#mydata  <- c(3,20,3,3,3)  #test data
#initval <- c(0.1,.15,.2)  #has constrain 0 <= pi < 0.25, i = 1,2,3

#get the pars from optim()
#optim(initval, fn = loglik, gr = loglikgr, x = mydata, 
#	lower = 0.00001, upper = 0.249999, method = 'L-BFGS-B')


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


#start to write a function getting the initial D distrance matrix
#the function is good for all dna sequences with first row is the outgroup. 
dmat <- function(dna){
	D <- matrix(0, nrow = nrow(dna), ncol = nrow(dna), 
			dimnames = list(rownames(dna), rownames(dna)))
	D <- D[-1, -1]	#here deduct the first column and row for the outgroup
	ls <- list()	#list to store the optimal solutions
	for(k in 1:(nrow(D) - 1)){
		for(i in 1:(ncol(D) - k)){
			ls[[i + k - 1]] <- optim(c(0.15, 0.15, 0.15),  fn = loglik, gr = loglikgr, 
							x = ext(dna[1,], dna[1+k,], dna[i+1+k,]), 
							lower = 0.0001, upper = 0.249999, method = 'L-BFGS-B')
			D[k, i + k]	 <- ls[[i + k - 1]]$par[2]
			D[i + k, k]	 <- ls[[i + k - 1]]$par[3]
		}
	}
	D <- -1/4*log(1-4*D)
return(D)
}
#the function returns Tij's which has been transformed by CJ69 model
#such that: T = alpha*t in the model


#this function finds indecies(leaf node names) for a specific entry in a matrix
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
		m[i, 2] <- matchnode(C, D, ls1[i])
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
