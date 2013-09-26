###test the optimizer
#source(file.choose())#new algorithm

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
dattest <- simn(1000)
D <- dmat(dattest)


matlist <- list()
for(i in 1:1000){
	dat <- simn(1000)
	matlist[[i]] <- dmat(dat)/alpha
}
matsum <- matlist[[1]]
for(i in 2:1000){
	matsum <- matsum + matlist[[i]]
}
matavg <- matsum/1000

#the actual distance matrix for my true tree
sp1 <- c(0, 1, 1.1, rep(1.2, 5))
sp2 <- c(0.1, 0, 0.2, rep(.3, 5))
sp3 <- c(1, 1, 0, rep(1.1 ,5))
sp4 <- c(rep(2, 3), 0, rep(1, 4))
sp5 <- c(rep(2.1, 3), 1.1, 0, rep(0.1, 3))
sp6 <- c(rep(2.2, 3), 1.2, 0.2, 0, .1, .1)
sp7 <- c(rep(3.2, 3), 2.2, 1.2, 1.1, 0, 1)
sp8 <- c(rep(3.2, 3), 2.2, 1.2, 1.1, 1, 0)
truemat <- rbind(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
colnames(truemat) <- rownames(truemat)

#function to check the mat accuracy.
check <- function(target, current, abstol){
		abs(current - target) <= abstol
}

#test
check(truemat, matavg, 0.02)

meanch <- function(matlist){
	meanmat <- matrix(0, nrow = nrow(matlist[[1]]), ncol = ncol(matlist[[1]]))
	colnames(meanmat) <- rownames(meanmat) <- colnames(matlist[[1]])
	entry <- as.numeric()
	for(i in 1:nrow(meanmat)){
		for(j in 1:ncol(meanmat)){
			for(k in 1:length(matlist)){
				entry[k] <- matlist[[k]][i, j]
				meanmat[i, j] <- mean(entry)
			}
		}
	}
return(meanmat)
}



###this function is to check the variance
varch <- function(matlist){
	varmat <- matrix(0, nrow = nrow(matlist[[1]]), ncol = ncol(matlist[[1]]))
	colnames(varmat) <- rownames(varmat) <- colnames(matlist[[1]])
	entry <- as.numeric()
	for(i in 1:nrow(varmat)){
		for(j in 1:ncol(varmat)){
			for(k in 1:length(matlist)){
				entry[k] <- matlist[[k]][i, j]
				varmat[i, j] <- var(entry)
			}
		}
	}
return(varmat)
}

means <- meanch(matlist)
vars <- varch(matlist)