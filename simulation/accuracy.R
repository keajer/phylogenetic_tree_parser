#####accuricy comparison
###check the accuracy for both NJ algorithm and new algorithm
#on 10 levels of alpha, comparing 1000 NJ trees on the true tree 
#versus 1000 new algorithm tree on the same true tree.
#each simulated dataset is based on 1000 sites dna sequence in length
#from a hard coded tree which follows the true tree under JC69 model.
library(ape)#package will be used for NJ tree
#source the functions for the new algorithm
source(file.choose())#newalgorithm.R
###start comparison
#here is the true tree we will compare on:
(mytree <- read.tree(text = "(((sp1:1, sp2:0.1):0.1, sp3:1):
		0.1,((((sp8:1, sp7:1):0.1, sp6:0.1):0.1, sp5:0.1):1, sp4:1):1);"))
plot(mytree)

#function to create transition matrix for Jukes and Cantor model
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

###functions to get the probablity of a node
#this one is for long edge length
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


#function to generate a nucleotide with probablity from TransMat
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


#function to get a list of tree (drops bad data set)
trtest <- function(n){
	njtree 	 <- list()
	newtree	 <- list()
	dat		 <- list()
	datnj	 <- list()
	dna		 <- list()
	distm	 <- list()
	D		 <- list()
	treelist <- list()
	for(i in 1:n){
		distm[[i]] <- NaN
		while(sum(is.nan(distm[[i]])) > 0 ){
			dat[[i]] <- simn(1000)
			D[[i]] <- dmat(dat[[i]])
			datnj[[i]] <- dat[[i]][-1, ]#drop the first row of the outgroup for NJ
			dna[[i]] <- as.DNAbin(tolower(datnj[[i]]))
			distm[[i]] <- dist.dna(dna[[i]], model = 'JC69')
		}
		njtree[[i]] <- nj(distm[[i]]) #nj tree
		newtree[[i]] <- algorithm(D[[i]])
	treelist[[i]] <- list(NJ = njtree[[i]], NEW = newtree[[i]])
	}
return(treelist)
}

#as tested in simulation.R previously, the alpha*t should 
#be less than 0.71 to get less than 10e-3 tolence of 
#the probablity for more than 3/4 difference in 1000 
#sites. with the max t = 4.4, our alpha should be around
#0.1613636, for testing accuracy of neighbor joining, i will
#increase alpha to be 0.25, so the prob. for exceed 3/4 is
pij <- 1 - (1/4 + 3/4*exp(-4*4.4*0.25))
(prob <- 1- pbinom(750, 1000, pij))#0.2427129
#in this case we will probably need to throw out some bad 
#data set from simulation, bad here means over 3/4 difference
#within 1000 sites for any of two sequences.
#split 0.25 into 10 levels of alpha, we will have
# alpha = 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2
#	0.225, 0.25

#for some reason R cannot pass variables for more than 
#3 funcions, so i didn't write a new function giving the
#tree list, but implemented this manually
##############################################
#here I will fix the edge length for the tree
##############################################
t1 = 1	#long edge length
t2 = 0.1#short edge length
t3 = 4	#edge length for the outgroup
#just want to be relistic ^ ^

#####ALPHA******
alpha = 0.025
#####change above alpha each time to get tree list
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls1 <- trtest(1000)
#####ALPHA******
alpha = 0.05
#####change above alpha each time to get tree list
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls2 <- trtest(10)
#####ALPHA******
alpha = 0.075
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls3 <- trtest(1000)
#####ALPHA******
alpha = 0.1
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls4 <- trtest(1000)
#####ALPHA******
alpha = 0.125
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls5 <- trtest(1000)
#####ALPHA******
alpha = 0.15
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls6 <- trtest(1000)
#####ALPHA******
alpha = 0.175
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls7 <- trtest(1000)
#####ALPHA******
alpha = 0.2
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls8 <- trtest(1000)
#####ALPHA******
alpha = 0.225
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls9 <- trtest(1000)
#####ALPHA******
alpha = 0.25
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)
P3 <- tranm(t3, alpha)
trls10 <- trtest(1000)

#once we get 10 list of 1000 nj and new trees, we will 
#try to compare every single tree with an unrooted true tree

#then calculate the percentage of acuuracy, and make
#a plot over each level of alpha
#this function is to get the accuracy rate for each list of 1000 trees
accuracy <- function(trls){
	accunj  <- as.numeric()
	accunew <- as.numeric()
	for(i in 1:length(trls)){
		accunj[i]  <- all.equal.phylo(unroot(mytree), trls[[i]][[1]], use.edge.length = F)
		accunew[i] <- all.equal.phylo(unroot(mytree), trls[[i]][[2]], use.edge.length = F)
	}
	ratenj  <- sum(accunj)/length(trls)
	ratenew <- sum(accunew)/length(trls)
return(c(ratenj, ratenew))
}

rate1 <- accuracy(trls1)
rate2 <- accuracy(trls2)
rate3 <- accuracy(trls3)
rate4 <- accuracy(trls4)
rate5 <- accuracy(trls5)
rate6 <- accuracy(trls6)
rate7 <- accuracy(trls7)
rate8 <- accuracy(trls8)
rate9 <- accuracy(trls9)
rate10 <- accuracy(trls10)

ynj  <- c(rate1[1], rate2[1], rate3[1], rate4[1], rate5[1], rate6[1], rate7[1], rate8[1], rate9[1], rate10[1])
ynew <- c(rate1[2], rate2[2], rate3[2], rate4[2], rate5[2], rate6[2], rate7[2], rate8[2], rate9[2], rate10[2])
x <- c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25)
plot(x, ynj, type = 'o', xlab = "alpha", ylab = "accuracy rate")
lines(x, ynew, type = 'o', col = 2)
title(main = "Neighbor-Joining vs. New Algorithm")


optim(c(0.15, 0.15, 0.15),  fn = loglik, gr = loglikgr, 
x = ext(dat[1,], dat[3,], dat[7,]), 
lower = 0.00001, upper = 0.249999, method = 'L-BFGS-B')

optim(c(0.15, 0.15, 0.15),  fn = loglik, gr = loglikgr, 
x = ext(dat[1,], dat[3,], dat[7,]), 
lower = 0.00001, upper = 0.249999, method = 'L-BFGS-B')