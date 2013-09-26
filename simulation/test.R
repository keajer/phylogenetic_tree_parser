#a test tree
#just to check how everything works so far
#for a esay tree, NJ can get 100% accuracy.
library(ape)

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

generate <- function(p){
	nucleotide <- sample(c('A', 'C', 'G', 'T'), 1, prob=p)
return(nucleotide)
}

(mytree <- read.tree(text = "((sp1:0.1, sp2:0.1):0.1, (sp3:0.1, sp4:0.1):0.1);"))
plot(mytree)


#test tree
P3 <- tranm(0.1, 1)
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
newtree <- function(){
	proot = c(0.25, 0.25, 0.25, 0.25)
	root  <- generate(proot)
	outgp <- generate(prob3(root))
	node1 <- generate(prob3(root))
	node2 <- generate(prob3(root))
	tip1  <- generate(prob3(node1))
	tip2  <- generate(prob3(node1))
	tip3  <- generate(prob3(node2))
	tip4  <- generate(prob3(node2))
	tips  <- c(outgp, tip1, tip2, tip3, tip4)
return(tips)
}

simn <- function(n){
	ls  <- list()
	outgp <- as.character()
	sp1 <- as.character()
	sp2 <- as.character()
	sp3 <- as.character()
	sp4 <- as.character()
	for(i in 1:n){
		ls[[i]] <- newtree()
		outgp[i] <- ls[[i]][1]
		sp1[i] <- ls[[i]][2]
		sp2[i] <- ls[[i]][3]
		sp3[i] <- ls[[i]][4]
		sp4[i] <- ls[[i]][5]
	}
	dataframe <- rbind(outgp,sp1, sp2, sp3, sp4)
return(dataframe)
}

#function to get a list of tree (drops bad data set)
trtest <- function(n){
	treelist <- list()
	dat		 <- list()
	dna		 <- list()
	distm	 <- list()
	for(i in 1:n){
		distm[[i]] <- NaN
		while(sum(is.nan(distm[[i]])) > 0){
			dat[[i]] <- simn(1000)
			dna[[i]] <- as.DNAbin(tolower(dat[[i]]))
			distm[[i]] <- dist.dna(dna[[i]], model = 'JC69')
		}
		treelist[[i]] <- nj(distm[[i]]) 
	}
return(treelist)
}

accuracy <- function(trls){
	accu <- as.numeric()
	for(i in 1:length(trls)){
		accu[i] <- all.equal.phylo(unroot(mytree), trls[[i]], use.edge.length = F)
	}
	rate <- sum(accu)/length(trls)
return(rate)
}

#here the max t is equal to 0.4, i set the max alpha to equal to 
pij <- 1 - (1/4 + 3/4*exp(-4*0.4*2))
(prob <- 1- pbinom(750, 1000, pij))#0.01362629
#here, i will split alpha to 10 levels upto alpha = 2

alpha <- 0.2
P <- tranm(0.1, alpha)
trls1 <- trtest(1000)

alpha <- 0.4
P <- tranm(0.1, alpha)
trls2 <- trtest(1000)

alpha <- 0.6
P <- tranm(0.1, alpha)
trls3 <- trtest(1000)

alpha <- 0.8
P <- tranm(0.1, alpha)
trls4 <- trtest(1000)

alpha <- 1
P <- tranm(0.1, alpha)
trls5 <- trtest(1000)

alpha <- 1.2
P <- tranm(0.1, alpha)
trls6 <- trtest(1000)

alpha <- 1.4
P <- tranm(0.1, alpha)
trls7 <- trtest(1000)

alpha <- 1.6
P <- tranm(0.1, alpha)
trls8 <- trtest(1000)

alpha <- 1.8
P <- tranm(0.1, alpha)
trls9 <- trtest(1000)

alpha <- 2
P <- tranm(0.1, alpha)
trls10 <- trtest(1000)

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

y <- c(rate1, rate2, rate3, rate4, rate5, rate6, rate7, rate8, rate9, rate10)
x <- c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25)
plot(x, y, type = 'o', xlab = "alpha", ylab = "accuracy rate")
title(main = "Neighbor-Joining accuracy check")


