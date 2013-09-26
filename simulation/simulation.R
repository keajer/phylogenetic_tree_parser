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

#two transition matrix with t1 = 1; t2 = 0.01; alpha = 0.15;
t1 = 1
t2 = 0.01
alpha = 0.15
P1 <- tranm(t1, alpha)
P2 <- tranm(t2, alpha)

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

library(ape)
#my tree:
(mytree <- read.tree(text = "(((sp1:1, sp2:0.01):0.01, sp3:1):
		0.01,((((sp8:1, sp7:1):0.01, sp6:0.01):0.01, sp5:0.01):1, sp4:1):1);"))
plot(mytree)

#start simulate the tree above:
simulate <- function(){
	#choose a root
	proot = c(0.25, 0.25, 0.25, 0.25)
	root  <- generate(proot)
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
	tips  <- c(tip1, tip2, tip3, tip4, tip5, tip6, tip7, tip8)
	#output <- as.character()
	#for(i in 1:8){
	#	output[i] <- paste('Tip', i, 'is: ', tips[i])
	#}
return(tips)
}

##function for n numbers of sites of the simulated tree
simn <- function(n){
	ls  <- list()
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
		sp1[i] <- ls[[i]][1]
		sp2[i] <- ls[[i]][2]
		sp3[i] <- ls[[i]][3]
		sp4[i] <- ls[[i]][4]
		sp5[i] <- ls[[i]][5]
		sp6[i] <- ls[[i]][6]
		sp7[i] <- ls[[i]][7]
		sp8[i] <- ls[[i]][8]
	}
	dataframe <- rbind(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
return(dataframe)
}

dat <- simn(1000) #create 1000 length sequence for 8 sp's
dna <- as.DNAbin(tolower(dat)) #create dnabin
distm <- dist.dna(dna, model = 'JC69') #create a dist matrix w/ model JC69
#hard code the distance matrix
#function to create the dist w/ jc69 algorithm
distjc <- function(seq1, seq2){
	diff <- length(seq1) - sum(seq1 == seq2)
	dis  <- -3/4 * log(1-4/3*diff/length(seq1))
return(dis)
}
dis <- matrix(0, 7, 7)
rownames(dis) <- c('sp2','sp3','sp4','sp5','sp6','sp7','sp8')
colnames(dis) <- c('sp1','sp2','sp3','sp4','sp5','sp6','sp7')
for(i in 1:7)dis[i, 1]   <- distjc(dat[1,], dat[i+1,])
for(i in 1:6)dis[i+1, 2] <- distjc(dat[2,], dat[i+2,])
for(i in 1:5)dis[i+2, 3] <- distjc(dat[3,], dat[i+3,])
for(i in 1:4)dis[i+3, 4] <- distjc(dat[4,], dat[i+4,])
for(i in 1:3)dis[i+4, 5] <- distjc(dat[5,], dat[i+5,])
for(i in 1:2)dis[i+5, 6] <- distjc(dat[6,], dat[i+6,])
dis[7, 7] <- distjc(dat[7,], dat[8,])
#this is identical to the distance matrix made by the dist.dan()
mm <- as.matrix(distm)[2:8, 1:7]
mm[row(mm)<col(mm)] <- 0
all.equal(mm, dis)
#checked

nj1 <- nj(distm) ### build-in NJ algorithm
plot(nj1, 'phy')
#test if same topology as mytree
all.equal.phylo(nj1, mytree, use.edge.length = F)
#they are different, the reason my due to the rooting
all.equal.phylo(nj1, unroot(mytree), use.edge.length = F)
#now, they have the same topology after i unrooted "mytree"


#in case of alpha = 1
#check the max alpha*t, think each site is a burnulli trail
#1000 sites will have a binomial(1000, pij)
#here the tree has the longest distance alpha*t = 4.04
#with alpha = 1 and t = 4.04 we can get
pij <- 1 - (1/4 + 3/4*exp(-4*4.04)) #the prob of a trans 
probs <- as.numeric()
for(i in 751:1000){
	probs[i-750] <- dbinom(i, 1000, pij)
}
sum(probs)
#or we can calculate directly from the cdf
(prob <- 1- pbinom(750, 1000, pij))
#which tells big probablity of having more than 750 site trans
#after some manual test i get the alpha*t should close to 0.71
#to get the prob < tol = 10e-3 with 1000 length in sequence

