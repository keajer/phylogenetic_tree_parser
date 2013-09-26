##easy tree test
(mytree <- read.tree(text = "((sp1:0.1, sp2:0.1):0.1, (sp3:0.1, sp4:0.1):0.1);"))
##run a test for the tree above:
#for nj first
alpha <- 0.8
P <- tranm(0.1, alpha)
dat <- simn(1000)
datnj <- dat[-1, ]
dna <- as.DNAbin(tolower(datnj)) #create dnabin
distm <- dist.dna(dna, model = 'JC69')
njtr <- nj(distm)
all.equal.phylo(njtr, unroot(mytree), use.edge.length = F)

#for new algorithm
D <- dmat(dat)/0.8
trnew <- algorithm(D)
all.equal.phylo(trnew, unroot(mytree), use.edge.length = F)
par(mfrow = c(1, 3))
plot(unroot(mytree))
plot(njtr)
plot(trnew)