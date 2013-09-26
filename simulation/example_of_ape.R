#phylogenetic tree example

## Speciation & Diversification ##
# in this session, we will concentrate on statistics and tests that can be run on the tree before we have species traits.

#first, lets load the relevant packages
library(laser)
library(geiger)
library(ape)
library(apTreeshape)

#now, load our tree and plot it, just to see what we have
#mytree<-read.tree("path/to/my/file.tre")
mytree = read.nexus("http://fishlab.ucdavis.edu/parrottree.nex")
plot(mytree)

#lets plot a species through time plot. we will use a plotting command from ape
ltt.plot(mytree)
#under constant birth/death rates, we expect an exponential increase in species numbers with time. this should translate to a straight line on a log scale. two plotting commands enable a qualitative examination of the pattern
ltt.plot(mytree,log='y')
# in ltt.plot (from the ape package), the input is the phylogeny, and the handle log='y' denotes a logarithmic scale for the y axis. this command enables control over titles, axes scale and graphic handles.

# an alternative command from laser produces a similar plot using a vector of branching times. most laser commands uses the vector of branching times as their input. this vector can be obtained using the command branching.times(mytree) 
plotLtt(branching.times(mytree))
#here we get a similar plot of lineage through time. 

# we can also look at the distribution of branching times in our tree
hist(branching.times(mytree))
# another qualitative index is tree balance. a fully bifurcated tree is "balanced" and a fully ladderized tree is "non-balanced". Ic, Colless index can be calculated for the tree using a matrix of # of daughter lineage in the "right" and "left" side of each bifurcation are there on a tree
bal.tree=balance(mytree)

#Clearly, this matrix depends on the number of nodes in the tree. Ic can be normalized to control for the number of species so that an Ic values of 0 represent perfectly bifurcated trees and Ic values of 1 are for ¡°ladderized¡± or comb-like trees 
tree.length=length(mytree$tip.label)
Ic=(sum(abs(bal.tree[,1]-bal.tree[,2])))/(((tree.length-1)/2)*(tree.length-2))
# Ic is useful when running simulations, to account for the effect of tree shape on simulation results.  

# Commonly, people are interested in the processes that generate the tree. the processes underlying speciation are reviewed by Nee (2006 Ann. Rev. Ecol. Evol). by and large, these are birth-death processes in which there are probabilities for speciation and extinction. speciation will result in the bifarcation of a branch, whereas extinction will lead to cutting of the branch.
#it is important to understand that these processes are highly stochastic, and running the same model twice can result in a vastly different result.

### Simulating trees ##
# simulating trees is a powerful feature of R which facilitates simulations and allow test of hypothesis and formulation of null models. Here we will use Geiger's birthdeath.tree function 
# Starting from a root node the function simulates the growth of a phylogenetic tree under a uniform, time-homogeneous birth-death process.  If birth is greater than death, then the number of lineages is expected to grow exponentially. 
# The arguments for the function include birth rate, death rate, the run time (time.stop = after how many years the simulation has to stop), the maximum number of taxa (taxa.stop), and a seed. the latter is useful in case you want to run your simulations again and get the exact same result. lastly, the function can be set up to return trees with 'all extinct' taxa if return.all.extinct=T

# We can start with a pure-birth tree (where d = 0) 
simulated.tree1<-birthdeath.tree(b=0.2, d=0, time.stop=10) 
plot(simulated.tree1) 

# its interesting to note the variation in tree shape and # species:
#we will plot 9 trees on one panel using the par command to prepare the graphic object
par(mfrow=c(3,3))
#then simulate 9 trees under the same parameters
for (i in 1:9){
	simulated.tree1<-birthdeath.tree(b=0.2, d=0, time.stop=10) 
	plot(simulated.tree1)
	}

# I like to think about this process as a repeated drawing of splitting condition with a binomial distribution ("split" vs "dont split") with a probability of b. you can simulate a tree yourself if you repeatedly type 
rbinom(1,1,0.2)
# to get a vector of size 1 with drawings from 0-1 pool, and a probability of 0.2 to get the 1s. all you have to do draw a line moving up, which bifurcates each time you get a 1.


# now do the same, but with a fixed tree size
for (i in 1:9){
	simulated.tree2<-birthdeath.tree(b=0.2, d=0, taxa.stop=15) 
	plot(simulated.tree2)
	}
 	
## note: you can drag the bottom right corner of the Quartz object and change the window size 	
# We can add a death rate and check out the tree. The return.all.extinct=F flag prevents trees with no survivors 
for (i in 1:9){
	simulated.tree2<-birthdeath.tree(b=0.2, d=0.05, time.stop=15, return.all.extinct=F) 
	plot(simulated.tree2)
	}
# note that not all taxa get to the current time- these are extinct taxa, that can be removed with 
simulated.tree3=prune.extinct.taxa(simulated.tree2) 
plot(simulated.tree3)
#its interesting to compare original and pruned trees side by side:
par(mfrow=c(4,2))
for (i in 1:4){
	simulated.tree2<-birthdeath.tree(b=0.2, d=0.05, time.stop=15, return.all.extinct=T) #note that we would get "all extinct" trees.
	plot(simulated.tree2) #to the left: plot original tree
	plot(prune.extinct.taxa(simulated.tree2)) #to the right: pruned tree
	}

# Other Tree Simulations: Trees can also be simulated using constant-rate birth-death with various constraints in TreeSim and PhySim. Random trees can be generated in ape by random splitting of edges (for non-parametric trees) or random clustering of tips (for coalescent trees). phybase can simulate coalescent trees as well.
    
#Different processes are assumed to leave a distinct "footprint" in the distribution of branching times along the tree. an array of methods have been developed to try and find the most likely process that generated the pattern in branching times. we will go over a few of these methods

#Gamma statistic [Pybus and Harvey (2000)] summarizes the information contained in the inter-node intervals of a phylogeny. It is assumed that the tree is ultrametric, and losses power rapidly with incomplete taxon sampling
gammaStat(mytree)
#under the assumption that the clade diversified with constant rates, a normal distribution with mean zero and standard-deviation unity is expected for the gamma statistic (Pybus and Harvey 2000). The null hypothesis that the clade diversified with constant rates is thus compared against this prediction using
2*(1 - pnorm(abs(gammaStat(mytree))))
# this is the probability (p-value) of the data under a constant-rate model for a two-tailed distribution
# a variant in Geiger does it all for you, assuming a one-tailed distribution
gamStat(branching.times(mytree), return.list=TRUE)


# to estimate the birth and death rates, you can use 
birthdeath(mytree)

#the capacity of R to simulate trees makes it a powerful platform for assesing the accuracy of the methods you are using. lete run it ot a tree we simulated before, one with known parameters, and see whether we recover the original parameters.
simulated.tree4<-birthdeath.tree(b=0.2, d=0.05, time.stop=10) 
simulated.tree4=prune.extinct.taxa(simulated.tree4)
plot(simulated.tree4)
birthdeath(simulated.tree4)
# did it recover the true rates?



# to get ML estimates for the extinction fraction  (b/d, the ratio between birth and death rates) and the net diversification rate (b-d).Confidence intervals are more robust than the stdErr.
# Non-linear optimization used to estimate the ML location can be difficult, leading to situations where the algorithm is trapped on local (rather than global) optima. Laser offers a function that tests alternative starting positions for the b/d parameter (termed a in Laser). this is given by the vector ai: 
bd(branching.times(mytree), ai = seq(from=0.05,to=0.95,by=0.05)) 
# the syntax above will attempt to re-start the algorithm 19 times with the rate b/d between 0.05 and 1. in the output we see the ML estimates of a, r1 (net diversification rate, b - d), and the aic score, useful for comparison with other models.
# in the case of the parrotfishes tree, the estimate for b/d is 0, hence the preferred model is expected to be the pure birth model. let¡¯s compare the two:
pureBirth(branching.times(mytree))
#Which model would you prefer?

# other models of speciation are available, including:

# density dependent models, where speciation rate decreases as the number of species in the clade, either with an exponential decay
DDX(branching.times(mytree))

#or with a logarithmic decay in rate
DDL(branching.times(mytree))

# an exponentially declining speciation rate through time and constant extinction
fitSPVAR(branching.times(mytree), init=c(.3, .2, .1))#initial parameters are net diversification rate, k, mu0 [see below]

# exponentially increasing extinction  and constant speciation
fitEXVAR(branching.times(mytree), init=c(.3, .01, 1)) #initial parameters are net diversification rate, mu0, z [see below]

#both speciation and extinction rates vary through time.
fitBOTHVAR(branching.times(mytree), init=c(.3, .5, .1, .5)) #initial parameters are net diversification rate, k, mu0, z [see below]

# in the above 3 models, the output includes the speciation rate, specified by parameters lam0 and k, and extinction through time is described by mu0 and z. lam0 and mu0 are the initial speciation and final extinction rates, respectively. k and z control the rate of decrease / increase in speciation and extinction, respectively and are 1 when rates are constant. 
#we set up the initial conditions for the algorithm to start searching. it is a good practice to re-iterate the search starting at different points. A starting point might be the parameter values under the pure birth or constant rate birth death model.

# we can easily compare these models using AIC scores, or use 
fitdAICrc(branching.times(mytree), modelset = c("pureBirth", "bd", "DDL", "DDX"), ints = NULL)
#to compare the nested models.




## detecting shifts in diversification rates using MEDUSA
#MEDUSA ?ts a series of diversi?cation models to a combination of phylogenetic and taxonomic data to use in cases of incomplete sampling (some other functions in Laser will calculate birth and death rates for such data). The input is a phylogenetic tree with branch lengths proportional to time showing the relationship among clades, and the diversity for living species in all of those clades. All taxa missing from the tree have to be assigned to one of the tip clades in the richness matrix. 
# we will use a mock taxonomic dataset:
richness=read.table('richness.txt',header=F,skip=1)

#The algorithm behind MEDUSA ?rst ?ts a Single diversification model to the entire dataset. Then, a series of breaks are added, so that different parts of the tree evolve with different parameter values (birth and or death rates). MEDUSA first compares all single-breakpoint models with the overall model, and selects the best one. Then all possible two-breakpoint models are compared with the best single-breakpoint model, and so on. this could take some time on a large tree... in the interest of our time, we will run up to 5 points

#arguments in runMedusa include estimateExtinction - If false, fits a series of pure-birth models; otherwise birth-death models are used. modelLimit- Maximum number of rate shifts to place on the tree, and cutAtStem- When rate shifts are placed on branches in Medusa, they could be placed at either the beginning or the end of the stem branch. cutAtStem=T cuts at the beginning of the stem branch, while cutAtStem=F cuts at the end.

out=runMedusa(mytree, richness, modelLimit=5)
# as the function progress, you will see the number of breakpoint currenlty calculated
out
# Output: a matrix where each row summarizes the next more complex model selected by the stepwise procedure. Columns record: the node where the next break was found, the likelihood, the number of parameters, the AIC score, and the AICc score (corrected for small sample size).

# I experienced a problem with the original package so I wrote the author (Luke Harmon) who sent me an update; download and save  the three .r files from  http://fishlab.ucdavis.edu/links.html and save them anywhere you want. then, we will source these files. this is a procedure that helps R find the source file for a command that you or someone else created. like the 'library" command, it is valid only for the current session 
source('/Users/Roi/Documents/Comparative_methods/R_source/medusa.R')
source('/Users/Roi/Documents/Comparative_methods/R_source/fitDiscrete.R')
source('/Users/Roi/Documents/Comparative_methods/R_source/node.leaves.R')

summaryMedusa(mytree, r,out)
#this command extracts the parameters for the selected model, providing likelihood estimates for each group, two rate parameters net diversification rate (b-d) and extinction fraction (d/b), and the number of taxa in the group. it also provides a figure of the tree with paintings of the various groups
 

############brownian trait simulations 
speciation=rbinom(1000,1,0.01)
species=1+cumsum(speciation)
plot(species)

bm.var=2.5
time=(1:1000)

states=matrix(NA,1000,species[1000])
states[1,1]=0

plot(time[1],states[1,1],xlim=c(0,1000),ylim=c(-60,60),pch=".", cex = 2.5)

for (i in 1:999){
	if (speciation[i]==1){
		a=sample(species[i]-1,1)
		print(a)
		states[i,species[i]]=states[i,a]
		}
	states[1+i,1:species[i]]=states[i,1:species[i]]+((rbinom(species[i],1,0.5)-0.5)*bm.var)
	
	for (r in 1:species[i]){
		#points(time[i],states[i,r],pch=".",col=r, cex = 2.5)
		lines(c(time[i-1],time[i]),c(states[i-1,r],states[i,r]),pch=".",col=r, cex = 2.5)
		}
		#slow down simulation
	for (z in 1:(10^5)){}
	}

#plot(time,states[,1],xlim=c(0,1000),ylim=c(min(states,na.rm=T),max(states,na.rm=T)),pch=".", cex = 2.5)
#for (r in 1:species[i]){
#	points(time,states[,r],pch=".",col=r, cex = 2.5)
#	}



