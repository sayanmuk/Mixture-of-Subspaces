source('embedding_functions.R')

breast = read.csv(file='wdbc.data',header=F)
summary(breast[,2]) #357 cases of benign, 212 cases of malignant
sum(breast=='NA')
#ambient dimension
m=6
#number of components to fit
K =4
#number of observations - make divisible by K
n=1000


#Get generic synthetic data
source('generate_toy_data.R')
#scatterplot3d(X,pch=20,color=rep(c(1,2),each=500),xlab="",ylab="",zlab="")

#Pick number of iterations
iter = 10000
#talk and report
talk = T
report = 100


#Set initial temperature for gibbs posterior. tuning parameter set through cross-validation
temperature  = .1

#metropolis hastings (initial) tuning
sub.met.sd = rep(1,K)

#Pick initial guess of dimensions for each subspace 
#needs to be of length K
R = rep(2,K)

#Find, through metropolis-hastings, the gibbs posterior for the subspaces U and the affine
#translations theta for each component. The loss function is distance to the union of subspaces,
# the random walk is a spherized normal walk. Note: Assumes at this point the data is called 'X'
source('gibbs_posterior_metrop.R')
diff.loss
prop.loss
curr.loss
