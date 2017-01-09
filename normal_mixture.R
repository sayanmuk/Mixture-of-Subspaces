#Fit a Mixture of Normal distributions without a subspace assumption
#input is a data matrix X where X is n x m, where m is the ambient dimension, and n is the number of observations
#we assume that this comes from a mixture of normals, and in order to learn classification we will use the mixture of normals
#the number of clusters to be desired is assumed to be the number of classes, which is determined before as K
 
#Load the library necessary for distributional functions
library("LCA")
library("bayesm")
library("MASS") 
library("RcppArmadillo")
library("inline")
library("mnormt")
#Intialize the parameters and also make objects to store the draws from MCMC.
#Prior Quantities
t0 = 1
Sigma0 = diag(m,m)/m
prior.mean = 0
d0 = m+1

#initialize and store the clustering
lat.clus = which(rmultinom(n=n, size = 1, prob=rep(1/K,K))==1,arr.ind=T)[,1]
lat.clus.mat = matrix(0,nrow=iter,ncol=n)
pi.mat = matrix(0,n,K)

#Intialize the multinomial weights
pi.mat = matrix(1/K,nrow=n,ncol=K)
dir.prior = rep(1/K,K)
n.vec = rep(0,K)
r0 = rep(1,K)
#initialize mean vectors, storage
mean.mat = matrix(0,nrow=m,ncol=K)
for(i in 1:K){
  mean.mat[,i] = rnorm(m)
}
mean.mat.store = array(0,dim=c(m,K,iter))

#Initialize variance matrix 
covar.array = array(0,dim=c(m,m,K))
for( i in 1:K){
  covar.array[,,i] = rwishart(nu=m+10,diag(m))$IW
}
covar.array.store = array(0,dim=c(m,m,K,iter))

#Need a function to evaluate likelihood ratios
q.finder = function(pi , Sigma, mu, x){
  return(pi*det(Sigma)^(-1/2)*exp(-1/2*t(t(mu-x))%*%solve(Sigma)%*%t(mu-x)))
}
  
pi.mat = matrix(0,n,K)

print(paste("Starting MCMC at ", Sys.time()))

for(i in 1:iter){
  #update on progress
  if(i%%1000 ==0){print(paste("Iteration ",i," at ", Sys.time()))}
  
  #Sample a mean and covariance for each of the clusters
  for(k in sample(1:K)){
    #for this cluster, find the number of people in this cluster
    n.vec[k] = sum(lat.clus == k)
    
    #Find the specific observations
    X.k = X[lat.clus==k,]
    mean.vec.k = as.vector(apply(X.k,2,mean))
    a.k = n.vec[k]*t0/(n.vec[k]*t0+1)
    post.mean.k = (1-a.k)*prior.mean + a.k*mean.vec.k
    d.k = max(n.vec[k] + d0, 2*m+1)
    #Replace this slow part with C
    #sample.S = matrix(0,m,m)
    #for(j in 1:n.vec[k]){
    #  sample.S = sample.S + t(X.k[j,]-mean.vec.k)%*%t(t((X.k[j,]-mean.vec.k)))
    #}
    sample.S = fastS(nS=n.vec[k],mS=m,XS=as.matrix(X.k),muS=mean.vec.k)
    #end of replacement
    post.S.k = Sigma0 + sample.S + a.k/t0*(mean.vec.k-prior.mean)%*%t(mean.vec.k-prior.mean)
    #Draw the variance parameter
    #print("Draw Covariance")
    #print(X.k)
    covar.array[,,k] = rwishart(d.k-m+1,solve(post.S.k))$IW
    #Draw the mean parameter
    #print("Draw Mean")
    mean.mat[,k] = mvrnorm(n=1,mu = post.mean.k,Sigma = a.k/t0*covar.array[,,k])
    #print("Draw Mean Success")
  }
  
  #Get the dirichlet draw for inclusion probabilities
  post.r = r0+n.vec
  pi.vec = rdirichlet(post.r)
  
  #Sort by inclusion probability, permute the mixture parameters accordingly
  ix= sort(pi.vec,index.return=T)$ix
  pi.vec = pi.vec[ix]
  #temporary parameters
  temp.covar.array = covar.array
  temp.mean.mat = mean.mat
  for(k in 1:K){
  covar.array[,,k]  = temp.covar.array[,,ix[k]]
  }
  for(k in 1:K){
    mean.mat[,k] = temp.mean.mat[,ix[k]]
  }
  #Calculate probability of inclusion for each observation
  #Also replace this slow part in C
  pi.mat = fast.pi.mat(nS=n, mS=m, KS = K, XS = as.matrix(X), muS = mean.mat, SigmaS = covar.array, piS = pi.vec)
  #for(j in 1:n){
   # for(k in 1:K){
    #  pi.mat[j,k] = q.finder(pi = pi.vec[k], Sigma = covar.array[,,k], mu = mean.mat[,k], x = X[j,])
  #  }
  #}
  #Done with replacement
  pi.mat = pi.mat/apply(pi.mat,1,sum)
  #draw latent allocations - this isn't so bad not in C so leave for now
  for(j in 1:n){
    lat.clus[j] = which(rmultinom(n=1, size = 1, prob=pi.mat[j,])==1,arr.ind=T)[,1]
  }
  lat.clus=lat.clus[1:n]
  #store the results of this iteration for analyzing the posterior
  lat.clus.mat[i,]=lat.clus
  mean.mat.store[,,i]=mean.mat
  covar.array.store[,,,i]=covar.array
  
}








