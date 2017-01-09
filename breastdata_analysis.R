source('embedding_functions.R')

breast = read.csv(file='wdbc.data',header=F)
summary(breast[,2]) #357 cases of benign, 212 cases of malignant
for(i in 3:32){
  breast[,i]= (breast[,i]-mean(breast[,i]))/sd(breast[,i])
}
X = breast[,3:32]
#ambient dimension
m=dim(X)[2]
#number of components to fit
K =2
#number of observations 
n=dim(X)[1]
hold.n = sample(1:n,floor(n/10),replace=F)

#Find true clustering for classification
true.cluss = as.numeric(breast[,2])-1; true.cluss = true.cluss[-hold.n]; true.hold.cluss = true.cluss[hold.n]
theta.true.list = vector("list",2)
theta.true.list[[1]] = apply(X[true.cluss==0,],2,mean)
theta.true.list[[2]] = apply(X[true.cluss==1,],2,mean)
#Pick number of iterations
iter = 10000
#talk and report
talk = T
report = 100
X = X[-hold.n,]
n = n - length(hold.n)
total.n = n + length(hold.n)

#Set initial temperature for gibbs posterior. Is this a tuning parameter? Too expensive to set through cross-validation
temperature  = 10^-7

#metropolis hastings (initial) tuning
sub.met.sd = rep(1,K)

#Pick initial guess of dimensions for each subspace (shouldnt really matter, but being close helps)
#needs to be of length K
R = rep(2,K)

#Find, through metropolis-hastings, the gibbs posterior for the subspaces U and the affine
#translations theta for each component. The loss function is distance to the union of subspaces,
# the random walk is a spherized normal walk. Note: Assumes at this point the data is called 'X'
source('gibbs_posterior_metrop.R')
#P.mat[[i]] = P.list  - access ith iteration, kth cluster iteration with P.mat[[i]][[k]]
#U.mat[[i]] = U.list
#theta.mat[[i]]=theta.list

fitted.submix = predict.submix(P = P.mat, theta = theta.mat, X = breast[hold.n,3:32],K=K)
submix.prediction =apply(apply(fitted.submix$cluster, 2, table),2,which.max)
1-sum(1+as.numeric(breast[,2])[hold.n]==submix.prediction)/length(hold.n) #.89

clus.mat = matrix(0,nrow=iter,ncol=n)
for(i in 1:iter){
  if(talk == T){
    if(i%%report==0){
      print(paste('iteration ', i, 'at ', Sys.time()))
      
    }
  }
  projections = P.mat[[i]]
  theta.mus = theta.mat[[i]]
  for(l in 1:K){
    if(is.null(dim(theta.mus[[l]]))){
      theta.mus[[l]]= matrix(theta.mus[[l]],nrow=m,ncol=1)
    }
  }
  clus.this = fast.log.loss(x=X,P.list = projections, mu = theta.mus,temperature=temperature)$clus
  clus.mat[i,]=clus.this 
}

clus.est=apply(clus.mat,2,Mode)
sum(!clus.est==true.cluss)/n #.872 // .875 // .971  - !clus.est because the components aren't necessarily ordered the same way
fit1 = glm(true.cluss ~ as.matrix(X), family="binomial") 
1-sum(round(fitted.values(fit1))==true.cluss)/n #.855
sum(round(predict(object = fit1, newdata = breast[hold.n,3:32],type="response" ))==breast[,2][hold.n])/length(hold.n)

#save.image(file='breast_analysis.Rdata')
#load(file='breast_analysis.Rdata')

#Calculate a non-subspace assumption normal mixture model
source('mixture_normal_functions.R')
source('normal_mixture.R')
#Objects are stored in lat.clus.mat, mean.mat.stor, covar.array.store
#save.image(file='breast_analysis.Rdata')
#load(file='breast_analysis.Rdata')
lat.clus.mean = apply(lat.clus.mat,2, mean)
lat.clus.est = round(lat.clus.mean)-1
sum(!lat.clus.est==true.cluss)/n #.694, .694, /// .682 -// .665 // .674 // .674

#Make predictions for the clustering of the hold out samples 

fitted.normix = predict.normmix(X = breast[hold.n,3:32],covar = covar.array.store, mean = mean.mat.store, K=K)
dim(fitted.normix)
predicted.normix = round(apply(fitted.normix,2,mean))
1-sum(as.numeric(breast[,2])[hold.n]==predicted.normix)/56 #// 0.642 // .696 // .661


class1 = breast[as.numeric(breast[,2])==1,3:32]
class2 = breast[as.numeric(breast[,2])==2,3:32]
sum(svd(class1)$d/max(svd(class1)$d))
sum(svd(class2)$d/max(svd(class2)$d))
