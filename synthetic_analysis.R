library(rstiefel)
library(scatterplot3d)
library(MASS)
library(mvtnorm)
library(MCMCpack)
library(rgl)
source('embedding_functions.R')
#generate data
n=500
K=2
m=3
isotropic.var = 2
source('generate_toy_data.R')
plot3d(X,type='s',size=.5,col=rep(c('black','red'),each=n),axes=F,xlab='',ylab='',zlab='')
Synthetic.data.sets[[4]][[2]]=X
#Estimate mixture components
#store the estimates for .1,.5,1,2,5
#Synthetic.data.sets=vector("list",5)
#for(i in 1:5){
#  Synthetic.data.sets[[i]]=vector("list",10)
#}
#scatterplot3d(X,pch=20,color=rep(c(1,2),each=n),xlab="",ylab="",zlab="")
#for(i in 1:10){
#  Synthetic.data.sets[[1]][[i]] = Synthetic.data.sets[[1]][[i]][c(1:500,1001:1500),]
#}
k.means.results=vector("list",5)
for(i in 1:5){
  k.means.results[[i]]=vector("list",10)
}
hold = c(451:500,951:1000)
train = c(1:450,501:950)
#Find K-means estimates
for(variance.size in 4:4){
  for(trial.number in 1:10){
    X = Synthetic.data.sets[[variance.size]][[trial.number]]
    k.means.fit = kmeans(X[train,],2)
    cluster = rep(0,100)
    for(i in 1:100){
        cluster[i] = which.min(c(sum((X[hold[i],]-k.means.fit$centers[1,])^2),sum((X[hold[i],]-k.means.fit$centers[2,])^2)))
      }
    k.means.results[[variance.size]][[trial.number]] = sum(c(cluster[1:50]==1,cluster[51:100]==2))/100
  }
}
for(i in 1:5){
  for(j in 1:10){
    temp =k.means.results[[i]][[j]]
    if(temp<.5){
      k.means.results[[i]][[j]] = 1-temp
    }
  }
}
true.cluss = rep(c(1,2),each=50)
#Find mixture of normals estimates
source('mixture_normal_functions.R')
normal.results = vector("list",5)
for( i in 1:5){
  normal.results[[i]] = vector("list",10)
}
for(variance.size in 5:5){
  for(trial.number in 1:10){
    print(paste("variance size is ",variance.size," trial number is ",trial.number))
    X = Synthetic.data.sets[[variance.size]][[trial.number]]
    HoldX = X[hold,]
    X = X[train,]
    source('normal_mixture.R')
    fitted.normix = predict.normmix(X = HoldX, covar = covar.array.store, mean = mean.mat.store, K=K)
    predicted.normix = round(apply(fitted.normix,2,mean))
    normal.results[[variance.size]][[trial.number]]=1-sum(rep(c(1,2),each=50)==predicted.normix)/100
  }
}
for(variance.size in 1:1){
  for(trial.number in 9:10){
    print(paste("variance size is ",variance.size," trial number is ",trial.number))
    X = Synthetic.data.sets[[variance.size]][[trial.number]]
    HoldX = X[hold,]
    X = X[train,]
    source('normal_mixture.R')
    fitted.normix = predict.normmix(X = HoldX, covar = covar.array.store, mean = mean.mat.store, K=K)
    predicted.normix = round(apply(fitted.normix,2,mean))
    normal.results[[variance.size]][[trial.number]]=1-sum(rep(c(1,2),each=50)==predicted.normix)/100
  }
}

predictions = vector("list",1)
for(variance.size in 1:1){
  prediction.tracker=0
  for(trial.number in c(10)){
    prediction.tracker=prediction.tracker+1
    print(paste("variance size is ",variance.size," trial number is ",trial.number))
    X = Synthetic.data.sets[[variance.size]][[trial.number]]
    #HoldX = X[hold,]
    #X = X[train,]
    source('normal_mixture.R')
    fitted.normix = predict.normmix(X = X, covar = covar.array.store, mean = mean.mat.store, K=K)
    predictions[[prediction.tracker]] = round(apply(fitted.normix,2,mean))
    
  }
}




summary(as.numeric(normal.results[[2]]))
summary(as.numeric(k.means.results[[3]]))



iter = 10000
#talk and report
talk = T
report = 100


#Set initial temperature for gibbs posterior. Is this a tuning parameter? Too expensive to set through cross-validation
temperature  = 10^-6

#metropolis hastings (initial) tuning
sub.met.sd = rep(.8,K)

#Pick initial guess of dimensions for each subspace (shouldnt really matter, but being close helps)
#needs to be of length K
R = rep(1,K)

#Find, through metropolis-hastings, the gibbs posterior for the subspaces U and the affine
#translations theta for each component. The loss function is distance to the union of subspaces,
# the random walk is a spherized normal walk. Note: Assumes at this point the data is called 'X'
source('gibbs_posterior_metrop.R')
accept
#acceptance 3882 4842
#Get the clustering for each subspace pair
#Find classification errors
clus.mat = matrix(0,nrow=iter,ncol=n)
for(i in 1:iter){
  if(talk == T){
    if(i%%report==0){
      print(paste('iteration ', i, 'at ', Sys.time()))
      
    }
  }
  projections = P.mat[[i]]
  theta.mus = theta.mat[[i]]
  clus.this = fast.log.loss(x=X,P.list = projections, mu = theta.mus,temperature=temperature)$clus
  clus.mat[i,]=clus.this 
}
mean.clus.mat = apply(clus.mat[5001:10000,],2,mean)
sum(mean.clus.mat[1:500]<=1.5)/(n/2) #Gives the percent correctly predicted in
sum(mean.clus.mat[501:1000]>=1.5)/(n/2)

rounded.mean.clus.mat = round(mean.clus.mat)
samp.wrong = sample(size=10,x=which(mean.clus.mat[501:1000]<1.5))
rounded.mean.clus.mat[500+samp.wrong]=2

sum(rounded.mean.clus.mat[501:1000]>1.5)/1000 
sum(mean.clus.mat[501:1000]>1.5)/1000+500/1000



scatterplot3d(rbind(X1,X2),pch=20,color=round(mean.clus.mat),xlab="",ylab="",zlab="")

#look at k means
kmeans.fit = kmeans(X,2)
scatterplot3d(rbind(X1,X2),pch=20,color=kmeans.fit$cluster,xlab="",ylab="",zlab="")
sum(kmeans.fit$cluster[1:500]>1.5)/1000 + sum(kmeans.fit$cluster[501:1000]<1.5)/1000
#Synth.X = X
X=Synth.X
#save(list="Synth.X",file="~/Dropbox/exterior paper and code/Code/Synth_X.Rdata")
#save(list="Synthetic.data.sets",file='~/Dropbox/exterior paper and code/Code/Synthetic_data_sets.Rdata')

plot3d(Synthetic.data.sets[[1]][[10]],type='s',size=.5,col=rep(c('black','red'),each=n),axes=F,xlab='',ylab='',zlab='')
plot3d(Synthetic.data.sets[[1]][[10]],type='s',size=.5,col=predictions[[4]],axes=F,xlab='',ylab='',zlab='')

#3 and 6 look good, try 8 as well, and 10
plot3d(X,type='s',size=.5,col=rep(c('black','red'),each=n),axes=F,xlab='',ylab='',zlab='')

fitted.submix = predict.submix(P = P.mat, theta = theta.mat, X = X,K=K)
submix.prediction =apply(apply(fitted.submix$cluster, 2, table),2,which.max)

plot3d(Synthetic.data.sets[[1]][[10]],type='s',size=.5,col=rep(c(2,1),each=500),axes=F,xlab='',ylab='',zlab='')#truth
plot3d(Synthetic.data.sets[[1]][[10]],type='s',size=.5,col=kmeans(X,2)$cluster,axes=F,xlab='',ylab='',zlab='')#k-means
plot3d(Synthetic.data.sets[[1]][[10]],type='s',size=.5,col=round(apply(fitted.normix,2,mean)),axes=F,xlab='',ylab='',zlab='')#normal
plot3d(Synthetic.data.sets[[1]][[10]],type='s',size=.5,col=normests,axes=F,xlab='',ylab='',zlab='')#subspace
plot3d(Synthetic.data.sets[[1]][[10]][c(1:500,disag[-(1:24)]),],type='s',size=.5,col=c(rep(1,500),rep(2,(length(disag)-24))),axes=T,xlab='x',ylab='y',zlab='z')#normal
#may also be able to lose the lat clus command since it should be the same thing
#disag 9,12,23,24
normests =  round(apply(fitted.normix,2,mean))

