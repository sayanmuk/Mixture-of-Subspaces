source('embedding_functions.R')

#load in entire data set,  the data is in the xa".".dat files but was compiled into one file
#that I named vehicle.dat. The first 18 are image attributes. THe last column gives the vehicle class
vehicle = read.table(file='vehicle.dat',header=F)
dim(vehicle)
#Data for fitting
X = vehicle[,1:18]
for(i in 1:18){
 X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i]) 
}

#ambient dimension
m=dim(X)[2]
#number of components to fit - we know before hand that there are 4 types of vehicle
K =4
#number of observations 
n=dim(X)[1]
hold.n = sample(1:n,floor(n/10),replace=F)


#true clustering for comparision
true.cluss = as.numeric(vehicle[,19]);true.cluss = true.cluss[-hold.n]; true.hold.cluss = true.cluss[hold.n]
X = X[-hold.n,]
n = n - length(hold.n)
total.n = n + length(hold.n)



#Pick number of iterations
iter = 10000
#talk and report
talk = T
report = 100


#Set initial temperature for gibbs posterior. Is this a tuning parameter? Too expensive to set through cross-validation
temperature  = 10^-8

#metropolis hastings (initial) tuning
sub.met.sd = rep(1,K)

#Pick initial guess of dimensions for each subspace (shouldnt really matter, but being close helps)
#needs to be of length K
R = rep(3,K)

#Find, through metropolis-hastings, the gibbs posterior for the subspaces U and the affine
#translations theta for each component. The loss function is distance to the union of subspaces,
# the random walk is a spherized normal walk. Note: Assumes at this point the data is called 'X'
source('embedding_functions.R')
source('gibbs_posterior_metrop.R')
fitted.submix = predict.submix(P = P.mat, theta = theta.mat, X = vehicle[hold.n,1:18],K=K)
submix.prediction =apply(apply(fitted.submix$cluster, 2, table),2,which.max)
1-sum(1+as.numeric(vehicle[,19])[hold.n]==submix.prediction)/length(hold.n) 


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

clus.est=apply(clus.mat,2,Mode)
sum(clus.est==true.cluss)/n #.480 //.475
library('glmnet')
fit1=glmnet(as.matrix(X),clus.est,family="multinomial")
apply(predict(fit1,as.matrix(X),type="response",s=0.01)[,,1],1,which.max)
sum(apply(predict(fit1,as.matrix(X),type="response",s=0.01)[,,1],1,which.max)==true.cluss)/n #.460 //.471
save.image(file='vehicle_analysis.Rdata')
#load(file="vehicle_analysis.Rdata")
#Calculate a non-subspace assumption normal mixture model
source('mixture_normal_functions.R')
source('normal_mixture.R')
#Objects are stored in lat.clus.mat, mean.mat.stor, covar.array.store
save.image(file='vehicle_analysis.Rdata')
#load(file="vehicle_analysis.Rdata")
lat.clus.mean = apply(lat.clus.mat,2, mean)
lat.clus.est = round(lat.clus.mean)-1
sum(!lat.clus.est==true.cluss)/n 

fitted.normix = predict.normmix(X = vehicle[hold.n,1:18],covar = covar.array.store, mean = mean.mat.store, K=K)
predicted.normix = apply(round(apply(fitted.normix,2,table)),2,which.max) 
1-sum(as.numeric(vehicle[,19])[hold.n]==predicted.normix)/length(hold.n)


class1 = vehicle[as.numeric(vehicle[,19])==1,1:18]
class2 = vehicle[as.numeric(vehicle[,19])==2,1:18]
class3 = vehicle[as.numeric(vehicle[,19])==3,1:18]
class4 = vehicle[as.numeric(vehicle[,19])==4,1:18]
sum(svd(class1)$d/max(svd(class1)$d))
sum(svd(class2)$d/max(svd(class2)$d))
sum(svd(class3)$d/max(svd(class2)$d))
sum(svd(class4)$d/max(svd(class2)$d))
