source('embedding_functions.R')

#load in data set,  the data is in the file heart.dat. 
#The first 13 are attributes from patient. THe last column indicates 
# 1 for absence of heart disease, and 2 for presence
heart = read.table(file='heart.dat',header=F)
n=dim(heart)[1]
#the attributes are described as
#Real: 1,4,5,8,10,12
#Ordered:11,
#Binary: 2,6,9
#Nominal:7,3,13
X = heart[,1:13]
for(i in c(1,4,5,8,10,12)){
  X[,i]=(X[,i]-mean(X[,i]))/sd(X[,i]) 
}
hold.n = sample(1:n,floor(n/10),replace=F)
#for comparison to the truth
true.cluss = heart[,14]-1;true.cluss = true.cluss[-hold.n]; true.hold.cluss = true.cluss[hold.n]
X = X[-hold.n,]
n = n - length(hold.n)
total.n = n + length(hold.n)
theta.true.list = vector("list",2)
theta.true.list[[1]] = apply(X[true.cluss==0,],2,mean)
theta.true.list[[2]] = apply(X[true.cluss==1,],2,mean)
#ambient dimension
m=dim(X)[2]
#number of components to fit - we know before hand that there are 4 types of vehicle
K =2
#number of observations 
n=dim(X)[1]



#Pick number of iterations
iter = 10000
#talk and report
talk = T
report = 100


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
source('embedding_functions.R')
source('mixture_normal_functions.R')
fitted.submix = predict.submix(P = P.mat, theta = theta.mat, X = heart[hold.n,1:13],K=K)
submix.prediction =apply(apply(fitted.submix$cluster, 2, table),2,which.max)
1-sum(1+as.numeric(heart[,14])[hold.n]==submix.prediction)/length(hold.n) #.63 /.778

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
sum(!clus.est==true.cluss)/n #.911 / .906  / .909 - !clus.est because the components aren't necessarily ordered the same way
fit1 = glm(true.cluss ~ as.matrix(X), family="binomial") 
summary(fit1)
sum(round(fitted.values(fit1))==true.cluss)/n #.855 / .851 /.852
sum(round(predict(object = fit1, newdata = as.matrix(heart[hold.n,1:13]),type="response" ))==heart[,14][hold.n])/length(hold.n)


save.image(file='heart_analysis.Rdata')
load(file = "heart_analysis.Rdata")

#Calculate a non-subspace assumption normal mixture model
source('mixture_normal_functions.R')
source('normal_mixture.R')
#Objects are stored in lat.clus.mat, mean.mat.stor, covar.array.store
save.image(file='heart_analysis.Rdata')
#load(file = "heart_analysis.Rdata")
lat.clus.mean = apply(lat.clus.mat,2, mean)
lat.clus.est = round(lat.clus.mean)-1
sum(!lat.clus.est==true.cluss)/n # .560/.58
 

fitted.normix = predict.normmix(X = heart[hold.n,1:13],covar = covar.array.store, mean = mean.mat.store, K=K)
predicted.normix = apply(round(apply(fitted.normix,2,table)),2,which.max) 
1-sum(as.numeric(heart[,14])[hold.n]==predicted.normix)/length(hold.n) #.52/.51




class1 = heart[as.numeric(heart[,14])==1,1:13]
class2 = heart[as.numeric(heart[,14])==2,1:13]
sum(svd(class1)$d/max(svd(class1)$d))
sum(svd(class2)$d/max(svd(class2)$d))

