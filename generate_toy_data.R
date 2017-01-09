#Toy example, of a line intersecting a plane in 3-space

k.true=sort(sample(1:(m-1),size = K, replace=T),decreasing=F)


#J.true=2
#Generate the subspaces
Utrue.list = list()
for(k in 1:K){
  Utrue.list[[k]] = rustiefel(m=m,R=k.true[k])
}
NU.true = list()
for(k in 1:K){
  NU.true[[k]]=Null(Utrue.list[[k]])
}
PNU.list = list()
for(k in 1:K){
  PNU.list[[k]] = NU.true[[k]]%*%t(NU.true[[k]])
}
#generate means in subspace coordinates
mutrue.list = list()
for(k in 1:K){
  mutrue.list[[k]] = rnorm(k.true[k])
}

#generate the residual space noise level
phitrue.list = rep(10,K)
sigmatrue.list = rep(isotropic.var,K)
#generate the subspace variances
sigma0.true.list = list()
for(k in 1:K){
  sigma0.true.list[[k]]=runif(k.true[k],isotropic.var,5.1)
}
Sigma0.true.list = list()
for(k in 1:K){
  Sigma0.true.list[[k]]=diag(sigma0.true.list[[k]],k.true[k])
}


#generate the euclidean space coordinate mean vector theta
theta.true.list = list()
for(k in 1:K){
  theta.true.list[[k]] = PNU.list[[k]]%*%rnorm(m)
}
X = matrix(0,nrow=n*k,ncol = m)
for(k in 1:K){
  X[(n*(k-1)+1):(n*k),]=mvrnorm(n,mu=Utrue.list[[k]]%*%mutrue.list[[k]]+theta.true.list[[k]],
                                Sigma = sigmatrue.list[k]^2*diag(m)/10+Utrue.list[[k]]%*%(Sigma0.true.list[[k]]-sigmatrue.list[k]*diag(k.true[k]))%*%t(Utrue.list[[k]]))
}
 
#generate the data
#X1 = mvrnorm(n/2,mu=U1.true%*%mu1.true+theta1.true,Sigma=sigma1.true^2*diag(3)/10+U1.true%*%(Sigma01.true-sigma1.true*diag(k.true[1]))%*%t(U1.true))
#X2 = mvrnorm(n/2,mu=U2.true%*%mu2.true+theta2.true,Sigma=sigma2.true^2*diag(3)/10+U2.true%*%(Sigma02.true-sigma2.true*diag(k.true[2]))%*%t(U2.true))
#true.cluss = c(rep(1,500),rep(2,500))
#X=rbind(X1,X2)
#n=dim(X)[1]
