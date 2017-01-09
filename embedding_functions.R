#setwd("~/Dropbox/exterior paper and code/Code")
library(rstiefel)
library(scatterplot3d)
library(RcppArmadillo)
library(inline)
library(MASS)

Trace = function(X){
  return(sum(diag(X)))
}

embed = function(P, subspace = F){
  #returns the conway embedding of the projection matrix P if subspace =F, 
  #or the subspace P if subspace =T
  m=dim(P)[1]
  if(subspace == T){
    P = P%*%t(P)
  }
  p=c()
  for(i in 1:m){
    p=c(p,P[i:m,i])
  }  
  return(p)
}

distance = function(U,V,subspace=T){
  #Find distance between m-dimensional subspaces
  m=dim(U)[1];
  d = sqrt(m)
  if(subspace==T){
  P1 = U%*%t(U);
  P2 = V%*%t(V);
  }else{
    P1=U;
    P2=V;
  }
  PZ = diag(.5,m)
  p1 = embed(P1)
  p2 = embed(P2)
  pz = embed(PZ)
  angle = acos(sum((p1-pz)*(p2-pz))/(sqrt(sum((p1-pz)^2))*sqrt(sum((p2-pz)^2))))
  return(d*angle/2)
}

con.distance = function(U,V,subspace=T){
  #Find distance between ambient m-dimensional subspaces
  #using the conway distance
  #aka projection distance
  m=dim(U)[1];
  d = sqrt(m)
  if(subspace==F){
    temp1 = eigen(U);
    U = temp1$vectors[,temp1$values>10^(-12)]
    temp2 = eigen(V);
    V = temp2$vectors[,temp2$values>10^(-12)]
  }
  sigma = svd(t(U)%*%V)$d
  sigma[sigma>1]=1
  princ = acos(sigma)
  return(sqrt(sum((sin(princ))^2)))
  
}



unembed = function(p){
  #Given an embedding of a subspace p,
  #return the Projection matrix P
  m = (sqrt(8*length(p)+1)-1)/2
  P = matrix(0,nrow=m,ncol=m);
  take = m:1
  place = 1
  for(i in 1:m){
    P[i:m,i]=p[place:(place+take[i]-1)]
    place=place+take[i]
  }
  P = P + t(P)-diag(diag(P),m)
  return(P)
}

con2sub = function(P,d=Inf,return.proj = T){
  #go from the embededd space to either the 
  #projection matrix or the subspace depending
  #on return.proj=T
 
  #get dimensions
  m = dim(P)[1]
  
  if(d == Inf){
    
    prj.rank = sum(diag(P))
    
    if(prj.rank%%1>.5){
      d=ceiling(prj.rank)
    }else{
      d=floor(prj.rank)
    }
    
  }
  temp = eigen(P)
  U = temp$vectors[,1:d]
  if(return.proj ==T){
    return(U%*%t(U))
  }else{
    return(U)
  }
}

runif.sphere = function(p=3){
  #sample a p-dimensional vector from the unit sphere in R^p
  v = rnorm(p,0,1)
  v=v/sqrt(sum(v^2))
  return(v)
}

runif.conway = function(n=1,m=3){
  
  p=m*(m+1)/2
  r = sqrt(m)/2
  PZ = diag(rep(.5,m),m)
  pz = embed(PZ)
  v =matrix(0,nrow=n,ncol=p)
  for(i in 1:n){
    v[i,]=r*runif.sphere(p)+pz
  }
  if(n==1){
  return(as.numeric(v))
  }else{
    return(v)
  }
}

conway.spherize = function(z,m=3){
  #Renormalize the point Z onto the conway sphere
  #by first mapping it to the unit sphere, and then 
  #translating and scaling to the Conway sphere
  z = z/sqrt(sum(z^2))
  pz = embed(diag(rep(.5,m),m))
  r = sqrt(m)/2
  return(z*r+pz)
}


embed.rustiefel = function(n,m,R){
  X = matrix(0,nrow=n,ncol=m*(m+1)/2)
  for(i in 1:n){
    X[i,]= embed(rustiefel(m=m,R=R),subspace=T)
  }
  return(X)
}


normal.walk = function(n=1000,m=3){
  v = matrix(0,nrow=n,ncol=m)
  v[1,]=rnorm(m)
  for( i in 2:n){
    v[i,]=rnorm(m)+v[i-1,]
  }
  return(v)
}

sphere.walk = function(n=1000,m=3){
  #generates a random walk on the unit sphere of length n
  v = matrix(0,nrow=n,ncol=m)
  s = matrix(0,nrow=n,ncol=m)
  v[1,]=rnorm(m)
  s[1,]=v[1,]/sum(v[1,]^2)
  for(i in 2:n){
    v[i,]=rnorm(m)+v[i-1,]
    s[i,] = v[i,]/sum(v[i,]^2)
  }
  return(list('v'=v,'s'=s))
}

sphere.step = function(z,sd=1,m=3){
  #given a specific location on the latent walk,
  #talk a step and spherize it. Return the new
  #location
  p = m*(m+1)/2
  v = rnorm(p,mean=z,sd=sd)
  s = v/sqrt(sum(v^2))
  return(list('z'=v,'s'=s))
  
  
}

conway.sphere.step = function(z, sd=1,m=3){
  #Given a specific location of the latent walk
  #take a step and spherize it. Return both the new location
  #of the latent walk and the spherized step. Now the 
  #'s' is given on the Conway sphere in m(m+1)/2 space
  r= sqrt(m)/2
  PZ = diag(rep(.5,m),m)
  pz = embed(PZ)
  step = sphere.step(z,sd=sd,m=m)
  s=r*step$s+pz
  return(list('s'=s,'z'=step$z))
}

conway.step = function(z,sd=1,m=3){
  #Given a location on the Conway sphere z,
  #a standard deviation sd, and the ambient dimension m
  # take a random step on the sphere, and return the sphere location
  r = sqrt(m)/2;
  p = m*(m+1)/2
  PZ = diag(rep(.5,m),m);
  pz= embed(PZ)
  tangent.step = rnorm(p,z,sd)
  tangent.step = tangent.step/sqrt(sum(tangent.step^2))
  step = tangent.step*r+pz
  return(list('s'=step,'z'=tangent.step))
}



gibbs.loss.prj=function(x , P, mu = NULL, subspace = F){
  #P is one a  m x m projection matrices (or a subspace)
  # x is the set of n observations, so is n x m. Finds the distance
  #between theclosest point on the subspace P and x
  n = dim(x)[1]
  if(subspace==T){
    P=P%*%P^T
  }
  m = dim(P)[1]
  d = sum(diag(P))
  
  if(is.null(mu)){
    mu = rep(0,m)
  }
  
  #Find the distance between the projection of the point and the point itself for each 
  #observation
  dist = rep(0, n)
  for(i in 1:n){
    dist[i] = d*sqrt(sum(m^2*(P%*%(t(x[i,])-mu) - (t(x[i,])-mu))^2))
  }
  
  return(dist)
  
}

log.gibbs.loss = function(x, P.list, mu = Inf,temperature=1, subspace =F){
  
  #P.list is a list of projection matrices, indexed by P.list[[k]]
  #mu is the list of theta, where each can be accessed as mu[[k]]
  #given a set of observations x
  K = length(P.list)
  n = dim(x)[1]
  dist = matrix(0,nrow = K, ncol = n)
  for(k in 1:K){
    #print('log gibbs loss')
   # print(k)
    dist[k,] = gibbs.loss.prj(x, P.list[[k]],mu = mu[[k]],subspace=subspace)
  }
  mindist = apply(dist,2,min)
  whichmindist = apply(dist,2,which.min)
  loss = -1*n*temperature *sum(mindist)
  return(list('loss'=loss,'clus'=whichmindist))
}


src = '

using namespace arma;

int n = as<int>(nS);
int d = as<int>(dS);
int m = as<int>(mS);
mat P = as<mat>(PS);
mat x = as<mat>(xS);
mat mu = as<mat>(muS);
mat dist = zeros(1,n);
int i;
for(i=0;i<n;i++){
dist(0,i) = d*pow(sum(pow(m,2)*pow(P*(trans(x.row(i))-mu) - (trans(x.row(i))-mu),2)),.5);
}
return wrap(dist);
'




fast.loss.prj = cxxfunction(signature(nS = "integer", dS = "integer",mS = "integer",
                                      PS = "double", muS="double", xS="double"),
                            plugin = "RcppArmadillo",
                            body = src)


fast.log.loss = function(x, P.list, mu, temperature = 1, subspace =F){
  #P.list is a list of projection matrices, indexed by P.list[[k]]
  #mu is the list of theta, where each can be accessed as mu[[k]]
  #given a set of observations x. Uses the Rcpp function fast.loss.prj
  #to iterate over the observations to get the loss
  K = length(P.list)
  n = dim(x)[1]
  m = dim(x)[2]
  x = as.matrix(x)
  dist = matrix(0,nrow = K, ncol = n)
  for(k in 1:K){
    dist[k,] = fast.loss.prj(xS = x, PS = P.list[[k]],muS = mu[[k]],nS = n, dS =dim(P.list[[k]])[2],mS =m )
  }
  mindist = apply(dist,2,min)
  whichmindist = apply(dist,2,which.min)
  loss = -1*n*temperature *sum(mindist)
  return(list('loss'=loss,'clus'=whichmindist))
  

}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

 


