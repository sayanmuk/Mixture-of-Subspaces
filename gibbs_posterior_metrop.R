#Intialize the list of K subspaces
U.list = list()
for(k in 1:K){
  U.list[[k]] = rustiefel(m=m, R = R[k])
}

#Initialize list of projection matrices
P.list = list()
for(k in 1:K){
  P.list[[k]] = U.list[[k]]%*%t(U.list[[k]])
}

#Initialize theta list
kmeans.init = kmeans(X,K)
theta.list = list()
for(k in 1:K){
  theta.list[[k]] =Null(U.list[[k]])%*%t(Null(U.list[[k]]))%*% kmeans.init$centers[k,]
} 

#Initialize theta storage
theta.mat = vector("list",iter)

#Intialize subspace storage. Each subspace iteration will store U.mat[[iter]] = U.list
#Allowing each subspace to be accessed as U.mat[[iter]][[k]]
U.mat = vector("list",iter)
P.mat = vector("list",iter)


#Initialize the distance matrix to be stored every iteration
distance = matrix(0,nrow=K,ncol = n)

#Initialize storage for the latent normal and sphere walks 
z.store = array(0,dim = c(iter,K,m*(m+1)/2))
s.store = array(0,dim = c(iter, K, m*(m+1)/2))

#initialize a random normal location
z = matrix(0,K,m*(m+1)/2)
s = matrix(0,K,m*(m+1)/2)
for(k in 1:K){
  temp = conway.sphere.step(z=rep(0,m*(m+1)/2),sd=1,m=m)
  z[k,]=temp$z
  s[k,]=temp$s
}

#Initialize acceptance counter
accept = rep(0,K)

#Get loss for initial estimates
curr.lossclus = fast.log.loss(x=X, P.list = P.list,mu = theta.list,temperature=temperature)
curr.loss = curr.lossclus$loss
curr.clus = curr.lossclus$clus

#initialize and store the clustering
lat.clus = which(rmultinom(n=n, size = 1, prob=rep(1/K,K))==1,arr.ind=T)[,1]
lat.clus.mat = matrix(0,nrow=iter,ncol=n)
pi.mat = matrix(0,n,K)

#Intialize the multinomial weights
pi.mat = matrix(1/K,nrow=n,ncol=K)
dir.prior = rep(1/K,K)
n.vec = rep(0,K)
r0 = rep(1,K)



#set up tuning parameter
tune = 0 
tune.accept = rep(0,K)
for(i in 1:iter){
 # print(paste('iter is ',i))
  tune = tune+1
  if(talk == T){
    if(i%%report==0){
      print(paste('iteration ', i, 'at ', Sys.time()))
      
    }
  }
  
  
  #For each component, generate a proposal subspace, and then
  #accept or reject it based on the gibbs posterior
  for(k in 1:K){
    #Get proposal projection matrix
    #print('get proposal')
    proposal = conway.step(z=z[k,],sd=sub.met.sd[k],m=m)
    #print('get Subspace')
    prop.sub = con2sub(P=unembed(proposal$s),return.proj = F)
    #restrict samples to m/2 ([m+1]/2 if m is odd ) to stay on lower half of sphere
    #print('Restrict')
    if(is.matrix(prop.sub)){
      if(dim(prop.sub)[2]>ceiling(m/2)){
        if(dim(prop.sub)[2]==m){
          #print('Proposal full dimension, replace with 1 dimension')
          prop.sub= rustiefel(R=1,m=m)
          proposal$z=embed(prop.sub,subspace=T)
        }else{
          #print('Proposal not full dimension')
          prop.sub = Null(prop.sub)
          proposal$z = embed(prop.sub,subspace=T)
        }
      }
    }
   # print('Get projection')
    prop.proj = prop.sub%*%t(prop.sub)
    
    #Set the proposal list
    prop.P.list = P.list
    prop.P.list[[k]] = prop.proj
    
    #Choose an appropriate theta in the null space
    prop.null = Null(prop.sub)
    prop.nullproj = prop.null%*%t(prop.null)
    prop.theta.list = theta.list
    prop.theta.list[[k]] = prop.nullproj%*%theta.list[[k]]
    for(l in 1:K){
      if(is.null(dim(prop.theta.list[[l]]))){
        prop.theta.list[[l]]=matrix(prop.theta.list[[l]],nrow=m,ncol=1)
      }
    }
    
    #Get the loss of the proposal
    prop.lossclus = fast.log.loss(x=X,P.list = prop.P.list, mu = prop.theta.list,temperature=temperature)
    prop.loss = prop.lossclus$loss
    #The coin toss
    toss = log(runif(1,0,1))
    diff.loss = prop.loss - curr.loss
    if(toss<diff.loss){
      accept[k] = accept[k] +1
      tune.accept[k] = tune.accept[k]+1
      #int.acc = int.acc +1
      P.list[[k]] = prop.proj
      U.list[[k]] = prop.sub
      theta.list[[k]] = prop.theta.list[[k]]
      z[k,] =proposal$z 
      curr.loss = prop.loss
      curr.clus = prop.lossclus$clus
      
    }
    
  }
   
  
  #Set new means based on closest clustering
  for(k in 1:K){
    theta.temp = apply(X[curr.clus==k,],2,mean);
    theta.list[[k]]=theta.temp
    theta.list[[k]] = Null(U.list[[k]])%*%t(Null(U.list[[k]]))%*%theta.list[[k]]  
  }
  
  #increase or decrease variance to adjust acceptance rates
  if(tune == 100){
    for(k in 1:K){
      if(tune.accept[k]<10){
        sub.met.sd[k] = .1*sub.met.sd[k]
      }else{
        if(tune.accept[k]<30){
          sub.met.sd[k]=.5*sub.met.sd[k]
        }else{
          if(tune.accept[k]<60){
            sub.met.sd[k]=2*sub.met.sd[k]
          }else{
            if(tune.accept[k]<90){
              sub.met.sd[k] = 5*sub.met.sd[k]
            }else{
              sub.met.sd[k]=10*sub.met.sd[k]
            }
          }
        }
      }
      tune.accept[k]=0
    }
    tune=0
  }
  
   
  #For storage at the end
  P.mat[[i]] = P.list
  U.mat[[i]] = U.list
  theta.mat[[i]]=theta.list
  
  
}
