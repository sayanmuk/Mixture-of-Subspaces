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

con.distance.het = function(U,V,subspace=T){
  #Find distance between m-dimensional subspaces
  #In embedded space
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

con.distance.hom = function(U,V,subspace=T){
  #Find distance between ambient m-dimensional subspaces using the chordal
  #distance (conway distance for same dimension subspaces)
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

ConDist = function(U,V, subspace = T){
  #Find conway distance between two subspaces
  #or projection distance if they are the same
  #dimension
  m=dim(U)[1];
  d = sqrt(m)
  if(subspace==F){
    temp1 = eigen(U);
    U = temp1$vectors[,temp1$values>10^(-12)]
    temp2 = eigen(V);
    V = temp2$vectors[,temp2$values>10^(-12)]
  }
  if(dim(U)[2]==dim(V)[2]){
    ret = con.distance.hom(U,V, subspace = T)
  }else{
    ret = con.distance.het(U,V, subspace = T)
  }
  
  return(ret)
  
}

AsimovDist = function(U,V){
  #Computes the asimov distance between subspaces
  k = dim(U)[2]
  return(acos(svd(t(U)%*%V)$d[k])) 
}

BinetCauchyDist = function(U,V){
  #Binet-Cauchy distance 
  temp = svd(t(U)%*%V)
  princ = acos(temp$d)
  k = dim(U)[2]
  prod = 1
  for(i in 1:k){
    prod = prod*cos(princ[i])^2
  }
  return(1 - prod)
}

ChordalDist = function(U, V){
  #Chordal distance
  temp = svd(t(U)%*%V)
  princ = acos(temp$d)
  return(sqrt(sum(sin(princ)^2)))
}

FubiniStudyDist = function(U,V){
  #Fubini - Study distance
  temp = svd(t(U)%*%V)
  princ = acos(temp$d)
  prod = 1
  for(i in 1:length(princ)){
    prod = prod * cos(princ[i])
  }
  return(acos(prod))
}
 
MartinDist = function(U,V){
  temp = svd(t(U)%*%V)
  princ = acos(temp$d)
  prod = 1
  for(i in 1:length(princ)){
    prod = prod*(1/cos(princ[i])^2)
  }
  return(sqrt(log(prod)))
}

ProcrustesDist = function(U, V){
  temp = svd(t(U)%*%V)
  princ = acos(temp$d)
  return(2*sqrt(sum(sin(princ/2)^2)))
}

SpectralDist = function(U, V){
  temp = svd(t(U)%*%V)
  k = dim(U)[2]
  princ = acos(temp$d)
  return(2*sin(princ[k]/2)) 
}


