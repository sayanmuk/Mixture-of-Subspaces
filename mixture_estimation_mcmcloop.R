#########################################START MCMC #########
for(iter in 1:MCMCiter){
  #Begin MCMC with the initialized values
  if(iter%%10==0){
    print(paste('iteration ',iter,' at ',Sys.time()))
  }
  #print('Start MCMC')
  #First update the U given the clusters for each thing
  #update Summary statistics needed to sample the BMF distributions
  #number of observations
  X1 = X[which(clus==1),]
  n1 = dim(X1)[1]
  X2 = X[which(clus==2),]
  n2 = dim(X2)[1]
  #Calculate sample means for each group
  x1.mean = apply(X1,2,mean)
  x2.mean = apply(X2,2,mean)
  x1.sum = apply(X1,2,sum)
  x2.sum = apply(X2,2,sum)
  
  
  #F1, convention is the first integer corresponds to which summary
  #it is and the second corresponds to the cluster
  F11 = matrix(0,3,k[1])
  F31 = matrix (0,3,3)
  ssqx1 = 0
  #F3, the data driven part
  for(i in 1:n1){
    F11 = F11 + X1[i,]%*%t(mu1)
    F31 = F31 + X1[i,]%*%t(X1[i,])
    ssqx1 = ssqx1 + t(X1[i,])%*%X1[i,]
  }
  F11=F11%*%(1/Sigma01)
  
  F12 = matrix(0,3,k[2])
  F32 = matrix(0, 3, 3)
  ssqx2 = 0;
  for( i in 1:n2){
    F12 = F12 + X2[i,]%*%t(mu2)
    F32 = F32 + X2[i,]%*%t(X2[i,])
    ssqx2 = ssqx2 + t(X2[i,])%*%X2[i,]
  }
  F12 = F12%*%diag(1/sigma02) 
  
  #F2, the variance part
  
  F21 = (1/2)*((1/sigma1)*diag(k[1]) - diag(1/sigma01,nrow=k[1]))
  F22 = (1/2)*((1/sigma2)*diag(k[2]) - diag(1/sigma02, nrow=k[2]))
  
  #Now the F statistics are collected, we need the null space of each associated theta 
  Nt1 = Null(theta1)
  Nt2 = Null(theta2)
  Ub1 = t(Nt1)%*%U1
  Ub2 = t(Nt2)%*%U2
  # Nt1%*% Ub1%*%t(theta1)
  #t(Nt1%*%t(Nt1)%*%U1)%*%theta1
  #print(paste('Sampling U1 at ',Sys.time()))
  Ub1 = rbmf.vector.gibbs(A=t(Nt1)%*%F31%*%Nt1,c=t(Nt1)%*%F11,x= Ub1)
  #B=F21,
  #C=t(Nt1)%*%F11
  U1 = Nt1%*%Ub1
  #print(paste('Done Sampling U1 at ',Sys.time()))
  #check 0's
  #t(U1)%*%theta1
  
  #print(paste('Sampling U2 at ',Sys.time()))
  count = 0
  for(j in sample(1,1:2,replace=F)){
    count=count+1
    Ub2[,j] = rbmf.vector.gibbs(A=t(Nt2)%*%F32%*%Nt2, c=t(Nt2)%*%F12[,j],  Ub2[,j])
    notj = setdiff(c(1,2),j)
    Ub2[,notj]=Null(Ub2[,j])
  }
  #,B=F22
  Ub2 = rbing.matrix.gibbs(A=t(Nt2)%*%F32%*%Nt2, B=F22,  Ub2)
  C=t(Nt2)%*%F12
  U2 = Nt2%*%Ub2
  #check 0's
  #t(U2)%*%theta2
  #print(paste('Done Sampling U2 at ',Sys.time()))
  #Now the U's are updated (add C^TX part - check rstief), now we want to update theta
  #Update the theta's
  M1 = Null(U1)
  M2 = Null(U2)
  
  #Given updated U's, get the new covariances
  Sigma1 = U1%*%(Sigma01-sigma1*diag(k[1]))%*%t(U1)+sigma1*diag(3)
  Sigma2 = U2%*%(Sigma02 - sigma2*diag(k[2]))%*%t(U2) + sigma2*diag(3)
  
  #And find their inverses since we will need to solve fr it multiple times
  Sigma1.inv = solve(Sigma1)
  Sigma2.inv = solve(Sigma2)
  
  #Find the Posterior covariances and means
  S.t1 = solve(n1*Sigma1.inv+diag(1/diag(S.theta1),3))
  S.t2 = solve(n2*Sigma2.inv+diag(1/diag(S.theta2),3))
  m.t1 = S.t1%*%(n1*Sigma1.inv%*%x1.mean+diag(1/diag(S.theta1),3)%*%m.theta1)
  m.t2 = S.t2%*%(n2*Sigma1.inv%*%x2.mean+diag(1/diag(S.theta2),3)%*%m.theta2)
  
  #Find the mean and covariance projected into the null space of U
  Sbt1 = solve(t(M1)%*%solve(S.t1,M1))
  Sbt2 = solve(t(M2)%*%solve(S.t2,M2))
  mbt1 = Sbt1%*%t(M1)%*%solve(S.t1,m.t1)
  mbt2 = Sbt2%*%t(M2)%*%solve(S.t2,m.t2)
  
  #Sample from the m-k dimensional normal
  thetab1 = mvrnorm(1,mbt1,Sbt1)
  thetab2 = mvrnorm(1,mbt2,Sbt2)
  
  #Project back into the normal coordinates
  theta1 = M1%*%thetab1
  theta2 = M2%*%thetab2
  
  #checking that they are still orthogonal
  #t(U1)%*%theta1
  #t(U2)%*%theta2
  
  
  #End Sampling Theta, begin sampling mu
  #Sample from a normal distribution for each cluster
  Sigma01.inv = diag(1/diag(Sigma01),k[1]);
  Sigma02.inv = diag(1/diag(Sigma02),k[2]);
  Sig.mu1 = solve(n1*Sigma01.inv+S.mu1inv)
  Sig.mu2 = solve(n2*Sigma02.inv+S.mu2inv)
  
  m.mu1.post = Sig.mu1%*%(t(U1)%*%Sigma1.inv%*%x1.sum + S.mu1inv%*%m.mu1)
  m.mu2.post = Sig.mu2%*%(t(U2)%*%Sigma2.inv%*%x2.sum + S.mu2inv%*%m.mu2)
  
  mu1 = mvrnorm(n=1,mu=m.mu1.post,Sigma = Sig.mu1);
  mu2 = mvrnorm(n=1, mu=m.mu2.post, Sigma= Sig.mu2);
  
  #End sampling mu, begin sampling the idiosyncratic noise sigma
  sum.id1=0
  sum.id2=0
  sum.sub1=rep(0,k[1])
  sum.sub2=rep(0,k[2])
  
  for(i in 1:n1){
    sum.id1 = sum.id1 +t(X1[i,])%*%U1%*%t(U1)%*%X1[i,]
    sum.sub1 = sum.sub1 + (t(U1)%*%X1[i,]-mu1)^2
  }
  
  for(i in 1:n2){
    sum.id2 = sum.id2 + t(X2[i,])%*%U2%*%t(U2)%*%X2[i,]
    sum.sub2 = sum.sub2 + (t(U2)%*%X2[i,]-mu2)^2
  }
  
  a.p.s1 = 1/2*(n1*(3-k[1]))+a.s1
  a.p.s2 = 1/2*(n2*(3-k[2]))+a.s2
  b.p.s1 = b.s1 + 1/2*ssqx1 + n1/2*sum(theta1^2)-1/2*sum.id1 - t(theta1)%*%x1.sum
  b.p.s2 = b.s2 + 1/2*ssqx2 + n2/2*sum(theta2^2)-1/2*sum.id2 - t(theta2)%*%x2.sum
  
  sigma1 = 1/rgamma(1,a.p.s1,b.p.s1)
  sigma2 = 1/rgamma(1,a.p.s2,b.p.s2)
  
  #subspace variances
  a.post.u1 = n1/2+a.u1
  a.post.u2 = n2/2+a.u2
  b.post.u1 = b.u1 + 1/2*sum.sub1
  b.post.u2 = b.u2 + 1/2*sum.sub2
  
  
  sigmo01 = 1/rgamma(k[1],a.post.u1,b.post.u1);
  sigma02 = 1/rgamma(k[2],a.post.u2,b.post.u2)
  
  #truncated distribution
  sigma01[sigma01<sigma1]=(runif(k[1],1,3)*sigma1)[sigma01<sigma1]
  sigma02[sigma02<sigma2]=(runif(k[2],1,3)*sigma2)[sigma02<sigma2]
  Sigma01 = diag(sigma01,k[1])
  Sigma02 = diag(sigma02,k[2])
  Sigma1 = U1%*%(Sigma01-sigma1*diag(k[1]))%*%t(U1)+sigma1*diag(3)
  Sigma2 = U2%*%(Sigma02 - sigma2*diag(k[2]))%*%t(U2) + sigma2*diag(3)
  #Sample the weights of pi
  pi.clus = as.numeric(rdirichlet(1,alpha+c(n1,n2)))
  #calculate the matrix of q_ij proportionality of likelihoods
  
  #Get the means
  clus1.mean = U1%*%mu1+theta1
  clus2.mean = U2%*%mu2+theta2
  
  q.clus[,1]=pi.clus[1]*dmvnorm(X,mean = clus1.mean,sigma=Sigma1)
  q.clus[,2]=pi.clus[2]*dmvnorm(X,mean=clus2.mean, sigma=Sigma2)
  #rescale weights to sum to 1
  q.clus = q.clus/apply(q.clus,1,sum)
  for(i in 1:n){
    clus[i] = which(as.numeric(rmultinom(1,1,prob = q.clus[i,]))==1)
  }
  #Now all the observations have been reassigned
  
  #Save after a burn in period
  if(iter>burn){
    save=iter-burn
    U1.save[,,save] = U1
    U2.save[,,save] = U2
    theta1.save[,save] = theta1
    theta2.save[,save] = theta2
    mu1.save[,save] = mu1
    mu2.save[,save] = mu2
    sigma01.save[,save] = sigma01
    sigma02.save[,save] = sigma02
    sigma1.save[save] = sigma1
    sigma2.save[save] = sigma2
    clus.save[,save] = clus
  }  
  #end MCMC
  
}
