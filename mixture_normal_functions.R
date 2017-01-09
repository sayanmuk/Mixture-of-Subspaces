#Utility functions to make fitting a normal mixture model faster

src = '

using namespace arma;

int n = as<int>(nS);
int m = as<int>(mS);
mat X = as<mat>(XS);
vec mu = as<vec>(muS);
mat S = zeros(m,m);
int i;
for(i=0;i<n;i++){
S  = S + (trans(X.row(i))-mu)*(X.row(i)-trans(mu)) ;
}
return wrap(S);
'



fastS= cxxfunction(signature(nS = "integer", mS = "integer",
                                      XS = "double", muS="double"),
                            plugin = "RcppArmadillo",
                            body = src)








src = '

using namespace arma;

int n = as<int>(nS);
int m = as<int>(mS);

//Rcpp::Rcout<<"Initialize K " <<std::endl;
int K = as<int>(KS);

//Rcpp::Rcout<<"Initialize X " <<std::endl;
mat X = as<mat>(XS);

//Rcpp::Rcout<<"Initialize mu " <<std::endl;
mat mu = as<mat>(muS);

//Rcpp::Rcout<<"Initialize Sigma " <<std::endl;
Rcpp::NumericVector vecArray(SigmaS);
Rcpp::IntegerVector arrayDims = vecArray.attr("dim");
arma::cube Sigma(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

//Rcpp::Rcout<<"Initialize pi " <<std::endl;
vec pi = as<vec>(piS);

//Rcpp::Rcout<<"Initialize pimat " <<std::endl;
mat pimat = zeros(n,K);
int j;
int k;
for(j=0;j<n;j++){
 for(k=0; k<K;k++){

//Rcpp::Rcout<<"Entering loop " <<std::endl;
//Rcpp::Rcout<<"Sigma Slice is "<< Sigma.slice(k) <<std::endl;
//Rcpp::Rcout<<"X mu trans part"<< solve(Sigma.slice(k),trans(X.row(j)-trans(mu.col(k))))<<std::endl;
//Rcpp::Rcout<<"solve portion working?"<< solve(Sigma.slice(k),trans(X.row(j)-trans(mu.col(k))))<<std::endl;

    pimat.submat(j,k,j,k)=pi(k)*pow(det(Sigma.slice(k)),-.5)*exp(-.5*(X.row(j)-trans(mu.col(k)))*solve(Sigma.slice(k),trans(X.row(j)-trans(mu.col(k)))));
  }
}
return wrap(pimat);
'


#Note that this assume mu matrix is transpose of how it is defined
fast.pi.mat= cxxfunction(signature(nS = "integer", mS = "integer", KS = "integer",
                                   XS = "double", muS="double",
                                   SigmaS = "numeric",
                                   piS = "double"
                                   ),
                         plugin = "RcppArmadillo",
                         body = src)



predict.normmix = function(X ,covar , mean , K){
  #covar is a m x m x K x iterates array, with the 3rd dimension indexing cluster, 4th dimension indexing MCMC draw
  #mean is a m x K x iterations
  #mean
  n = dim(X)[1]
  iterates = dim(covar)[4]
  predicted.value.mat = matrix(0, iterates, n)
  #to get the predicted values, for each vector we're going to get the probability
  #probability.mat  = mat(0,n,K)
  probability.vec = rep(0,K)
  print(paste("Iteration ", 0," at ", Sys.time()))
  #for each of the x variates, get a predicted value f
  for(i in 1:iterates){
    if(i%%100==0){
    print(paste("Iteration ", i," at ", Sys.time()))
    }
    #for each of the iterates, predict the value corresponding to each of the holdout values (maybe this should be getting done in C?)
    for(j in 1:n){
      #print(paste("obseration ", j," at ", Sys.time()))
        for(k in 1:K){
          probability.vec[k] = dmnorm(X[j,], mean[,k,i],covar[,,k,i])       
        } 
      if(sum(probability.vec)==0){
        for(k in 1:K){
          probability.vec[k] = dmnorm(X[j,], mean[,k,i],covar[,,k,i],log=T) 
        }
        probability.vec[which.max(probability.vec)]=1
        probability.vec = exp(probability.vec)
      }
      #print("Unscaled probability vector: ");print(probability.vec)
      probability.vec = probability.vec/sum(probability.vec)
      #print("Scaled probability vector: ");print(probability.vec)
      predicted.value.mat[i,j] = which(rmultinom(n=1,size=1,prob = probability.vec)==1,arr.ind=T)[1]
    }
  }
  return(predicted.value.mat)
}



 
predict.submix = function(P,theta,X,K){
  #Given the list of lists of P[[iter]][[k]] as the projection matrices
  #and theta[[iter]][[k]] as the list of thetas
  #assign observations X[n,m] to K clusters. 
  #While we are here, give a confidence interval for the distance of subspaces by calculating the conway distance at each iter
  
  n = dim(X)[1];
  m = dim(X)[2];
  iterations = length(P);
  #store cluster assignments
  cluster.mat = matrix(0,nrow=iterations, ncol=n);
  #store distances for each iteration
  dist.mat= matrix(0,nrow=K,ncol=n);
  #store each pairwaise distance
  n.pairs = K*(K-1)/2
  subdists = matrix(0,nrow = iterations, ncol = n.pairs);
  print(paste("Iteration ", 0," at ", Sys.time()))
  for(i in 1:iterations){
    if(i%%100==0){
      print(paste("Iteration ", i," at ", Sys.time()))
    }
    subdist.col = 1 #getting pairwise entry for subdist.col 1 first, incrementing from there
    for(j in 1:(K-1)){
      for(h in (j+1):K){
        subdists[i,subdist.col] = distance(P[[i]][[j]],P[[i]][[h]],subspace=F)
        subdist.col=subdist.col+1;
      }
    }
    #Now get the distance for each observation from the subspace
    for(k in 1:K){
      dist.mat[k,] = gibbs.loss.prj(x=X,P=P[[i]][[k]],mu=theta[[i]][[k]],subspace=F)
    }
    #now determine which cluster was closest to the point
    
    cluster.mat[i,]=apply(dist.mat,2,which.min)
  }
  return(list("cluster"=cluster.mat,"sub.dist"=subdists))
  
}


