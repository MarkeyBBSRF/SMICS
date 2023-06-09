library(trust)
deconvolute <- function(trueprob, groupdata, n.ini){
  ####the mle part######
  ###write my function
  loglikelihood<-function(v){
    
    transprop<-(c(exp(v)/(sum(exp(v))+1),1/(sum(exp(v))+1)))
    sumexprob <- rowSums(t(transprop*t(trueprob)))
    pmatrix<-matrix(c(log(sumexprob),c(log(1-sumexprob))),ncol = 1)
    result=0
    for (i in 1:length(groupdata)) {
      assign(paste("nsmatrix",i,sep = ''),matrix(c(groupdata[[i]]$mutcount,
                                                   groupdata[[i]]$DP-groupdata[[i]]$mutcount),nrow = 1))
      
      result=result+as.numeric((-1)*get(paste("nsmatrix",i,sep = ''))%*%pmatrix)
    }
    dl = calc.dl(v)
    ddl = calc.ddl(v)
#    print(result)
#    print(transprop)
    return(list(value=result, gradient=dl, hessian=ddl))
    #    return(result)
  }
  
  calc.dl<-function(v){
    q = trueprob           #i * k
    r<-(c(exp(v)/(sum(exp(v))+1),1/(sum(exp(v))+1)))  #transprop   k vector
    K = length(r)
    p <- c(t(r) %*% t(q))            #sumexprob    i vector
    dl = 0
    for (j in 1:length(groupdata)){
      x = groupdata[[j]]$mutcount
      n = groupdata[[j]]$DP
      tempa = x/p - (n-x)/(1-p)
      
      part1 = t(tempa) %*% (q %*% r)    #1*1 
      part2 = c(t(tempa) %*% q)           # k vector
      
      vec1 = matrix(part1, K, 1)
      vec2 = matrix(part2, K, 1)
      dl = dl + r * (vec2 - vec1)
    }
    return((-1)*dl[1:(K-1),])   # for -loglikelihood
  }
  calc.ddl <- function(v){
    q = trueprob           #i * k
    r<-(c(exp(v)/(sum(exp(v))+1),1/(sum(exp(v))+1)))  #transprop   k vector
    K = length(r)
    p <- c(t(r) %*% t(q))            #sumexprob    i vector
    ddl = 0
    for (j in 1:length(groupdata)){
      x = groupdata[[j]]$mutcount
      n = groupdata[[j]]$DP
      tempa = x/p - (n-x)/(1-p)
      tempb = (-1) * (x/p/p + (n-x)/(1-p)/(1-p))
      
      part1 = t(tempa) %*% (q %*% r)    #1*1 
      part2 = c(t(tempa) %*% q)           # k vector
      part3 = c(t(tempb * p) %*% q)    #k vector 
      part4 = t(tempb * p) %*% p    #1*1
      
      mat1 = matrix(part1, K, K)
      mat2 = matrix(part2, K, K)
      mat3 = matrix(part3, K, K)
      mat4 = matrix(part4, K, K)
      mat5 = t(tempb * q) %*% q  
      
      rmat = matrix(r, K, K)
      dpart = rmat * (mat5 - mat3 - t(mat3) + mat4) * t(rmat)
      ddpart = rmat * (2*mat1 - mat2 - t(mat2)) * t(rmat)
      vec1 = matrix(part1, K, 1)
      vec2 = matrix(part2, K, 1)
      diag(ddpart) = r * (1-2*r) * (vec2 - vec1)
      
      ddl = ddl + dpart + ddpart  
    }
    return((-1)*ddl[1:(K-1), 1:(K-1)])   #for -loglikelihood
  }
  
  ff = list()
  ml = NULL
  set.seed(1234568)
  ini.value = matrix(runif(n.ini*(ncol(trueprob)-1)),nrow=n.ini)
  for (j in 1:n.ini){
    #    ff[[j]] = optim(ini.value[j,],loglikelihood,hessian=TRUE)
    ff[[j]] = trust(loglikelihood, ini.value[j,],rinit=1, rmax=100, iterlim=1000 )
    ml[j] = ff[[j]]$value
  }
  
  fit = ff[[which.min(ml)]]
  fisher_info<-solve(fit$hessian)
  prop_sigma<-sqrt(diag(fisher_info))
  prop_sigma<-diag(prop_sigma)
  #  v<-fit$par
  v <- fit$argument
#  print(fit)
  
  ##############do the transformation
  ##set a matrix to contain the estimators and var
  pro<-matrix(0,nrow =ncol(trueprob) ,ncol = 3)
  
  ##get the estimators of proportion
  for (i in 1:(ncol(trueprob)-1)) {
    pro[i,1]=exp(v[i])/(sum(exp(v))+1)
  }
  pro[ncol(trueprob),1]=1/(sum(exp(v))+1)
  return(list(pro=pro,  v=v, prop_sigma=prop_sigma, hessian=fit$hessian,fml=ml))
}
