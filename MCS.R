##01/06/2022:change line 89 "result[[2]]"

source("deconvolute.R")
MCS <- function(data){
  ##################################get the treat and dmso group data##################################
  
  ind=which(data$group==1)
  trtfreqdata<-list()
  for (i in 1:length(ind)) { 
    
    trtfreqdata[[i]]<-data.frame(DP=data$freqdata[,paste("DP", data$samples[ind[i]], sep=".")],
                                 mutcount=data$freqdata[,paste("AD", data$samples[ind[i]], sep=".")])
    
  }
  ind=which(data$group==0)
  dmsofreqdata<-list()
  for (i in 1:length(ind)) { 
    
    dmsofreqdata[[i]]<-data.frame(DP=data$freqdata[,paste("DP", data$samples[ind[i]], sep=".")],
                                  mutcount=data$freqdata[,paste("AD", data$samples[ind[i]], sep=".")])
    
  }
  ##################################get the treat and dmso group data##################################
  ###################################get the cell_ proportion results and corresponding hessian matrix#################################
  result = list()
  
  groupname<-c("trtfreqdata","dmsofreqdata")
  
  for (j in 1:2){
    
    result[[j]] = deconvolute(as.matrix(data$trueprob),get(groupname[j]),1)
  }
  ###################################get the cell_ proportion results and corresponding hessian matrix#################################
  ###################################calculate the sd matrix of inhibition###################################
  mycellline<-colnames(tt[[1]])
  
  mysum1<-sum(exp(result[[1]]$v))+1
  mysum2<-sum(exp(result[[2]]$v))+1
  
  r1<-exp(result[[1]]$v)/mysum1
  r2<-exp(result[[2]]$v)/mysum2
  
  ##g1j
  giw1j<-matrix(nrow = length(mycellline)-1,ncol = length(mycellline)-1)
  for (i in 1:(length(mycellline)-1)) {
    for (j in 1:(length(mycellline)-1)) {
      if (i==j){ 
        giw1j[i,j]<-r1[i]*(1-r1[i])/r2[i]
        
      }else{
        giw1j[i,j]<- -r1[j]*r1[i]/r2[i]}
      
    }
    
  }
  ##g2j
  giw2j<-matrix(nrow = length(mycellline)-1,ncol = length(mycellline)-1)
  for (i in 1:(length(mycellline)-1)) {
    for (j in 1:(length(mycellline)-1)) {
      if(i==j){
        giw2j[i,j]<-r1[i]*(-1/r2[i]+1)
        
      }else{
        giw2j[i,j]<- r1[i]*(r2[j]/r2[i])
        
      }
      
    }
    
  }
  
  
  ##gK
  gK<-rep(0,2*(length(mycellline)-1))
  for (i in 1:(length(mycellline)-1)) {
    gK[i]<- -r1[i]*(mysum2/mysum1)
    gK[i+(length(mycellline)-1)]<- r2[i]*(mysum2/mysum1)
    
  }
  
  
  delg<-cbind(giw1j,giw2j)
  delg<-rbind(delg,gK)
  
  
  ##get the original matrix
  myoriginal<-matrix(rep(0,(2*(length(mycellline)-1))*2*(length(mycellline)-1)),
                     nrow = 2*(length(mycellline)-1),ncol = 2*(length(mycellline)-1))
  myoriginal[1:(length(mycellline)-1),1:(length(mycellline)-1)]<-solve(result[[1]]$hessian)
  myoriginal[length(mycellline):(2*(length(mycellline)-1)),length(mycellline):(2*(length(mycellline)-1))]<-solve(result[[2]]$hessian)
  
  
  
  ##sd matrix
  newsd<-delg%*%myoriginal%*%t(delg)
  #####get sd of inhibition#######
  ###################################calculate the sd matrix of inhibition###################################
  ###################################calculate the sd matrix of proportion###################################
  
  ##sd of group 1
  gii_1<-matrix(nrow = length(mycellline)-1,ncol = length(mycellline)-1)
  for (i in 1:(length(mycellline)-1)) {
    for (j in 1:(length(mycellline)-1)) {
      if (i==j){ 
        gii_1[i,j]<- r1[i]*(1-r1[i])
        
      }else{
        gii_1[i,j]<- -r1[j]*r1[i]}
      
    }
    
  }
  
  gK_1<-rep(0,(length(mycellline)-1))
  for (i in 1:(length(mycellline)-1)) {
    gK_1[i]<- -r1[i]/mysum1
    
  }
  delg_1<-rbind(gii_1,gK_1)
  
  myoriginal_1<-solve(result[[1]]$hessian)
  ##sd matrix
  sd_prop_1<-delg_1%*%myoriginal_1%*%t(delg_1)
  
  ##sd of group 2
  gii_2<-matrix(nrow = length(mycellline)-1,ncol = length(mycellline)-1)
  for (i in 1:(length(mycellline)-1)) {
    for (j in 1:(length(mycellline)-1)) {
      if (i==j){ 
        gii_2[i,j]<- r2[i]*(1-r2[i])
        
      }else{
        gii_2[i,j]<- -r2[j]*r2[i]}
      
    }
    
  }
  
  gK_2<-rep(0,(length(mycellline)-1))
  for (i in 1:(length(mycellline)-1)) {
    gK_2[i]<- -r2[i]/mysum2
    
  }
  delg_2<-rbind(gii_2,gK_2)
  
  myoriginal_2<-solve(result[[2]]$hessian)
  
  ##sd matrix
  sd_prop_2<-delg_2%*%myoriginal_2%*%t(delg_2)
  
  ###################################calculate the sd matrix of proportion###################################
  
  allresult=list(esti_prop=cbind(result[[1]]$pro[,1],result[[2]]$pro[,1]),esti_inhibition=result[[1]]$pro[,1]/result[[2]]$pro[,1],
                 sd_prop=cbind(sqrt(diag(sd_prop_1)),sqrt(diag(sd_prop_2))),sd_inhibition=sqrt(diag(newsd)),
                 tempresult=result)
  
  return(allresult)
  
  
}



