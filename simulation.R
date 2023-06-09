############# simulation data function ################

simudata <- function(trueprob, trueprop,data, depth, seed){
  set.seed(seed)
  ##get the sum probability of each position of all cell lines
  sumprob <- rowSums(t(trueprop*t(trueprob)))
  mutcount= rbinom(n=rep(1,length(data[,depthname[j]])), size=depth, prob=sumprob)
  
  ##get my simulated data
  
  return(mutcount)
}
############# simulation data function ################
source("MCS.R")
seed=12345
nsimu=1000
##get the simulation results
load("inputdata.RData") 
ttwhole = tt$freqdata
bb = tt$trueprob
depthname<-c("DP.treat1","DP.treat2","DP.treat3","DP.DMSO1","DP.DMSO2","DP.DMSO3")
celllines = c("C3A", "Huh7", "SNU449", "SNU475", "PLC", "SK-Hep1")
load("estiresult.RData") 

simulationdata<-list()
for (i in 1:nsimu) {
  seed=seed+i
  freqdata = matrix(nrow=nrow(ttwhole), ncol=length(celllines)*2)
  for (j in 1:(length(celllines)/2)) {
    
    freqdata[,(j-1)*2+1]<-ttwhole[,depthname[j]]
    freqdata[,j*2]<-simudata(bb,estiresult$esti_prop[,1],ttwhole,ttwhole[,depthname[j]],seed)
    freqdata[,(j+2)*2+1]<-ttwhole[,depthname[j+3]]
    freqdata[,(j+3)*2]<-simudata(bb,estiresult$esti_prop[,2],ttwhole,ttwhole[,depthname[(j+3)]],seed)
    
  }
  colnames(freqdata) = c("DP.treat1", "AD.treat1", "DP.treat2","AD.treat2", "DP.treat3","AD.treat3", 
                         "DP.DMSO1","AD.DMSO1", "DP.DMSO2","AD.DMSO2", "DP.DMSO3", "AD.DMSO3")
  
  simulationdata[[i]]<-list(trueprob=bb,freqdata=freqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
finalresult = list()

for (i in 1:nsimu){
  
  finalresult[[i]] = MCS(simulationdata[[i]])
  
}
save(finalresult,file = "finalresult0106.RData")