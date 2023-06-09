
source("MCS.R")

library(trust)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(goeveg)

load("inputdata.RData") 
ttwhole = tt$freqdata
bb = tt$trueprob
depthname<-c("DP.treat1","DP.treat2","DP.treat3","DP.DMSO1","DP.DMSO2","DP.DMSO3")
celllines = c("C3A", "Huh7", "SNU449", "SNU475", "PLC", "SK-Hep1")
load("estiresult.RData") 

##basic info
nsnp<-c(100,200,400,800,1000,1500)
set.seed(123)
mysample<-sample(1:nrow(ttwhole),1500,replace = F)
set.seed(NULL)
samplett<-ttwhole
samplebb<-bb
subsamplett<-list()
subsamplebb<-list()
##create subsamples to do estimation for different sample sizes
for (i in 1:length(nsnp)) {
  samplett<-ttwhole[mysample[1:nsnp[i]],]
  samplebb<-bb[mysample[1:nsnp[i]],]
  subsamplett[[i]]<-samplett
  subsamplebb[[i]]<-samplebb
}

####
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
seed=12345
nsimu=1000
#########sample size 100
##get the simulation results
k=1
simulationdata<-list()
subfreqdata<-ttwhole[1:nsnp[k],]
for (i in 1:nsimu) {
  seed=seed+i
  for (j in 1:(length(celllines)/2)) {
    
    subfreqdata[,(j-1)*2+1]<-subsamplett[[k]][,depthname[j]]
    subfreqdata[,j*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,1],subsamplett[[k]],subsamplett[[k]][,depthname[j]],seed)
    subfreqdata[,(j+2)*2+1]<-subsamplett[[k]][,depthname[j+3]]
    subfreqdata[,(j+3)*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,2],subsamplett[[k]],subsamplett[[k]][,depthname[(j+3)]],seed)
    
  }
  
  simulationdata[[i]]<-list(trueprob=subsamplebb[[k]],freqdata=subfreqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
subresult = list()

for (i in 1:nsimu){
  
  subresult[[i]] = MCS(simulationdata[[i]])
  
}
save(subresult,file = "subresult100.RData")
#########sample size 100
#########sample size 200
##get the simulation results
k=2
simulationdata<-list()
subfreqdata<-ttwhole[1:nsnp[k],]
for (i in 1:nsimu) {
  seed=seed+i
  for (j in 1:(length(celllines)/2)) {
    
    subfreqdata[,(j-1)*2+1]<-subsamplett[[k]][,depthname[j]]
    subfreqdata[,j*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,1],subsamplett[[k]],subsamplett[[k]][,depthname[j]],seed)
    subfreqdata[,(j+2)*2+1]<-subsamplett[[k]][,depthname[j+3]]
    subfreqdata[,(j+3)*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,2],subsamplett[[k]],subsamplett[[k]][,depthname[(j+3)]],seed)
    
  }
  
  simulationdata[[i]]<-list(trueprob=subsamplebb[[k]],freqdata=subfreqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
subresult = list()

for (i in 1:nsimu){
  
  subresult[[i]] = MCS(simulationdata[[i]])
  
}
save(subresult,file = "subresult200.RData")
#########sample size 200
#########sample size 400
##get the simulation results
k=3
simulationdata<-list()
subfreqdata<-ttwhole[1:nsnp[k],]
for (i in 1:nsimu) {
  seed=seed+i
  for (j in 1:(length(celllines)/2)) {
    
    subfreqdata[,(j-1)*2+1]<-subsamplett[[k]][,depthname[j]]
    subfreqdata[,j*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,1],subsamplett[[k]],subsamplett[[k]][,depthname[j]],seed)
    subfreqdata[,(j+2)*2+1]<-subsamplett[[k]][,depthname[j+3]]
    subfreqdata[,(j+3)*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,2],subsamplett[[k]],subsamplett[[k]][,depthname[(j+3)]],seed)
    
  }
  
  simulationdata[[i]]<-list(trueprob=subsamplebb[[k]],freqdata=subfreqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
subresult = list()

for (i in 1:nsimu){
  
  subresult[[i]] = MCS(simulationdata[[i]])
  
}
save(subresult,file = "subresult400.RData")
#########sample size 400
#########sample size 800
##get the simulation results
k=4
simulationdata<-list()
subfreqdata<-ttwhole[1:nsnp[k],]
for (i in 1:nsimu) {
  seed=seed+i
  for (j in 1:(length(celllines)/2)) {
    
    subfreqdata[,(j-1)*2+1]<-subsamplett[[k]][,depthname[j]]
    subfreqdata[,j*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,1],subsamplett[[k]],subsamplett[[k]][,depthname[j]],seed)
    subfreqdata[,(j+2)*2+1]<-subsamplett[[k]][,depthname[j+3]]
    subfreqdata[,(j+3)*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,2],subsamplett[[k]],subsamplett[[k]][,depthname[(j+3)]],seed)
    
  }
  
  simulationdata[[i]]<-list(trueprob=subsamplebb[[k]],freqdata=subfreqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
subresult = list()

for (i in 1:nsimu){
  
  subresult[[i]] = MCS(simulationdata[[i]])
  
}
save(subresult,file = "subresult800.RData")
#########sample size 800
#########sample size 1000
##get the simulation results
k=5
simulationdata<-list()
subfreqdata<-ttwhole[1:nsnp[k],]
for (i in 1:nsimu) {
  seed=seed+i
  for (j in 1:(length(celllines)/2)) {
    
    subfreqdata[,(j-1)*2+1]<-subsamplett[[k]][,depthname[j]]
    subfreqdata[,j*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,1],subsamplett[[k]],subsamplett[[k]][,depthname[j]],seed)
    subfreqdata[,(j+2)*2+1]<-subsamplett[[k]][,depthname[j+3]]
    subfreqdata[,(j+3)*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,2],subsamplett[[k]],subsamplett[[k]][,depthname[(j+3)]],seed)
    
  }
  
  simulationdata[[i]]<-list(trueprob=subsamplebb[[k]],freqdata=subfreqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
subresult = list()

for (i in 1:nsimu){
  
  subresult[[i]] = MCS(simulationdata[[i]])
  
}
save(subresult,file = "subresult1000.RData")
#########sample size 1000
#########sample size 1500
##get the simulation results
k=6
simulationdata<-list()
subfreqdata<-ttwhole[1:nsnp[k],]
for (i in 1:nsimu) {
  seed=seed+i
  for (j in 1:(length(celllines)/2)) {
    
    subfreqdata[,(j-1)*2+1]<-subsamplett[[k]][,depthname[j]]
    subfreqdata[,j*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,1],subsamplett[[k]],subsamplett[[k]][,depthname[j]],seed)
    subfreqdata[,(j+2)*2+1]<-subsamplett[[k]][,depthname[j+3]]
    subfreqdata[,(j+3)*2]<-simudata(subsamplebb[[k]],estiresult$esti_prop[,2],subsamplett[[k]],subsamplett[[k]][,depthname[(j+3)]],seed)
    
  }
  
  simulationdata[[i]]<-list(trueprob=subsamplebb[[k]],freqdata=subfreqdata,group=c(1,1,1,0,0,0),
                            samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))
  
}
set.seed(NULL)
##get the estimated results for 1000 simulations
subresult = list()

for (i in 1:nsimu){
  
  subresult[[i]] = MCS(simulationdata[[i]])
  
}
save(subresult,file = "subresult1500.RData")
#########sample size 1500
