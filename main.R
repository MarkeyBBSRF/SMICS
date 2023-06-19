## Import the cell line mixture sequencing data. The data should be a list with the following 
## four elements
## trueprob: a data frame listing the frequency of each SNP, i.e. 0, 0.5 or 1, in each cell line
## freqdata: a data frame for observed SNP data, where there are two columns for each sample.
##           The first column (DP) is the read depth and the second column (AD) is the number 
##           of reads supporting the alternative allele
## group: a vector of treatment group indicators of the samples
## samples: a vector of sample names
## Copyright (2023) University of Kentucky.
load("inputdata.RData") 

## Load the function for performing cell line mixture deconvolution
source("MCS.R")

library(trust)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(goeveg)
############################################## Begin of fig4  ##############################################
estiresult<-MCS(tt)
#save(estiresult, file="estiresult.RData")

##figure A - get confidence interval
ngroup<-2
celllines = c("C3A", "Huh7", "SNU449", "SNU475", "PLC", "SK-Hep1")
upperbound_prop<-matrix(nrow = length(celllines),ncol = ngroup)
for (i in 1:length(celllines)) {
  for (j in 1:ngroup) {
    upperbound_prop[i,j]<-estiresult$esti_prop[i,j]+1.96*estiresult$sd_prop[i,j]
    
  }
  
}
upperbound_prop<-round(upperbound_prop,3)

lowerbound_prop<-matrix(nrow = length(celllines),ncol = ngroup)
for (i in 1:length(celllines)) {
  for (j in 1:ngroup) {
    lowerbound_prop[i,j]<-estiresult$esti_prop[i,j]-1.96*estiresult$sd_prop[i,j]
    
  }
  
}

lowerbound_prop<-round(lowerbound_prop,3)

bound_inhib<-matrix(nrow = length(celllines),ncol = 2)
for (i in 1:length(celllines)) {
  
  bound_inhib[i,1]<-estiresult$esti_inhibition [i]-1.96*estiresult$sd_inhibition[i]
  bound_inhib[i,2]<-estiresult$esti_inhibition [i]+1.96*estiresult$sd_inhibition[i]
  
}
#bound_inhib<-round(bound_inhib,3)


mytable<-matrix(nrow = length(celllines),ncol = 3)
rownames(mytable)<-celllines
for (i in 1:length(celllines)) {
  
  mytable[i,1]<-paste(round(estiresult$esti_prop[i,1],3),"(",lowerbound_prop[i,1],",",upperbound_prop[i,1],")",sep = " ")
  mytable[i,2]<-paste(round(estiresult$esti_prop[i,2],3),"(",lowerbound_prop[i,2],",",upperbound_prop[i,2],")",sep = " ")
  mytable[i,3]<-paste(round(estiresult$esti_inhibition[i]*359500/450000,3),
                      "(",round(bound_inhib[i,1]*359500/450000,3),
                      ",",round(bound_inhib[i,2]*359500/450000,3),")",sep = " ")
}

print(mytable)

##figure B - plot of cell count ratio
realdata<-c(0.49,0.56,0.93,0.93,0.86,1.18)

plotdata<-data.frame(sequenced=estiresult$esti_inhibition*359500/450000,
                     observed=realdata,cellline=celllines)


ggplot(data=plotdata,aes(x=sequenced,y=observed))+
  geom_point(aes(colour=cellline),size=3)+
  scale_color_discrete(name = "Cell line")+
  xlim(0.45,1.18)+
  ylim(0.45,1.18)+
  xlab("Ratio of cell counts from SMICS")+
  ylab("Ratio of cell counts from \n cell proliferation assay")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.background = element_rect(fill="transparent"),text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),plot.margin = unit(c(0.1,0,0.1,0.2), "cm"),
        legend.title = element_text(size=12),legend.text = element_text(size=12),legend.key.height = unit(0.4,'cm'))+
  guides(fill= guide_legend(nrow = 1))+
  
  geom_abline(intercept = 0, slope = 1, linetype="dashed")


##figure C 
estidata<-estiresult$esti_inhibition*359500/450000
realsd<-c(0.01, 0.06, 0.19, 0.23, 0.09, 0.21)
seqse<-estiresult$sd_inhibition*359500/450000

dt<-data.frame(X=c(1:6),estratio=estidata,estsd=seqse,realratio=realdata,realsd=realsd)
dt$X=c("C3A", "Huh7", "SNU449","SNU475", "PLC", "SK-Hep1")
dt$X = factor(dt$X,levels=dt$X)

dt
colnames(dt) = c("X","estratio","estsd","realratio","realsd")
dt2 = dt %>% dplyr::select(1,3,5)
dt2 = reshape2::melt(dt2,id.vars="X")

ggplot(data=dt2, aes(x=X, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+scale_y_continuous(expand = expansion(mult = c(0, .15)))+
  scale_fill_manual(name=NULL,values=c("#F0E442", "#0072B2"),breaks=c("estsd", "realsd"),labels=c("SMICS", "Cell proliferation assay"))+
  ylab("Standard error")+
  xlab(NULL)+
  theme_bw() +
  theme(legend.position = c(0.5, .95),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.background = element_rect(fill="transparent"),text=element_text(size=14))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size = 14),plot.margin = unit(c(0,0,0,0), "cm"),
        legend.title = element_text(size=14),legend.text = element_text(size=14),legend.key.height = unit(0.4,'cm'))+
  guides(fill= guide_legend(nrow = 1))

############################################## End of fig4  ############################################


############################################## Begin of fig5  ##############################################

## Figure 5a
## First load the simulation results. The original simulation code is provided in simulation.R

load("finalresult0106.RData")

mycellline = c("C3A", "Huh7", "SNU449", "SNU475", "PLC", "SK-Hep1")
nsimu=1000
###check cp
inhibprop<-estiresult$esti_prop[,1]/estiresult$esti_prop[,2]
inhibitionresult<-matrix(nrow = nsimu,ncol = length(inhibprop))
for (i in 1:nsimu) {
  
  inhibitionresult[i,]<-finalresult[[i]]$esti_prop[,1]/finalresult[[i]]$esti_prop[,2]
  
}

coverage<-matrix(0,nrow=nsimu, ncol=length(inhibprop))
bias = matrix(0,nrow=nsimu, ncol=length(inhibprop))
#########change from "estiresult$sd_inhibition[j]" to "finalresult[[i]]$sd_inhibition[i,j]"
for (i in 1:nsimu){
  bias[i,] = inhibitionresult[i,] - inhibprop
  for (j in 1:length(inhibprop)) {
    if((inhibprop[j]<(inhibitionresult[i,j]+1.96*finalresult[[i]]$sd_inhibition[j]))&(inhibprop[j]>(inhibitionresult[i,j]-1.96*finalresult[[i]]$sd_inhibition[j])))
    {coverage[i,j]=1}
  }
}

mybias_inhib<-abs(colMeans(bias))
trust_mycp_inhib<-colMeans(coverage)
trust_mycp_inhib

##rmse
rmse<-rep(0,length(inhibprop))
for (j in 1:length(inhibprop)) {
  
  rmse[j]<-sqrt(mean((inhibitionresult[,j]-inhibprop[j])^2))
  
}

##treat
##sd of group 1
coverage_trt<-matrix(0,nrow=nsimu, ncol=length(inhibprop))
bias_trt = matrix(0,nrow=nsimu, ncol=length(inhibprop))

##cp

for (i in 1:nsimu) {
  bias_trt[i,] =finalresult[[i]]$esti_prop[,1]-estiresult$esti_prop[,1]
  
  for (j in 1:length(mycellline)) {
    if((estiresult$esti_prop[j,1]<finalresult[[i]]$esti_prop[j,1]+1.96*finalresult[[i]]$sd_prop[j,1])&(estiresult$esti_prop[j,1]>finalresult[[i]]$esti_prop[j,1]-1.96*finalresult[[i]]$sd_prop[j,1]))      
    {coverage_trt[i,j]=1}
  }
  
}
mybias_trt<-(abs(colMeans(bias_trt)))
mycp_trt<-(colMeans(coverage_trt))
mycp_trt
##rmse
rmse_trt<-rep(0,length(inhibprop))
for (j in 1:length(inhibprop)) {
  
  rmse_trt[j]<-sqrt(mean((bias_trt[,j])^2))
  
}

##dmso
##sd of group 2
coverage_dmso<-matrix(0,nrow=nsimu, ncol=length(inhibprop))
bias_dmso = matrix(0,nrow=nsimu, ncol=length(inhibprop))

##cp
for (i in 1:nsimu) {
  bias_dmso[i,] =finalresult[[i]]$esti_prop[,2]-estiresult$esti_prop[,2]
  for (j in 1:length(mycellline)) {
    if((estiresult$esti_prop[j,2]<finalresult[[i]]$esti_prop[j,2]+1.96*finalresult[[i]]$sd_prop[j,2])&(estiresult$esti_prop[j,2]>finalresult[[i]]$esti_prop[j,2]-1.96*finalresult[[i]]$sd_prop[j,2]))      
    {coverage_dmso[i,j]=1}
  }
  
}
mybias_dmso<-(abs(colMeans(bias_dmso)))
mycp_dmso<-(colMeans(coverage_dmso))
mycp_dmso
##rmse
rmse_dmso<-rep(0,length(inhibprop))
for (j in 1:length(inhibprop)) {
  
  rmse_dmso[j]<-sqrt(mean((bias_dmso[,j])^2))
  
}
mytable<-data.frame(biastrt=mybias_trt,cptrt=mycp_trt,rmsetrt=round(rmse_trt,4), biasdmso=mybias_dmso,cpdmso=mycp_dmso,rmsedmso=round(rmse_dmso,4),
                    biasinhib=mybias_inhib,cpinhib=trust_mycp_inhib,rmseinhib=round(rmse,4))
rownames(mytable)<-celllines
print(mytable)

## Figure 5b

## Load simulation results. The original simulation code is provided in subsample of 6 celllines.R
path<-"./subresult"
nsimu=1000
trueprop<-c(0.2,0.35,0.1,0.1,0.2,0.05)
nsnp<-c(100,200,400,800,1000,1500)



trtlist<-list()
for (j in 1:length(nsnp)) {
  load(paste(path,nsnp[j],".RData",sep = ""))
  trtlist[[j]]<-matrix(nrow = nsimu,ncol = length(celllines))
  for (i in 1:nsimu) {
    trtlist[[j]][i,]<-subresult[[i]][["esti_prop"]][,1]
    
  }
  
}

dmsolist<-list()
for (j in 1:length(nsnp)) {
  load(paste(path,nsnp[j],".RData",sep = ""))
  dmsolist[[j]]<-matrix(nrow = nsimu,ncol = length(celllines))
  for (i in 1:nsimu) {
    dmsolist[[j]][i,]<-subresult[[i]][["esti_prop"]][,2]
    
  }
  
}


inhib<-list()
for (j in 1:length(nsnp)) {
  load(paste(path,nsnp[j],".RData",sep = ""))
  inhib[[j]]<-matrix(nrow = nsimu,ncol = length(celllines))
  for (i in 1:nsimu) {
    inhib[[j]][i,]<-subresult[[i]][["esti_prop"]][,1]/subresult[[i]][["esti_prop"]][,2]
    
  }
  
}

######a567
mycv1<-matrix(ncol=length(celllines),nrow = length(nsnp))

for (i in 1:length(nsnp)) {
  mycv1[,i]<-apply(trtlist[[i]],2,cv)
  
}

mycv1<-as.vector(mycv1)

mygroup<-matrix(ncol=length(celllines),nrow=length(nsnp))

for (i in 1:length(nsnp)) {
  mygroup[,i]<-rep(nsnp[i],length(celllines))
  
}
mygroup<-as.vector(mygroup)


myratio1<-data.frame(cv=mycv1,
                     cellline=celllines,
                     group=mygroup)
plotdata<-matrix(rep(0,length(nsnp)*(length(celllines)+1)),ncol = length(celllines)+1)
for (i in 1:length(celllines)) {
  plotdata[,i]<-myratio1[which(myratio1[,2]==celllines[i]),1]
  
}
plotdata[,length(celllines)+1]<-1:length(nsnp)
plotdata<-as.data.frame(plotdata)
colnames(plotdata)<-c("DEDA1","DEDA2","DEDA3","DMSO1","DMSO2","DMSO3","group")
myratio1$cellline<-as.factor(myratio1$cellline)
library(ggplot2)
suppressWarnings(library(ggplot2))
##plot parameter setting
ggplot()+ 
  geom_line(plotdata,mapping=aes(group,DEDA1),color="#F8766D")+
  geom_line(plotdata,mapping=aes(group,DEDA2),color="#B79F00")+ 
  geom_line(plotdata,mapping=aes(group,DEDA3),color="#619CFF")+
  geom_line(plotdata,mapping=aes(group,DMSO1),color="#F564E3")+
  geom_line(plotdata,mapping=aes(group,DMSO2),color="#00BA38")+
  geom_line(plotdata,mapping=aes(group,DMSO3),color="#00BFC4")+
  geom_point(data=myratio1,mapping = aes(x=factor(group,levels = c("100","200","400","800","1000","1500")),
                                         y=cv,colour=cellline),size=2)+
  ylim(0,0.32)+
  ggtitle("DEDA")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),legend.position = "none",axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=17,face="bold"))+
  scale_x_discrete()+
  labs(x="Number of SNPs",y="Coefficient of Variation")


####dmso
mycv2<-matrix(ncol=length(celllines),nrow = length(nsnp))

for (i in 1:length(nsnp)) {
  mycv2[,i]<-apply(dmsolist[[i]],2,cv)
  
}
mycv2<-as.vector(mycv2)

myratio2<-data.frame(cv=mycv2,
                     cellline=celllines,
                     group=mygroup)

plotdata<-matrix(rep(0,length(nsnp)*(length(celllines)+1)),ncol = length(celllines)+1)
for (i in 1:length(celllines)) {
  plotdata[,i]<-myratio2[which(myratio2[,2]==celllines[i]),1]
  
}
plotdata[,length(celllines)+1]<-1:length(nsnp)
plotdata<-as.data.frame(plotdata)
colnames(plotdata)<-c("DEDA1","DEDA2","DEDA3","DMSO1","DMSO2","DMSO3","group")
myratio2$cellline<-as.factor(myratio2$cellline)

##plot parameter setting
ggplot()+ 
  geom_line(plotdata,mapping=aes(group,DEDA1),color="#F8766D")+
  geom_line(plotdata,mapping=aes(group,DEDA2),color="#B79F00")+ 
  geom_line(plotdata,mapping=aes(group,DEDA3),color="#619CFF")+
  geom_line(plotdata,mapping=aes(group,DMSO1),color="#F564E3")+
  geom_line(plotdata,mapping=aes(group,DMSO2),color="#00BA38")+
  geom_line(plotdata,mapping=aes(group,DMSO3),color="#00BFC4")+
  geom_point(data=myratio2,mapping = aes(x=factor(group,levels = c("100","200","400","800","1000","1500")),
                                         y=cv,colour=cellline),size=2)+
  ylim(0,0.32)+
  ggtitle("DMSO")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),legend.position = "none",axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=17,face="bold"))+
  scale_x_discrete()+
  labs(x="Number of SNPs",y="Coefficient of Variation")


###ratio
mycv3<-matrix(ncol=length(celllines),nrow = length(nsnp))

for (i in 1:length(nsnp)) {
  mycv3[,i]<-apply(inhib[[i]],2,cv)
  
}

mycv3<-as.vector(mycv3)

myratio3<-data.frame(cv=mycv3,
                     cellline=celllines,
                     group=mygroup)

plotdata<-matrix(rep(0,length(nsnp)*(length(celllines)+1)),ncol = length(celllines)+1)
for (i in 1:length(celllines)) {
  plotdata[,i]<-myratio3[which(myratio3[,2]==celllines[i]),1]
  
}
plotdata[,length(celllines)+1]<-1:length(nsnp)
plotdata<-as.data.frame(plotdata)
colnames(plotdata)<-c("DEDA1","DEDA2","DEDA3","DMSO1","DMSO2","DMSO3","group")
myratio3$cellline<-as.factor(myratio3$cellline)

##plot parameter setting
ggplot()+ 
  geom_line(plotdata,mapping=aes(group,DEDA1),color="#F8766D")+
  geom_line(plotdata,mapping=aes(group,DEDA2),color="#B79F00")+ 
  geom_line(plotdata,mapping=aes(group,DEDA3),color="#619CFF")+
  geom_line(plotdata,mapping=aes(group,DMSO1),color="#F564E3")+
  geom_line(plotdata,mapping=aes(group,DMSO2),color="#00BA38")+
  geom_line(plotdata,mapping=aes(group,DMSO3),color="#00BFC4")+
  geom_point(data=myratio3,mapping = aes(x=factor(group,levels = c("100","200","400","800","1000","1500")),
                                         y=cv,colour=cellline),size=2)+
  ylim(0,0.32)+
  ggtitle("Ratio of cell counts")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),legend.position = "none",axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=17,face="bold"))+
  scale_x_discrete()+
  labs(x="Number of SNPs",y="Coefficient of Variation")


##############################################  End of fig5  ##############################################


############################################# Begin of Figure 6 ###########################################

## Load simulation results. The original simulation code is provided in NCI60 sd011022.R
library(ggplot2)
dt = read.csv("plot of NCI60_0109.csv",stringsAsFactors = F)
ggplot(dt,aes(x=esti,y=true,color=cancertype)) + geom_point(size=.5) +
  geom_errorbar(aes(ymin=(true-sd), ymax=true+sd), width=.02,size=.3,position=position_dodge(.9)) +
  labs(x = "Estimated Ratio of cell counts", y = "True Ratio of cell counts",color= "Cell line type") + theme_classic()
############################################# End of Figure 6 ###########################################













