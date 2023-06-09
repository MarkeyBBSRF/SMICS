## DMSO
# proportion = 1/60
# assume 1000 cells, cell count for a cell line is 1000/60

## Treated
# cell count 1000/60 * 0.206 for BC:MCF7 cell line
# cell count 1000/60 * 0.4526 for BC:MDA-MB-231/ATCC cell line
# proportion = 1000/60 * 0.206 / (1000/60 * 0.206 + 1000/60 * 0.4526 + ...)
# for negative value, set it to 0.01

# Simulate a dataset with 3 replicates in each group. Set read depth = 100 at each SNP site
# SNP genotype: A = 0, H = 0.5, B =1

##import genotype for each cell line
aa = read.csv("DNA__Affy_500K_SNP_CRLMM.txt.txt", stringsAsFactors = F, 
              sep='\t', skip=13,header = F)
bb = aa[,9:68]
####so each cell line has missing value and cell line #45 OVCAR-3 has level "u"

mycellline<-bb[1,]
bb<-bb[-1,]
colnames(bb)<-mycellline

bb[bb==""|bb=="U"]<-NA
bb<-na.omit(bb)

removed = c(which(rowSums(bb=='A')==60), which(rowSums(bb=='H')==60), which(rowSums(bb=='B')==60))
bb = bb[-removed,]
##convert genotype to probability to get trueprobability

bb[bb=="A"]="0"
bb[bb=="H"]="0.5"
bb[bb=="B"]="1"


bb<-apply(bb,2,as.numeric)
##turn dataframe to be matrix
checkbb<-as.matrix(bb)


##get all the positions that provides info
mut1<-which((checkbb[,9]-checkbb[,11])!=0)
mut2<-which((checkbb[,16]-checkbb[,48])!=0)
mut3<-which((checkbb[,16]-checkbb[,50])!=0)
mut4<-which((checkbb[,48]-checkbb[,50])!=0)
mut5<-which((checkbb[,27]-checkbb[,33])!=0)
mut6<-which((checkbb[,27]-checkbb[,34])!=0)
mut7<-which((checkbb[,33]-checkbb[,34])!=0)

##take the union
myvec<-list(mut1,mut2,mut3,mut4,mut5,mut6,mut7)
mymut<-Reduce(union,myvec)

##get the filtered dataset of bb
cc<-bb[mymut,]
bb<-bb[mymut,]

##check na in bb
length(which(is.na(bb)))

##################################################0926:what is new###########################################
##get the cancer type for the 58 cell lines  and we will take care of the new cell lines below
mygroup<-mycellline

for (i in 1:length(mycellline)) {
  mygroup[i]<-gsub(":.*", "", mycellline[i])
  
}

newmygroup<-c()
for (i in 1:length(mygroup)) {
  newmygroup=append(newmygroup,mygroup[,i])
  
}
##set the MDA-N as breast cancer
newmygroup[34]<-"BR"

##substitute with full name of cancer type

mycancertype<-c("Breast cancer","CNS cancer","Colon cancer","Leukemia","Melanoma","Lung cancer"
                ,"Ovarian cancer","Prostate cancer","Renal cancer")
myshort<-unique(newmygroup)

for (i in 1: length(myshort)) {
  newmygroup[which(newmygroup==myshort[i])]<-mycancertype[i]
  
}

##################################################0926:what is new###########################################

##import growth percent for each cell line

growthpercent<-read.csv("growth percent.csv",header  = F)
growthpercent$V1<-as.character(growthpercent$V1)
##get rid of cancer type
for (i in 1:length(mycellline)) {
  mycellline[i]<-gsub(".*:", "", mycellline[i])
  
}
growthpercent[which(growthpercent[,1]=="MDA-MB-231/ATCC"),1]<-"MDA-MB-231"
growthpercent[which(growthpercent[,1]=="MDA-MB-468"),1]<-"MDA-N"


mytemp<-data.frame(V1=c("K-562","CAKI-1"),V2=rep(colMeans(growthpercent[,2,drop=F]),2))
growthpercent<-rbind(growthpercent,mytemp)
growthpercent$V2<-growthpercent$V2/100
##check coverage of cell lines
length(which(mycellline %in% growthpercent$V1  ==T))

### For subsampling cell lines
#bb = bb[,1:6]
#q1 = c(which(rowSums(bb==1)==ncol(bb)), which(rowSums(bb==0.5)==ncol(bb)), which(rowSums(bb==0)==ncol(bb)))
#bb = bb[-q1,]
#mycellline = mycellline[,1:6]

##then we will get true proportion for dmso
dmsoprop<-rep(1/length(mycellline),length(mycellline))

treatprop<-growthpercent$V2/colSums(growthpercent[,2,drop=F])

##we would like to simulate 3 replicated samples for each group with different proportion
simudata <- function(trueprob, trueprop, dp=100, seed){
  set.seed(seed)
  ##get the sum probability of each position of all cell lines
  sumprob <- rowSums(t(trueprop*t(trueprob)))
  
  ##get my simulated data
  freqdata<-data.frame(DP =rep(dp,nrow(trueprob)),
                       mutcount = rbinom(n=rep(1,nrow(trueprob)), size=dp, prob=sumprob))
  return(freqdata=freqdata)
}

k=3
treatfreqdata<-list()
dmsofreqdata<-list()
for (i in 1:k) {
  dmsofreqdata[[i]]<-simudata(bb,dmsoprop,seed=122+i)
  treatfreqdata[[i]]<-simudata(bb,treatprop,seed=122+i)
}

##make it a standard data table
mydmso<-data.frame(matrix(unlist(dmsofreqdata), nrow=nrow(bb), byrow=F),stringsAsFactors=FALSE)
mytrt<-data.frame(matrix(unlist(treatfreqdata), nrow=nrow(bb), byrow=F),stringsAsFactors=FALSE)
mydp<-data.frame(mytrt,mydmso)
colnames(mydp)<-c("DP.treat1","AD.treat1","DP.treat2","AD.treat2","DP.treat3","AD.treat3",
                  "DP.DMSO1","AD.DMSO1","DP.DMSO2","AD.DMSO2","DP.DMSO3","AD.DMSO3")

tt<-list(trueprob=bb,freqdata=mydp,group=c(1,1,1,0,0,0),
           samples = c("treat1", "treat2", "treat3", "DMSO1", "DMSO2", "DMSO3"))

####estimation
source("MCS.R")

library(trust)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(goeveg)

estiresult<-MCS(tt)

myplotdata<-data.frame(esti=estiresult[[2]],true=treatprop/dmsoprop,sd=estiresult[[4]],cancertype=newmygroup)
write.csv(myplotdata,file = "plot of NCI60_0109.csv")

