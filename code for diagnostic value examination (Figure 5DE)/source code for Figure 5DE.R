# ROC

rm(list = ls()) 
options(stringsAsFactors = F)

if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library(pacman)
p_load(data.table,pROC,limma)

inputFile="gse47472.txt"     

biof=read.table(inputFile, header=T, sep="\t", check.names=F)
biof=as.matrix(biof)
rownames(biof)=biof[,1]
exp=biof[,-1]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

sample=read.table("clinical47472.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)


group <- c(rep("0",14),rep("1",8)) 
group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])


biof2<-merge(y, x,by = "row.names", all = T)
colnames(biof2)[1]="id"
colnames(biof2)[2]="status"
Veendata<-read.table("hubgene.txt",header=F,sep="\t",check.names = F)
veenExp=biof2[,c("id","status",as.vector(Veendata[,1]))]


for(wenyun in colnames(veenExp[,3:ncol(veenExp)])){
  myroc=roc(veenExp$status, veenExp[,wenyun])
  myci=ci.auc(myroc, method="bootstrap")
  mycinum=as.numeric(myci)
  pdf(file=paste0(wenyun,"_","ROC",".pdf"))   
  plot(myroc, print.auc=T, col="#00AFBB", legacy.axes=T, main=wenyun)   
  text(0.40, 0.45, paste0("95% CI: ",sprintf("%.03f",mycinum[1]),"-",sprintf("%.03f",mycinum[3])), col="#00AFBB")
  dev.off()
}





# nomograph

library(haven)
data <- read.csv("ht_hubgene.csv")
data<-as.data.frame(data)
data<-na.omit(data)

View(data)
names(data)
str(data)

library(rms)
dd=datadist(data)
options(datadist="dd")

form<-as.formula("Sur_status~ALPL+BLNK+JPH2+HLA.DQB1+HLA.DMA+HLA.DRA")

model <- lrm(form, data=data) 

fun1 <- function(x) plogis(x-model$coef[1]+model$coef[1])
fun2 <- function(x) plogis(x-model$coef[1]+model$coef[2])

f <- Newlabels(model, c(Sur_status='group'))  

g <- nomogram(f, fun=list('Prob Y=1'=plogis, 'Prob Y=0'=fun2), 
              fun.at=c(.01,.05,seq(.1,.9,by=.1),.95,.99))
plot(g)

# export nomograph as Figure 5E