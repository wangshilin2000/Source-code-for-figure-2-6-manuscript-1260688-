# lasso

rm(list = ls())
options(stringsAsFactors = F)

if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library(pacman)
p_load(survival,glmnet,ggplot2,ggsci,patchwork,limma)



inputFile="ht_GSE140947.txt"       
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)   


sample=read.table("clinical140947.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)


group <- c(rep("0",12),rep("1",12))

group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])

set.seed(8888)

if(T){
  cvfit = cv.glmnet(x, y,family = "binomial", nlambda=100, alpha=1,nfolds = 10) 
  fit <- glmnet(x,y,family = "binomial")
  cvfit$lambda.min
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)  
}


write.table(geneCoef, file="geneCoef140947.xls", sep="\t", quote=F, row.names=F)
write.table(file="lassoset140947.txt",lassoGene,sep="\t",quote=F,col.names=F,row.names=F) 

pdf("Figure4DE.pdf",height = 5,width = 7)
layout(matrix(c(1,1,2,2), 2, 2, byrow = F))   
plot(fit,xvar = 'lambda')
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")  
dev.off()


# SVM

rm(list = ls()) 
options(stringsAsFactors = F)

if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library(pacman)
p_load(glmnet,tidyverse,VennDiagram,sigFeature,e1071,caret,randomForest,limma)
source('biofsvmRFE.R')   

inputFile="ht_GSE140947.txt" 
biof=read.table(inputFile, header=T, sep="\t", check.names=F)
biof=as.matrix(biof)
rownames(biof)=biof[,1]
exp=biof[,-1]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)


sample=read.table("clinical140947.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]

group <- c(rep("0",12),rep("1",12))
group=as.matrix(as.numeric(group))
rownames(group)=rownames(data)
colnames(group)="Type"
input <- as.data.frame(cbind(group,data))
input$Type=as.factor(input$Type)

set.seed(8888)

if(T){
  svmRFE(input, k = 10, halve.above = 20) 
  nfold = 10
  nrows = nrow(input)
  folds = rep(1:nfold, len=nrows)[sample(nrows)]
  folds = lapply(1:nfold, function(x) which(folds == x))
  results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) 
  top.features = WriteFeatures(results, input, save=F) 
  head(top.features)
  write.csv(top.features,"feature_svm140947.csv")
  featsweep = lapply(1:20, FeatSweep.wrap, results, input) 
  no.info = min(prop.table(table(input[,1])))
  errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
}

pdf("Figure4H.pdf",width = 10,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()

pdf("Figure4I.pdf",width = 10,height = 5)
Plotaccuracy(1-errors,no.info=no.info) 

dev.off()

which.min(errors) 


# Random forest

rm(list = ls())
options(stringsAsFactors = F)
#Sys.setenv("VROOM_CONNETION_SIZE"=131072*6000)

options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

library("pacman")
p_load(randomForest,limma,ggpubr,dplyr,limma)
inputFile="ht_gse140947.txt"      


biof=read.table(inputFile, header=T, sep="\t", check.names=F)
biof=as.matrix(biof)
rownames(biof)=biof[,1]
exp=biof[,-1]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

sample=read.table("clinical140947.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
colnames(data)=gsub("-", "_", colnames(data))

group <- c(rep("normal",12),rep("AA",12))


set.seed(8888)
bdf=randomForest(as.factor(group)~.,
                 data=data,impotance=ture,na.action = na.pass,
                 ntree=20)
bdf
pdf(file="Figure4L.pdf", width=6, height=6)
plot(bdf, main="Random forest", lwd=2)
dev.off()


optionTrees=which.min(bdf$err.rate[,1])
optionTrees
bdf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

importance=importance(x=bdf2)
importance=as.data.frame(importance)
importance$size=rownames(importance)
importancedata=importance[,c(2,1)]
names(importancedata)=c("Gene","importance")

importancedata= arrange(importancedata,desc(importance))
bfsci=importancedata[1:20,]

if(T){
  p=ggdotchart(bfsci, x = "Gene", y = "importance",
               color = "importance", 
               sorting = "descending",                      
               add = "segments",                            
               add.params = list(color = "lightgray", size = 2), 
               dot.size = 6,                        
               font.label = list(color = "white", size = 9,
                                 vjust = 0.5),               
               ggtheme = theme_bw()         ,               
               rotate=TRUE                                       )
  p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
    gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +
    grids()   
}


pdf(file="Figure4N.pdf", width=6, height=6)
print(p1)
dev.off()

bdfGenes=arrange(importancedata,desc(importance))
write.table(bdfGenes, file="bdfGenes140947.xls", sep="\t", quote=F, col.names=T, row.names=F)
