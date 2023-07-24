
# Analysis of immune cell infiltration ( use GSE7084 database )
rm(list = ls()) 
options(stringsAsFactors = F)


if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(e1071,parallel,devtools,preprocessCore,limma,impute)


source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "gse7084.txt", perm=100, QN=T)
write.table(results,"CIBERSORT-Results.txt",sep="\t",quote = F,row.names = T)


rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,RColorBrewer,ggpubr)


biof<- read.table("CIBERSORT-Results.txt",header=T,sep="\t",check.names=F,row.names=1)  
ncol(biof)
biof2=biof %>% 
  filter(`P-value`<=0.05) %>%    
  select(-c(23:25))    

biof2=cbind(id=row.names(biof2), biof2)
row.names(biof2)=NULL


cln=read.table("clinical7084.txt",header = T,sep = "\t")
WSL=merge(cln,biof2,by="id")

mydata<-WSL %>% as.data.frame() %>% melt(variable.name="imc",value.name="pecent")
mypalette <- colorRampPalette(brewer.pal(8,"Set3"))


pdf("Figure6A.pdf",15,12)
ggplot(mydata,aes(x=id,y=pecent,fill=imc))+  
  geom_point(aes(color =status,y=-Inf),shape = 2,size = 4)+    
  geom_bar(stat = 'identity')+ xlab("") +
  ylab("Estiamted Relative Percentage") +      
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+ 
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())

dev.off()


rm(list = ls()) 
options(stringsAsFactors = F)
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,RColorBrewer,ggpubr)


biof<- read.table("CIBERSORT-Results.txt",header=T,sep="\t",check.names=F,row.names=1)
ncol(biof)
biof2=biof %>% 
  filter(`P-value`<=0.05) %>% 
  select(-c(23:25))

if(T){
  mypalette <- colorRampPalette(brewer.pal(8,"Set3"))
  biofdata <- biof2 %>% as.data.frame() %>%
    rownames_to_column("Sample") %>% 
    gather(key = Cell_type,value = Proportion,-Sample)
  a = biofdata %>% 
    group_by(Cell_type) %>%     
    summarise(m = median(Proportion)) %>%    
    arrange(desc(m)) %>%    
    pull(Cell_type)
  
  biofdata$Cell_type = factor(biofdata$Cell_type,levels = a)   
  
  ggplot(biofdata,aes(Cell_type,Proportion,fill = Cell_type)) + 
    geom_boxplot(outlier.shape = 21,color = "black",alpha=0.5) + 
    theme_bw() + 
    labs(x = "Cell_type", y = "Estiamted Relative Percentage") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") + 
    scale_fill_manual(values = mypalette(22))
}

# export figure as Figure6B



rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,corrplot,vioplot)


biof<- read.table("CIBERSORT-Results.txt",header=T,sep="\t",check.names=F)
ncol(biof)
colnames(biof)[1]="id"
ncol(biof)

biof2=biof %>% 
  filter(`P-value`<=0.05) %>% 
  select(-c(23:25))   

write.table(biof2,"biof.txt",sep="\t",quote = F,row.names = F)
# Only one row name is retained when the row name is repeated
mydata=read.table("biof.txt",sep="\t",header=T,row.names=1,check.names=F)


pdf("Figure6C.pdf",12,10)              
corrplot(corr=cor(mydata),
         method="circle",
         type="upper",
         addCoefasPercent=T,
         order = "AOE", 
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,   
         col=colorRampPalette(c("navy", "white", "darkred"))(100)  
)
dev.off()


group1=3 
group2=7 
pdf("Figure6D.pdf",10,10)              
par(las=1,mar=c(10,6,3,3))
plot(c(1:ncol(mydata)),c(1:ncol(mydata)),
     xlim=c(0,63),ylim=c(min(mydata),max(mydata)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="green",
     xaxt="n")

for(i in 1:ncol(mydata)){
  group1d=mydata[1:group1,i]
  group2d=mydata[(group1+1):(group1+group2),i]
  vioplot(group1d,at=3*(i-1),lty=1,add = T,col = 'navy')
  vioplot(group2d,at=3*(i-1)+1,lty=1,add = T,col = 'firebrick3')
  WT=wilcox.test(group1d,group2d)
  p=round(WT$p.value,4)
  mline=max(c(group1d,group2d))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mline,mline))
  text(x=3*(i-1)+0.5, y=mline+0.02,labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.6)
  text(seq(1,64,3),-0.09,xpd = NA,labels=colnames(mydata),cex = 0.8,srt = 65,pos = 1)
}
dev.off()



rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,corrplot,vioplot)

biof<- read.table("gse7084.txt",header=T,sep="\t",check.names=F,row.names=1)
biof=t(biof)
sample=read.table("clinical7084.txt",sep="\t",header=T,check.names=F,row.names = 1)

biof=biof[rownames(sample),]  
x=as.matrix(biof)

biof2<-merge(sample, x,by = "row.names", all = T)
colnames(biof2)[1]="id"
colnames(biof2)[2]="status"
veengene<-read.table("hubgene.txt",header=F,sep="\t",check.names = F)   
veenExp=biof2[,c("id",as.vector(veengene[,1]))]

CIB_data=read.table("CIBERSORT-Results.txt",header = T,sep = "\t",check.names = F)
CIB_data=CIB_data[,-c(24:26)]
CIB_data=cbind(id=row.names(CIB_data), CIB_data)
row.names(CIB_data)=NULL


mydata=merge(veenExp,CIB_data,by="id")

ncol(veenExp) 

bfsci=data.frame()
for(b in colnames(mydata[,2:ncol(veenExp)])){    
  for(a in colnames(mydata[,(ncol(veenExp)+1):ncol(mydata)])){     
    v1=as.numeric(mydata[,b])
    v2=as.numeric(mydata[,a])
    mycor=cor.test(v1,v2)
    cor=mycor$estimate
    p=mycor$p.value
    bfsci=rbind(bfsci,cbind(gene=b,imc=a,cor,p))
    
  }
}

write.table(bfsci,"Venn_imc.txt",row.names = F,sep = "\t",quote = F)




rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,patchwork,ggplot2)


mydata=read.table("Venn_imc.txt",header = T,sep = "\t")
mydata=mydata %>%  drop_na()  
table(mydata$gene)    


venn_gene="JPH2"  # Change the name of each gene

if(T){
  
  mydata=mydata[mydata$gene==venn_gene,]
  mydata$p2 <- cut(mydata$p,
                   breaks = c(0, 0.05,1), 
                   labels = c("< 0.05","> 0.05"),right=F)
  p3=mydata[,"p"]
  p3=ifelse(p3<0.001, "<0.001", sprintf("%.4f", p3)) 
  
  
  P1=ggplot(mydata, aes(x=cor, y=imc,color=p),colour=NA) +
    geom_segment(aes(x=0, xend=cor, y=imc, yend=imc))+
    geom_point(aes(size=abs(cor)))+
    scale_colour_gradient2(low="#009900",mid="#006600",high = "#FF3300", midpoint = 0)+
    theme_linedraw()+
    labs(y=NULL,title =venn_gene,
         x="Correlation coefficient",
         size = "abs(cor)")
  P1
  P2=ggplot()+geom_text(mydata,mapping = aes(x = 0, y =imc, 
                                             color = p2, 
                                             label = p3))+
    scale_color_manual(name="",
                       values = c("darkred", "grey85"))+
    theme_void()+
    guides(color=F)
  P1|P2
}
ggsave(paste0(venn_gene,"_","Correlation_coefficient",".pdf"),width = 6,height = 6)


