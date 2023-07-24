# miRNA_mRNA

rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,ggvenn)

venn_gene=read.table("hubgene.txt",header = F,sep = "\t",check.names = F)
venn_gene=unique(venn_gene)
colnames(venn_gene)[1]= "gene"


biof1=fread("Starbase.txt",header = T,sep = "\t",check.names = F)
biof2=fread("miRTarBase.txt",header = T,sep = "\t",check.names = F)
biof3=fread("TargetScan.txt",header = T,sep = "\t",check.names = F)


conb1=merge(venn_gene,biof1,by="gene")
conb2=merge(venn_gene,biof2,by="gene")
conb3=merge(venn_gene,biof3,by="gene")

write.table(conb1,"miRNA_Starbase.txt",sep = "\t",quote = F,row.names = F)
write.table(conb2,"miRNA_miRTarBase.txt",sep = "\t",quote = F,row.names = F)
write.table(conb3,"miRNA_TargetScan.txt",sep = "\t",quote = F,row.names = F)


Starbase_miRNA_list = conb1[,2] 
Starbase_miRNA_list=unique(Starbase_miRNA_list)

miRTarBase_miRNA_list = conb2[,2] 
miRTarBase_miRNA_list=unique(miRTarBase_miRNA_list)

TargetScan_miRNA_list = conb3[,2] 
TargetScan_miRNA_list=unique(TargetScan_miRNA_list)

datalist <- list("Starbase"=Starbase_miRNA_list,
                 "miRTarBase"=miRTarBase_miRNA_list,
                 "TargetScan"=TargetScan_miRNA_list)


opar <- par(family = "Roboto Condensed")
biocolor<- c("blue", "yellow", "green", "red")
ggvenn(datalist,
       fill_color=biocolor,
       fill_alpha = .7,
       stroke_linetype = "longdash",
       set_name_size = 5,
       text_size=5) 

miR_con1 <- inner_join(conb1,conb2,by ="miRNA")
miR_con2 <- inner_join(miR_con1,conb3,by ="miRNA")
dim(miR_con2)


fwrite(data.table(unique(miR_con2[,2])), 
       "miRNA_veen.txt", 
       col.names=F,
       row.names = F, 
       sep = "\t", 
       quote = F)


miRNA_veen = fread("miRNA_veen.txt",header = F,sep = "\t",check.names = F)  
colnames(miRNA_veen) ="miRNA"

Sveen = fread("miRNA_Starbase.txt",header = T,sep = "\t",check.names = F)  

miRNA_mRNA= inner_join(miRNA_veen,Sveen,by="miRNA")

dim(miRNA_mRNA)

miRNA_mRNA=unique(miRNA_mRNA)
dim(miRNA_mRNA)

fwrite(miRNA_mRNA, 
       "miRNA_mRNA.txt",     
       col.names=T,
       row.names = F, 
       sep = "\t", 
       quote = F)



# TF - miRNA - mRNA network

rm(list = ls())
options(stringsAsFactors = F)
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(stringi,data.table, ggplot2,tidyverse,tidyr)


biof = read.table("TRANSFAC_and_JASPAR_PWMs_table.txt",header=T,sep="\t", stringsAsFactors=FALSE,check.names=F)
colnames(biof)
biof=biof[,c("Term","P-value","Genes")]    
biof$Term

if(T) {
  biof2 <- biof %>% 
    separate(Term,c("TF","type"), sep = " \\(")%>%    
    filter(type=="human)")%>%
    filter(`P-value` <= 0.5)%>%
    select(c(1,3,4))
  
  biodate <-biof2 %>% as_tibble() %>% 
    separate_rows(Genes, sep = ";")     
  colnames(biodate) }

biodate= biodate %>% select(c(1,3))


fwrite(biodate, 
       "TF_mRNA.txt",     
       col.names=T,
       row.names = F, 
       sep = "\t", 
       quote = F)


if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(stringi,data.table, ggplot2,tidyverse,tidyr)

biof1 = read.table("miRNA_mRNA.txt",header=T,sep="\t", stringsAsFactors=FALSE,check.names=F)
colnames(biof1)=c("nod1","nod2")     
biof2 = read.table("TF_mRNA.txt",header=T,sep="\t", stringsAsFactors=FALSE,check.names=F)
colnames(biof2)=c("nod1","nod2")
biof=rbind(biof1,biof2)


miRNAtype=data.frame(gene=unique(biof1[,1]),type="miRNA")
TFtype=data.frame(gene=unique(biof2[,1]),type="TF")
mRNAtype=data.frame(gene=unique(biof[,2]),type="mRNA")

type=rbind(miRNAtype,TFtype,mRNAtype)


miRNAlist=data.frame(miRNAtype[,1])
TFlist=data.frame(TFtype[,1])
mRNAlist=data.frame(mRNAtype[,1])


fwrite(miRNAlist, "miRNAlist.txt",  col.names=F, row.names = F,  sep = "\t", quote = F)
fwrite(TFlist, "TFlist.txt",  col.names=F, row.names = F,  sep = "\t", quote = F)
fwrite(mRNAlist, "mRNAlist.txt",  col.names=F, row.names = F,  sep = "\t", quote = F)


fwrite(biof, 
       "miRNA_TF_mRNA.txt", 
       col.names=T,
       row.names = F, 
       sep = "\t", 
       quote = F)


fwrite(type, 
       "type.txt", 
       col.names=T,
       row.names = , 
       sep = "\t", 
       quote = F)
