
rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,pheatmap,ggplot2,limma)


## heatmap

biof<-read.csv("ht_GSE140947.csv")

biof=as.matrix(biof)
rownames(biof)=biof[,1]
biof=biof[,-1]
biof=matrix(as.numeric(as.matrix(biof)),nrow=nrow(biof),dimnames=list(rownames(biof),colnames(biof)))

gse_N <- c(rep(1,12),rep(2,12))
type <- c(rep("normal",12),rep("AA",12))

design <- model.matrix(~ 0+factor(type))
colnames(design) <- c("AA", "normal")

fit <- lmFit(biof, design)

contrast.matrix <- makeContrasts(AA-normal, levels=design)
bioffit <- contrasts.fit(fit,contrast.matrix)

bioffit <- eBayes(bioffit)

allgene<-topTable(bioffit, coef=1, adjust="fdr",number = "all")

allgene1 <- data.frame(id=row.names(allgene),allgene)

write.table(allgene1,"allgene140947.txt",sep = "\t",quote = F)
colnames(allgene)


FCfilter=1
pfilter=0.05

diffgene<-allgene[abs(allgene$logFC)>=FCfilter & allgene$P.Value<pfilter,]

diffgene <- data.frame(id=row.names(diffgene),diffgene)
fwrite(diffgene,"diffgene140947.txt",sep = "\t",quote = F)

upgene<-allgene[allgene$logFC>=FCfilter & allgene$P.Value<pfilter,]   
upgene <- data.frame(id=row.names(upgene),upgene)
fwrite(upgene,"upgene140947.txt",sep = "\t",quote = F)

downgene<-allgene[allgene$logFC<= (-FCfilter) & allgene$P.Value<pfilter,] 
downgene <- data.frame(id=row.names(downgene),downgene)
fwrite(downgene,"downgene140947.txt",sep = "\t",quote = F)

diff_list=diffgene[,1]

heatmapexp<-biof[diff_list,]


Geo=c(rep("GSE140947",24))
Type=c(rep("normal",12),rep("AA",12))
names(Geo)=colnames(heatmapexp)
annotation=cbind(Geo,Type)
annotation=as.data.frame(annotation)


p=pheatmap(heatmapexp, 
           scale = "row",    
           annotation=annotation, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),   
           cluster_cols =F,    
           fontsize_row=4,   
           fontsize_col=6,    
           show_rownames = T,   
           border=FALSE)  
p
pdf("Figure3K.pdf",10,8)
p
dev.off()

# For GSE7084 dataset

biof<-read.csv("ht_gse7084_2507.csv")

biof=as.matrix(biof)
rownames(biof)=biof[,1]
biof=biof[,-1]
biof=matrix(as.numeric(as.matrix(biof)),nrow=nrow(biof),dimnames=list(rownames(biof),colnames(biof)))

gse_N <- c(rep(1,8),rep(2,7))
type <- c(rep("normal",8),rep("AA",7))

design <- model.matrix(~ 0+factor(type))
colnames(design) <- c("AA", "normal")

fit <- lmFit(biof, design)

contrast.matrix <- makeContrasts(AA-normal, levels=design)
bioffit <- contrasts.fit(fit,contrast.matrix)

bioffit <- eBayes(bioffit)

allgene<-topTable(bioffit, coef=1, adjust="fdr",number = "all")

allgene1 <- data.frame(id=row.names(allgene),allgene)

write.table(allgene1,"allgene7084.txt",sep = "\t",quote = F)
colnames(allgene)


FCfilter=1
pfilter=0.05

diffgene<-allgene[abs(allgene$logFC)>=FCfilter & allgene$P.Value<pfilter,]

diffgene <- data.frame(id=row.names(diffgene),diffgene)
fwrite(diffgene,"diffgene7084.txt",sep = "\t",quote = F)

upgene<-allgene[allgene$logFC>=FCfilter & allgene$P.Value<pfilter,]   
upgene <- data.frame(id=row.names(upgene),upgene)
fwrite(upgene,"upgene7084.txt",sep = "\t",quote = F)

downgene<-allgene[allgene$logFC<= (-FCfilter) & allgene$P.Value<pfilter,] 
downgene <- data.frame(id=row.names(downgene),downgene)
fwrite(downgene,"downgene7084.txt",sep = "\t",quote = F)

diff_list=diffgene[,1]

heatmapexp<-biof[diff_list,]


Geo=c(rep("GSE7084",15))
Type=c(rep("normal",8),rep("AA",7))
names(Geo)=colnames(heatmapexp)
annotation=cbind(Geo,Type)
annotation=as.data.frame(annotation)


p=pheatmap(heatmapexp, 
           scale = "row",    
           annotation=annotation, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),   
           cluster_cols =F,    
           fontsize_row=4,   
           fontsize_col=6,    
           show_rownames = T,   
           border=FALSE)  
p
pdf("Figure3L.pdf",10,8)
p
dev.off()
