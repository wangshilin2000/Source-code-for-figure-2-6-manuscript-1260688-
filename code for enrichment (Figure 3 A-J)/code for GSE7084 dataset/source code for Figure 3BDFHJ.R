# GSEA enrichment

rm(list = ls()) 
options(stringsAsFactors = F)

library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputdocument="diffgene7084.txt"         
gmtdocument="c2.cp.kegg.v7.4.symbols.gmt"      


WSL<-read.table(inputdocument, header=T, sep="\t", check.names=F)
logFC<-as.vector(WSL[,2])
names(logFC)=as.vector(WSL[,1])
logFC<-sort(logFC, decreasing=T)

mygmt<-read.gmt(gmtdocument)

BB<-GSEA(logFC, TERM2GENE=mygmt, pvalueCutoff = 1)
BBTab<-as.data.frame(BB)
BBTab<-BBTab[BBTab$p.adjust<0.05,]
write.table(BBTab,file="GSEA7084.txt",sep="\t",quote=F,row.names = F)

ABC="GSEA7084_sellect.txt"
BBTab<-read.table(ABC, header=T, sep="\t", check.names=F)
rownames(BBTab)=BBTab[,1]

termNum<-6      
BBUp=BBTab[BBTab$NES>0,]
if(nrow(BBUp)>=termNum){
  showTerm=row.names(BBUp)[1:termNum]
  gseaplot=gseaplot2(BB, showTerm, base_size=8, title="Enriched in Treat Group")
  pdf(file="Figure3B.pdf", width=8, height=6)
  print(gseaplot)
  dev.off()
}

dev.off()

# GO enrichment

rm(list = ls()) 
options(stringsAsFactors = F)
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(colorspace,stringi,ggplot2,digest,GOplot,clusterProfiler,org.Hs.eg.db,DOSE,enrichplot)


outputfile1="GO_up.txt"
outputfile2= "Figure3D.pdf"

pFilter=0.05     
qFilter=0.05     


biof =read.table("7084_2507upgene_name.txt", header=T, sep="\t", check.names=F)


colnames(biof)[1]="Gene"
genes=as.vector(biof[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        


colorSel="qvalue"
if(qFilter>0.05)
{colorSel="pvalue"}


bio=enrichGO(gene=gene,OrgDb = org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff = 1, ont="all", readable =T)
GORich=as.data.frame(bio)
GORich=GORich[(GORich$pvalue<pFilter & GORich$qvalue<qFilter),]

write.table(GORich,file=outputfile1,sep="\t",quote=F,row.names = F)


showNum=15
if(nrow(GORich)<15){
  showNum=nrow(GORich)
}

pdf(file=outputfile2,width =12,height =10)
Gd=dotplot(bio,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel,label_format = 100) + facet_grid(ONTOLOGY~., scale='free')
print(Gd)
dev.off()

#downgene

rm(list = ls()) 
options(stringsAsFactors = F)

library(pacman)
p_load(colorspace,stringi,ggplot2,digest,GOplot,clusterProfiler,org.Hs.eg.db,DOSE,enrichplot)


outputfile1="GO_down.txt"
outputfile2= "Figure3F.pdf"

pFilter=0.05     
qFilter=0.05 
biof =read.table("7084_2507downgene_name.txt", header=T, sep="\t", check.names=F)


colnames(biof)[1]="Gene"
genes=as.vector(biof[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        


colorSel="qvalue"
if(qFilter>0.05)
{colorSel="pvalue"}


bio=enrichGO(gene=gene,OrgDb = org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff = 1, ont="all", readable =T)
GORich=as.data.frame(bio)
GORich=GORich[(GORich$pvalue<pFilter & GORich$qvalue<qFilter),]

write.table(GORich,file=outputfile1,sep="\t",quote=F,row.names = F)


showNum=15
if(nrow(GORich)<15){
  showNum=nrow(GORich)
}

pdf(file=outputfile2,width =12,height =10)
Gd=dotplot(bio,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel,label_format = 100) + facet_grid(ONTOLOGY~., scale='free')
print(Gd)
dev.off()

#  KEGG enrichment

p_load(clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2, circlize, RColorBrewer, dplyr, ggpubr, ComplexHeatmap)


pFilter <- 1
qFilter <- 1

colorSel <- ifelse(qFilter > 0.05, "pvalue", "qvalue")

rt <- read.table("7084_2507upgene_name.txt", header = TRUE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
rt <- data.frame(genes, entrezID = as.character(entrezIDs))
gene <- entrezIDs[entrezIDs != "NA"]
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) {
  match_idx <- match(strsplit(x, "/")[[1]], as.character(rt$entrezID))
  paste(rt$genes[match_idx], collapse = "/")
}))

KEGG <- KEGG[(KEGG$pvalue < pFilter & KEGG$qvalue < qFilter), ]

write.table(KEGG, file = "KEGG_7084up.txt", sep = "\t", quote = FALSE, row.names = FALSE)
showNum <- min(15, nrow(KEGG))
pdf("Figure3H.pdf", width = 9, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 100, color = colorSel)
dev.off()

# downgene

rt <- read.table("7084_2507downgene_name.txt", header = TRUE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
rt <- data.frame(genes, entrezID = as.character(entrezIDs))
gene <- entrezIDs[entrezIDs != "NA"]
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) {
  match_idx <- match(strsplit(x, "/")[[1]], as.character(rt$entrezID))
  paste(rt$genes[match_idx], collapse = "/")
}))

KEGG <- KEGG[(KEGG$pvalue < pFilter & KEGG$qvalue < qFilter), ]

write.table(KEGG, file = "KEGG_7084down.txt", sep = "\t", quote = FALSE, row.names = FALSE)
showNum <- min(15, nrow(KEGG))
pdf("Figure3J.pdf", width = 9, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 100, color = colorSel)
dev.off()

