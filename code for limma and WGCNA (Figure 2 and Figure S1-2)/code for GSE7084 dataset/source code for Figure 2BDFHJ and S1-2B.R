
# Limma

rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
library(pacman)
p_load(tidyverse,data.table,pheatmap,ggplot2,limma)

biof<-read.csv("gse7084_2507.csv")

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


# Volcano plot
xmax<-max(allgene$logFC)
ymax<-max(-log10(allgene$P.Value))
xmax
ymax

allgene$sig = as.factor(ifelse(allgene$P.Value < pfilter & abs(allgene$logFC) > FCfilter,   
                               ifelse(allgene$logFC > FCfilter,"Up", "Down"), "Not"))    

p1=ggplot(allgene,aes(logFC,-log10(P.Value)))

p2=p1+ geom_point(aes(color =sig))+xlim(-xmax,xmax) + ylim(0,ymax)+    
  labs(title="Volcano plot",x="log2FC", y="-log10(P.Value)")+    
  theme(plot.title=element_text(hjust=1))+
  theme_classic() +                                                   
  geom_vline(xintercept = c(-FCfilter,FCfilter), linetype ="dashed") +   
  geom_hline(yintercept = -log10(pfilter), linetype ="dashed") +    
  scale_color_manual(values =c("#00AFBB","grey", "#FC4E07"))       


p2
pdf("Figure2B.pdf")
p2
dev.off()

dev.off()

# heatmap

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
pdf("Figure2D.pdf",10,8)
p
dev.off()

dev.off()

# WGCNA

rm(list = ls()) 
options(stringsAsFactors = F)

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
install.packages("WGCNA")
BiocManager::install("GO.db")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
library(pacman)
library(WGCNA)
p_load(tidyverse,data.table,sva,WGCNA)


inputdata1="gse7084_2507.txt" 


biof=fread(inputdata1,sep="\t",header=T,check.names=F,quote="!",data.table = F)
biof=as.matrix(biof)
rownames(biof)=biof[,1]
biof=biof[,-1]
biof=matrix(as.numeric(as.matrix(biof)),nrow=nrow(biof),dimnames=list(rownames(biof),colnames(biof)))

biofvar=apply(biof,1,var)

biof=biof[which(biofvar>quantile(biofvar, probs = seq(0, 1, 0.25))[4]),]   
datExpr = t(biof)   
dim(datExpr)   


gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average")

plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")

cutset = 140
abline(h=cutset,col="red")
clust = cutreeStatic(sampleTree, cutHeight = cutset, minSize = 10)

table(clust)

ks = (clust==1)
ks

datExpr = datExpr[ks, ]


powers = c(c(1:10), seq(from = 12, to=20, by=2))
SoftTh= pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# export Clustering dendrogram as Figure S1B

pdf("FigureS2B.pdf",8,6)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(SoftTh$fitIndices[,1], -sign(SoftTh$fitIndices[,3])*SoftTh$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(SoftTh$fitIndices[,1], -sign(SoftTh$fitIndices[,3])*SoftTh$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.8,col="red")
plot(SoftTh$fitIndices[,1], SoftTh$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(SoftTh$fitIndices[,1], SoftTh$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

dev.off()

SoftTh$powerEstimate
mypower = SoftTh$powerEstimate
mypower


adjacency = adjacency(datExpr, power = mypower)

TOM = TOMsimilarity(adjacency)

dissTOM = 1-TOM
hierTOM = hclust(as.dist(dissTOM),method="average")


minModuleSize = 100    
geneTree = hclust(as.dist(dissTOM), method = "average")

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = F,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors(GEO)")



MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2
abline(h= MEDissThres, col = "red")

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# export Cluster Dendrogram as Figure 2F

clinical=read.table("clinical7084_2507.txt",header = T,sep = "\t",row.names = 1)
rownames(clinical)
traitrr = match(rownames(datExpr), rownames(clinical))
datcli = clinical[traitrr, ]


moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleTraitCor = cor(MEs, datcli, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

pdf("Figure2H.pdf",8,8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

dev.off()

modNames = substring(names(MEs), 3)
traitNames = names(clinical)


MMcor= as.data.frame(cor(datExpr, MEs, use = "p"))
MMP = as.data.frame(corPvalueStudent(as.matrix(MMcor), 
                                     nSamples))
names(MMcor) = paste("MM", modNames, sep="")
names(MMP) = paste("p.MM", modNames, sep="")


GScor = as.data.frame(cor(datExpr,datcli, use = "p"))
GSP = as.data.frame(corPvalueStudent(as.matrix(GScor), 
                                     nSamples))
names(GScor) = paste("GS.", traitNames, sep="");
names(GSP) = paste("p.GS.", traitNames, sep="");

MMGS<-cbind(MMcor, MMP,GScor, GSP)

write.table(MMGS, "MMGS7084.txt", sep = "\t", quote = F)


module="turquoise"    
column = match(module, modNames)
moduleGenes = moduleColors==module

trait="AA"    
traitColumn=match(trait,traitNames)


par(mfrow = c(1,1))
verboseScatterplot(abs(MMcor[moduleGenes, column]),
                   abs(GScor[moduleGenes, traitColumn]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ",trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

a=which(merge$colors=="turquoise")  
b=datExpr[,a]
c=colnames(b)
write.table(c,"WGCNAgene7084.txt",quote = F,row.names = F,col.names = F)

# export gene significance plot as Figure 2J
