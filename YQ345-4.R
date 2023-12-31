


setwd("")


###
##############################GEO数据处理#########

library(tidyr)
library(tibble)
library(dplyr)
library(data.table)
library(GEOquery)

matrix<-read.delim2('GSE95233_series_matrix.txt',sep = '\t',check.names = F,header = T)
anno <-data.table::fread("GSE95233_family.soft",skip ="ID",header = T)
# write.csv(anno,"anno.csv")
View(anno)
anno$`Gene Symbol`

# anno <- read.csv("anno.csv",check.names = F)
#
# symb = merge(matrix,anno,by ="ID_REF")
# write.csv(symb,"GSE20680_expr.csv")


probe2symbol <- anno %>%
  
  dplyr::select("ID","Gene Symbol") %>% dplyr::rename(probeset = "ID",symbol="Gene Symbol") %>%
  
  filter(symbol!= "") %>%#去重
  tidyr::separate_rows( `symbol`,sep="///")


# fix(matrix)
# fix(probe2symbol)
# AB = merge(matrix,probe2symbol,by = "ID")
# AB = AB[!duplicated(AB$ID),]
# write.csv(AB,"GSE63492_expr.csv")
matrix = matrix[!duplicated(matrix$ID_REF),]
rownames(matrix)=matrix[,1]
matrix=matrix[,-1]
exprSet <- matrix %>% as.data.frame() %>%
  rownames_to_column(var="probeset") %>%
  #合并探针的信息
  inner_join(probe2symbol,by="probeset") %>%
  #去掉多余信息
  dplyr::select(-probeset) %>%
  #重新排列
  dplyr::select(symbol,everything())%>%
  #求出平均数(这边的点号代表上一步产出的数据)
  # mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% #报错的话就从duplicated开始跑
  #去除symbol中的NA
  filter(symbol != "NA") %>%
  #把表达量的平均值按从大到小排序
  # arrange(desc(rowMean)) %>%
  # symbol留下第一个
  distinct(symbol,.keep_all = T) %>%
  #反向选择去除rowMean这一列
  dplyr::select(-rowMean) %>%
  # 列名变成行名
  column_to_rownames(var = "symbol")

duplicated(exprSet$symbol)#看重复
exprSet1 = exprSet[!duplicated(exprSet$symbol),]
rownames(exprSet1) = exprSet1[,1]
exprSet1=exprSet1[,-1]
out<-exprSet1
out=na.omit(out)


write.csv(out,"GSE95233_expr.csv",row.names = T)#输出csv文件

a <- read.csv(file = "GSE65682_expr.csv",check.names = F)
a=a[!duplicated(a$id),]
rownames(a)=a[,1]
a=a[,-1]
b <- read.csv("GSE65682-sample.csv",check.names = F)
id = b$ID
ab = a[,id]
write.table(cbind(Symbol=rownames(ab),ab),'Estimate.txt',sep = '\t',quote = F,row.names = F)

write.csv(ab,'GSE65682_expr.csv')


gene = "ELANE"
ab = out[gene,]
ab = t(ab)
write.csv(ab,"GSE95233-ELANE.csv")
#############################1. CIBERSORT############

library("e1071")
library("preprocessCore")
library("GSVA")
library("estimate")
library("ggpubr")
library("tidyr")
library("ggplot2")
library("pheatmap")

source("CIBERSORT.R")


#CIBERSORT肿瘤免疫微环境评估，样本行名基因，列名样本

results=CIBERSORT("LM22.txt", "Estimate.txt", perm=100, QN=F)

#记得筛选p<0.05的样本

#barplot
output <- read.table('CIBERSORT-Results.txt',sep='\t',check.names = F,header = T)

output <- output[output$`P-value`<0.05,]
write.csv(output,"CIBERSORT(过滤后).csv")

data=output[,-c(24,25,26)]
rownames(data) <- data$Mixture
data <- data[,-1]
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)

pdf('CIBERSORT免疫细胞比例堆积条形图.pdf',height=8,width=12)

# png('CIBERSORT免疫细胞比例heatmap.png',height=700,width=1200)
par(las=1,mar=c(8,8,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=90,xpd=T);text(a1,-0.04,colnames(data),adj=1,cex=1);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),
       col=col,pch=10,bty="n",cex=1)
dev.off()

#横轴不写样本名

data <- read.csv("CIBERSORT(过滤后).csv",check.names=F,row.names=1)
data=data[,-c(23,24,25)]
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)
pdf('CIBERSORT-免疫细胞比例.pdf',height=10,width=18)
#png('infi_heatmap.png',height=800,width=1200)
par(las=1,mar=c(8,8,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F,tick = F)
par(xpd=T);text(100,0,'TCGA-LUAD Patients (N = 431, High-risk = 212, Low-risk = 219)',adj=0.9,cex=1.7);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),
       col=col,pch=15,bty="n",cex=1.2)
dev.off()

###################################1. CIBERSORT画图新式########

#调包日常
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot","ggpubr","bslib","ggthemes")
# BiocManager::install("ggthemes")
library(ggthemes)
library(dplyr)
library(tidyr)
library(tibble)
# Read in results 
####输入文件不用删除P值等倒数三列
cibersort_raw <- read.csv("CIBERSORT(过滤后).csv",row.names = 1,check.names = F)

dd1 <- cibersort_raw %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = 2:23,
               names_to = "CellType",
               values_to = "Composition")

# View(dd1)
plot.info <- dd1[,c(5,1,6)]          


#####箱线图
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition") +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))        
ggsave(filename = "CIBERSORT免疫细胞比例箱线图.pdf",width = 11,height = 8)


####堆积条形图

ggbarplot(
  plot.info,
  x = "sample",
  y = "Composition",
  size = 0,
  fill = "CellType",
  color = "CellType"
  
) +
  theme_base() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 1
    ),
    legend.position = "bottom"
  )
ggsave(filename = "CIBERSORT免疫细胞比例堆积条形图.pdf",width = 14,height = 8)


##################################1. 检测哪些细胞是差异的####

a <- read.csv("CIBERSORT(过滤后).csv",check.names=F)
a = a[,-c(24:26)]
b <- read.csv("GSE65682-sample.csv",check.names = F)
ab = merge(a,b,by="ID")
write.csv(ab,"CIBERSORT-testinput.csv")

ssgseaOut<-read.csv('CIBERSORT-testinput.csv',check.names = F,row.names = 1)##直接用上一步输出文件即可
data= ssgseaOut[,-c(23)]
data$Type<- factor(c(rep('high',760),rep('low',42)),levels = c('high','low'))
data[,1:(ncol(data)-1)] <- apply(data[,1:(ncol(data)-1)],2,function(x){log2(x+1)})
#t检验
tdata <- data.frame()
for (i in colnames(data)[1:(ncol(data)-1)]){
  rt = data[,c(i,'Type')]
  colnames(rt) <- c('immunocyte','Type')
  test <- wilcox.test(immunocyte~Type,rt)
  tdata <- rbind(tdata,data.frame(Gene = i,Pvalue = test$p.value))
}
write.csv(tdata, "CIBERSORT差异免疫细胞检验结果.csv")


#################################1. 差异免疫细胞的箱线图##########
library(tidyr)
library(ggplot2)
library(ggpubr)


dbox <- read.csv("CIBERSORT-testinput.csv",row.names = 1,check.names = F)
dbox1 = dbox[,-c(23)]
dbox1$group <- c(rep('Sepsis',760),rep('Healthy',42))
dbox1  <- gather(dbox1,immune,Score,1:(ncol(dbox1)-1))
View(dbox1)

mycol <- c("turquoise","coral1")



pdf('CIBERSORT差异免疫细胞箱线图.pdf',width = 8,height = 6)

p <- ggboxplot(dbox1,x='immune',y='Score',fill ='group',
               ylab = 'Proportion',
               xlab ='',palette = mycol,
               # add ='jitter',
               size =0.4)+
  rotate_x_text(45)+
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        axis.line = element_line(size=0.3))+
  stat_compare_means(size = 3.5,aes(group=group),
                     label = 'p.signif',method = 'wilcox.test')

p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
  font('x.text',face = 'bold')+font('y.text',face = 'bold')+
  font('legend.title',face = 'bold')+font('legend.text',face = 'bold')

dev.off()

################################1. 中性粒细胞的比较箱线图######

data <- read.csv("CIBERSORT-testinput.csv",
                 check.names = F,row.names=1)
dbox1=data
group <-  c(rep('Sepsis',760),rep('Healthy',42))

length(group) = dim(dbox1)[[1]]
dbox1$group <- group
table(dbox1$group)
my_comparition <- list(c("Sepsis","Healthy"))

mycol <- c("coral1","turquoise")
p <- dbox1 %>%
  ggboxplot(x= "group",y = c(colnames(dbox1)[22]),
            fill = "group",combine = T,
            palette = mycol,
            ylab = "Proportion")

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition)
ggsave(file="中性粒细胞的箱线图.pdf",width = 5,height = 5)

####################1. WGCNA筛选中性粒细胞相关模块############
library(WGCNA)
library(tidyr)


enableWGCNAThreads()   #多线程工作

dataExpr <- read.csv('GSE65682_expr.csv',check.names = F,row.names = 1) %>%
  t() %>% 
  as.data.frame()

sample <- read.csv("GSE65682-sample.csv", row.names = 1, check.names = F)
sample = sample[, c("Sepsis", "Healthy", "Neutrophils")]
# 转置常规表达量，保存dataExpr.Rdata
# 检查所有基因和样本的缺失值是否足够低。
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples) > 0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
#再过滤一遍
dataExpr <- dataExpr[, colMeans(dataExpr) > 1]
save(dataExpr,file = 'dataExpr.Rdata')
sampleTree = hclust(dist(dataExpr), method = "complete")

pdf(file = "Fig1.Sample_cluster.pdf", width = 12, height = 8)
par(mar = c(0,4,2,0), xpd = F)
plot(sampleTree, main = "Sample clustering to detect outliers",cex.axis = 0.5, cex.main = 1, sub = "", xlab = "")
abline(h = 160, col = "red")
dev.off()

##############所以这一步不一定能够做，剪切高度问题,这个根据实际设置后可用
### Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
keepSamples = (clust==1)
# dataExpr = dataExpr[keepSamples,]

# 将疾病与健康作为性状，数字从0-9
dataTraits <- sample
dataTraits <- dataTraits[rownames(dataExpr), ]
# rownames(dataTraits) <- rownames(dataExpr)
save(dataTraits, file = 'dataTraits.Rdata')
head(dataTraits)
sampleTree2 = hclust(dist(dataExpr), method = "complete")
traitColors = numbers2colors(dataTraits, signed = FALSE)
pdf("Fig2.cluster.dendrogram.pdf",width=14,height=8)
# par(pin=c(14, 8), mai=c(10,500,5,5), mar = c(10,500,5,5))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(dataTraits), cex.dendroLabels=1,
                    marAll = c(1, 7, 3, 1),
                    # cex.colorLabels = 0.6,
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# save(dataExpr,file='dataExpr.Rdata')
##power值散点图
powers =seq(from = 1, to=20, by=1)  #幂指数范围1:20
sft = pickSoftThreshold(dataExpr,  RsquaredCut = 0.85, # networkType = "signed",
                        powerVector = powers, verbose = 5)
sft$powerEstimate
# save(sft,file='sft.Rdata')
pdf('Fig3.Soft.Threshold.pdf',width = 10,height = 6)
par(mfrow = c(1,2))
cex1 = 0.85
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.87,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


## 检验选定的β值下记忆网络是否逼近 scale free 
# ADJ1_cor <- abs(WGCNA::cor( dataExpr,use = "p" ))^softPower
# 基因少（<5000）的时候使用下面的代码：
# k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))
# 基因多的时候使用下面的代码：

k <- softConnectivity(datE=dataExpr, power= 4) 
sizeGrWindow(10, 10)
pdf("scaleFreePlot.pdf", height = 6, width = 15)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

###邻接矩阵转换
softPower = sft$powerEstimate #最佳power值
#服务器跑
adjacency = adjacency(dataExpr, power = softPower)
softPower
###TOM矩阵,服务器跑
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
# save(dissTOM,file='dissTOM.Rdata')
#load('dissTOM.Rdata')
###基因聚类，服务器跑
geneTree = hclust(as.dist(dissTOM), method = "average")
save(geneTree, file='geneTree.Rdata')
#load('geneTree.Rdata')
pdf(file="Gene.cluster.pdf",width=6,height=5)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

###动态剪切模块识识别
#服务器跑
minModuleSize =  100  #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
save(dynamicMods,file='dynamicMods.Rdata')
# load('dynamicMods.Rdata')
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="Fig4.DynamicTree.pdf",width=6,height=5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


############################1.WGCNA本地部分#########
library(WGCNA)

load("dataExpr.Rdata")
load('dynamicMods.Rdata')
dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# save(MEList, file = "MEList.RDdata")
#如果不合并，运行下两行
# moduleColors = MEList$validColors
# table(moduleColors)
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf('Fig5.merge.cluster.pdf',width = 9,height = 5)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2 #剪切高度可修改
abline(h=MEDissThres, col = "red")
dev.off()

###相似模块合并
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
load('geneTree.Rdata')
pdf(file="Fig6.DynamicTree.pdf",width=6,height=5)
plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
# save(moduleColors,file = 'moduleColors.Rdata')

table(moduleColors)
colorOrder = c("grey", standardColors(30))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

###模块与性状数据热图
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
load("dataTraits.Rdata")
moduleTraitCor = cor(MEs, dataTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
write.table(moduleTraitPvalue,'moduleTraitPvalue.txt',sep = '\t',quote = F,row.names = T)
write.table(moduleTraitCor,'moduleTraitCor.txt',sep = '\t',quote = F,row.names = T)

pdf(file="Fig7.Module_trait.pdf",width=12,height=12)
par(mar = c(10, 10, 3, 3))

#必须至少是两列，如果是一列就如下处理
# dataTraits1 <- cbind(dataTraits, dataTraits)
# labeledHeatmap(Matrix = cbind(moduleTraitCor,moduleTraitCor),
#                xLabels = colnames(dataTraits1),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(500),
#                textMatrix = cbind(textMatrix,textMatrix),
#                setStdMargins = FALSE,
#                cex.text =0.8,
#                main = paste("Module-trait relationships"))
# dev.off()

#大于两列：
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(500),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =0.8,
               main = paste("Module-trait relationships"))
dev.off()

## 所有基因GS绝对值的平均数来表征该module与表型之间的相关性
dir.create("across_modules")
for (i in 1:ncol(dataTraits)) {
  which.trait <- colnames(dataTraits)[i]
  
  moduleTraitCor[, which.trait]
  
  y <- dataTraits[, which.trait]
  
  GS <- as.numeric(cor(y ,dataExpr, use="p"))
  
  GeneSignificance <- abs(GS)
  
  ModuleSignificance <- tapply(
    
    GeneSignificance,
    
    moduleColors, mean, na.rm=T)
  pdf(file=paste0("./across_modules/", which.trait, "-", "Gene significance across modules.pdf"),width=12,height=5)
  plotModuleSignificance(GeneSignificance, moduleColors)
  dev.off()
}


###计算MM和GS值
#modNames = substring(names(MEs), 3)
#modNames <- c('brown','greenyellow','yellow','cyan','grey60','pink')
modNames <- unique(moduleColors)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(MMPvalue) = paste("p.", colnames(MMPvalue), sep="")

traitNames=colnames(dataTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, dataTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

# 基因与性状的相关性
# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析

pdf('Fig8.Gene.trait_MEcyan-1.pdf',width = 10,height = 10)
#par(mfrow=c(6,6))
# pheno = "Sepsis"
for (i in colnames(geneTraitSignificance)){
  moduleGenes = moduleColors=='cyan'
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, 'MEcyan']),
                     abs(geneTraitSignificance[moduleGenes,i]),
                     xlab = ("Module Membership in module"),
                     ylab = paste("Gene significance",pheno),
                     main = paste0('Trait ',substr(i,4,nchar(i)),"\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = 'cyan')
  abline(v=0.6,h=0.2,col="red")
}
dev.off()


pdf('Fig8.Gene.trait_MEgreen.pdf',width = 10,height = 10)
#par(mfrow=c(6,6))
pheno = "Sepsis"
for (i in colnames(geneTraitSignificance)){
  moduleGenes = moduleColors=='green'
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, 'MEgreen']),
                     abs(geneTraitSignificance[moduleGenes,i]),
                     xlab = ("Module Membership in module"),
                     ylab = paste("Gene significance",pheno),
                     main = paste0('Trait ',substr(i,4,nchar(i)),"\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = 'green')
  #abline(v=0.6,h=0.2,col="red")
}
dev.off()


### 输出GS_MM数据
probes <- colnames(dataExpr)
geneInfo0 <- data.frame(
  probes = probes,
  moduleColor = moduleColors
)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0, geneTraitSignificance[, Tra],
    GSPvalue[, Tra]
  )
  names(geneInfo0) <- c(
    oldNames, names(geneTraitSignificance)[Tra],
    names(GSPvalue)[Tra]
  )
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0, geneModuleMembership[, mod],
    MMPvalue[, mod]
  )
  names(geneInfo0) <- c(
    oldNames, names(geneModuleMembership)[mod],
    names(MMPvalue)[mod]
  )
}
geneOrder <- order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.csv", sep = ",", row.names = F)


########################2. 差异表达基因筛选#######

library(data.table)
library(limma)
library(ggplot2)
library(pheatmap)

exprSet <- read.csv("GSE65682_expr.csv",check.names = F,row.names = 1)
# View(exprSet)
# exprSet = apply(exprSet,2,function(x){log2(x+1)})

# exprSet = apply(ab,2,function(x){log2(x+1)})

exprSet <- as.matrix(exprSet)


group_list=c(rep('Healthy',42),rep('Sepsis',760))
group_list <- factor(group_list,levels = c("Sepsis","Healthy"))

# boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)

######差异表达
pvalue <-0.05
logFoldChange <- 0.5
dat <- exprSet
design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dat)
design
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(Sepsis-Healthy, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.csv(cbind(Symbol=rownames(allDiff),allDiff),file="limmaOut.csv")

diffSig = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.csv(cbind(Symbol=rownames(diffSig),diffSig), file="diffSig(0.5).csv")

diffUp = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC>logFoldChange)),]
write.csv(cbind(Symbol=rownames(diffUp),diffUp), file="up.csv")

diffDown = allDiff[(allDiff$adj.P.Val < pvalue & (allDiff$logFC<(-logFoldChange))),]
write.csv(cbind(Symbol=rownames(diffDown),diffDown), file="down.csv")

#火山图##

pvalue <-0.05
logFC <- 0.5
allDiff <- read.csv('limmaOut.csv',check.names = F,row.names = 1)
allDiff$Significant <- ifelse(allDiff$adj.P.Val<pvalue & abs(allDiff$logFC)>= logFC,
                              ifelse(allDiff$logFC> logFC,'up','down'),'stable')


mycol <- c("turquoise","#3D3D3D","coral1")
# mycol <- c("#20b2aa","#696969","#bc6e7f")

pdf(file="Sepsis vs. Healthy DEGs_volcano.pdf",width=5.5,height=5.5)
p <- ggplot(allDiff, aes(logFC, -log10(adj.P.Val), colour= Significant))+
  geom_point(size=2.2)+xlim(-4,4)+theme_bw()+
  scale_color_manual(values = mycol,name='Significant')+
  labs(title="TCGA-Sepsis-Healthy-DEGs",x="Log2(Fold Change)",y="-Log10 (adjust.P-Value)")+
  geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 0.5)+
  geom_vline(xintercept = c(-logFC, logFC), linetype=4,lwd = 0.5)+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 10),text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
print(p)
dev.off()



###############volcano########

table<-read.csv('2. 中性粒细胞相关基因功能分析/差异表达基因筛选/limmaOut.csv',check.names = F,row.names = 1)


library('ggplot2')
library('ggrepel')

table$group<-factor(ifelse(table$P.Value < 0.05 & abs(table$logFC) >=0.5, 
                           ifelse(table$logFC>=0.5 ,'Up','Down'),'None'),
                    levels=c('Up','Down','None'))
table(table$group)

list<-c('PLCG1',
        'CASP6',
        'GSDMB',
        'CASP4',
        'ELANE',
        'NLRP3')

pdf("2. 中性粒细胞相关基因功能分析//volcano.pdf",width = 8,height = 6)

ggplot(table,aes(x=logFC,y=-log10(P.Value),color=group))+
  geom_point(size=2)+
  geom_label_repel(data=table[rownames(table) %in% list,],min.segment.length = 0.001, 
                   segment.size = 0.5,box.padding = 0.15, #字或框和点的距离
                   max.overlaps=9999,aes(label=list),
                   size=3,segment.color='black',show.legend=F)+
  theme_classic(base_size = 18)+ 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_blank())+
  ylab('-log10(pvalue)')+
  xlab('log2(FoldChange)')+
  geom_vline(xintercept=c(-0,0),lty=3,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
  scale_color_manual(values=c('#E41A1C','#377EB8',"#b9b9b8"),name = "DEGs",
                     breaks=c("Up", "Down", "None"),
                     labels=c(paste0('Up',' (',nrow(table[table$group=='Up',]),')'),
                              paste0('Down',' (',nrow(table[table$group=='Down',]),')'),
                              paste0('None',' (',nrow(table[table$group=='None',]),')')))
dev.off()



###
########################2.差异关键模块基因韦恩图########

a <- read.csv("diffSig(0.5).csv",check.names = F)
b <- read.csv("关键模块基因.csv",check.names = F)


library(VennDiagram)
# summary(lasso_fea %in% top.features[1:which.min(errors), "FeatureName"]) 

pdf("DE-WGCNA.pdf", width = 6, height = 6)
grid.newpage()

venn.plot <- venn.diagram(list(Drug_targets = a$Symbol,b = b$probes),NULL, 
                          fill = c("turquoise","coral1"), 
                          alpha = c(0.8,0.8), cex = 4, cat.fontface=3, 
                          category.names = c("Sepsis-Healthy-DEGs","WGCNA"), 
                          main = "Overlap")
grid.draw(venn.plot)
dev.off()

# d <- read.csv("file:///E:/22年5月项目/51. XA0107/1. 差异NMRGs筛选及其功能/DE-NMRGs.csv",check.names = F)
# ad = b[b$Symbol %in% d$Symbol,]
# write.csv(ad,"DE-NMRGs.csv")


##############2.差异关键模块基因富集分析########

setwd("")

library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

###GO富集分析################
# .ID转换，ENTREZID是进行GO分析最好的ID类型（针对clusterProfiler）
deg.aging.related <- read.csv('差异关键模块基因.csv') %>% .$Symbol

target_gene_id<- bitr(deg.aging.related,fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
write.csv(target_gene_id,file='gene_entriz.csv')

result_go_aging_genes <- enrichGO(target_gene_id$ENTREZID,
                                  OrgDb = "org.Hs.eg.db",
                                  keyType = "ENTREZID",
                                  pvalueCutoff = 0.05,
                                  ont = "ALL",readable = T)

write.csv(result_go_aging_genes,"enrichGO.csv")


##画图必备函数
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

##展示term数
display_number = c(10, 10)
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$ENTREZID,
                   keyType = "ENTREZID",
                   pvalueCutoff = 0.2,
                   ont = "MF",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
# write.csv(ego_result_MF,"GO_MF.csv")
# ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
#                    gene = target_gene_id$ENTREZID,
#                    pvalueCutoff = 0.2,
#                    ont = "CC",
#                    readable=TRUE)
# ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]

ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id$ENTREZID,
                   pvalueCutoff = 0.2,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[2], ])
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]

go_enrich_df <- data.frame(
  ID = c(ego_result_BP$ID, ego_result_MF$ID),
  Description = c(ego_result_BP$Description,  ego_result_MF$Description),
  GeneNumber = c(ego_result_BP$Count,ego_result_MF$Count),
  type = factor(c(
    rep("biological process", display_number[2]),
    rep("molecular function", display_number[1])
  ), levels = c("molecular function", "biological process"))
)



## numbers as data on x axis

go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

labels <- (sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names
))

# names(labels) <- rev(1:nrow(go_enrich_df))
labels <- as.factor(rev(go_enrich_df$Description))



## colors for bar // green, blue, orange
CPCOLS <- c("#02B1e6", "#E81D22", "#F9BC15")

p <- ggplot(data = go_enrich_df, aes(x = number, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS) +
  theme_bw() +
  scale_x_discrete(labels = labels) +
  xlab("") +
  theme(axis.text = element_text(face = "bold", color = "gray50")) +
  labs(title = "The Most Enriched GO Terms")

p

pdf("go_barplot.pdf", width = 9, height = 6)
p
dev.off()

########按照p值上色
kk = result_go_aging_genes
pdf(file="GO_barplot_p.pdf",width = 9,height = 6)
# png(file="GO_barplot.png",width = 900,height = 600)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY",font.size = 10) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=80))+
  facet_grid(ONTOLOGY~., scale='free')

dev.off()


pdf(file="GO_bubble.pdf",width = 9,height = 7)
png(file="GO_bubble.png",width = 900,height = 700)
enrichplot::dotplot(kk,showCategory = 10,split="ONTOLOGY",font.size = 10)  #scale_x_discrete(labels=function(x) str_wrap(x, width=80))+
# facet_grid(ONTOLOGY~.,scale='free')
dev.off()

##按照p值
##用pvalue做填充
# kegg= data.frame(kk)
# View(kegg)


kegg <- read.csv("enrichGO.csv",row.names = 1)
kegg = kegg[kegg$ONTOLOGY == "BP",]
#对富集结果按照p.adjust进行从小到大排序，保证最显著的通路在前
kegg <- kegg[order(kegg$pvalue),]
#这里画图只展示top10的通路
kegg <- kegg[1:10,]
#提取每条通路里面差异表达的基因数
top10 <- data.frame(kegg$Description,kegg$Count ,kegg$pvalue)
colnames(top10) <- c("Description","count","pvalue")
#fill=padj fill颜色填充，使用连续值padj
p <- ggplot(data=top10,aes(x=Description,y=count,fill=pvalue))
#coord_flip()颠倒坐标轴
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))
#ylim(0,30) 更改横坐标的范围这里坐标轴颠倒了，虽然看起来是x轴，但其实是y轴
p3 <- p2 + ylim(0,32) + scale_fill_gradient(low="firebrick1",high="royalblue")
p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="",title="GO Term")
#输出为png格式的图片
png("GO_BP_barplot.png",width=500,height=480)
print(p4)
dev.off()
#输出为pdf的文件
pdf("GO_BP_barplot.pdf",width=9,height = 7)
print(p4)
dev.off()

###GO好看的气泡图
diffSig1 <- read.csv("diffSig(0.5).csv")
data <- read.csv('enrichGO.csv',check.names = F)
a <- data[data$ONTOLOGY == "BP",]
a=a[1:10,]
b = data[data$ONTOLOGY == "CC",]
b = b[1:10,]

c = data[data$ONTOLOGY == "MF",]
c = c[1:10,]

kkid_file = rbind(a,b,c)

kkid <- data.frame(Category = kkid_file$ONTOLOGY,
                   ID = kkid_file$ID,
                   Term = kkid_file$Description,
                   Genes = gsub('/',',',kkid_file$geneID),
                   adj_pval = kkid_file$p.adjust)

diffSig1 <- data.frame(ID = diffSig1$Symbol, logFC = diffSig1$logFC) 
cirGO <- circle_dat(kkid,diffSig1) 
#过滤基因重叠度等于或超过0.75的GO term
reduced_circ <- reduce_overlap(cirGO, overlap = 0.95)
GOBubble(reduced_circ, labels = 2.8)

pdf('enrichGO.Bubble.pdf',width = 12,height = 7)
GOBubble(cirGO, title = 'Gene Ontology Bubble plot', display = 'multiple',
bg.col = T, labels = 4)

# GOBubble(reduced_circ, labels = 3)
dev.off()

##KEGG
kegg <- read.csv("enrichKEGG.csv",row.names = 1)

#对富集结果按照p进行从小到大排序，保证最显著的通路在前
kegg <- kegg[order(kegg$pvalue),]
#这里画图只展示top10的通路
kegg <- kegg[1:5,]
#提取每条通路里面差异表达的基因数
top10 <- data.frame(kegg$Description,kegg$Count ,kegg$pvalue)
colnames(top10) <- c("Description","count","pvalue")
#fill=padj fill颜色填充，使用连续值padj
p <- ggplot(data=top10,aes(x=Description,y=count,fill=pvalue))
#coord_flip()颠倒坐标轴
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))
#ylim(0,30) 更改横坐标的范围这里坐标轴颠倒了，虽然看起来是x轴，但其实是y轴
p3 <- p2 + ylim(0,4) + scale_fill_gradient(low="firebrick1",high="royalblue")
p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="",title="KEGG Pathway")
#输出为png格式的图片
png("KEGG_barplot.png",width=500,height=480)
print(p4)
dev.off()
#输出为pdf的文件
pdf("KEGG_barplot.pdf",width=9,height = 7)
print(p4)
dev.off()

##KEGG富集分析#################

deg.aging.related <- read.csv('差异关键模块基因.csv') %>% .$Symbol


target_gene_id<- bitr(deg.aging.related,fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
result_KEGG <- enrichKEGG(target_gene_id$ENTREZID,organism = "hsa", pvalueCutoff =0.05)

write.csv(result_KEGG,"enrichKEGG.csv")

###将KEGG富集文件的entrize改为Symbol

# enrichOutput<- read.csv("enrichKEGG.csv",check.names = F)
# biotype <- read.csv("gene_entriz.csv",check.names = F)
# 
# enrichOutput$geneSymbol <- unlist(lapply(enrichOutput$geneID, function(v)
# 
#   paste(biotype$SYMBOL[match(strsplit(as.character(v), '/', fixed=TRUE)[[1]],
#                              biotype$ENTREZID)], collapse = '/')))
# write.csv(enrichOutput,"enrichKEGG.csv")

result_KEGG <- as.data.frame(result_KEGG)
result_KEGG = result_KEGG[1:5,]
result_KEGG$number <- factor(rev(1:nrow(result_KEGG)))

KEGG_enrich_df <- data.frame(
  ID = c(result_KEGG$ID),
  Description = c(result_KEGG$Description),
  GeneNumber = c(result_KEGG$Count),
  type = factor(c(
    rep("KEGG Pathway"))
  ), levels = c("KEGG Pathway"))


labels <- sapply(
  levels(KEGG_enrich_df$Description)[as.numeric(KEGG_enrich_df$Description)],
  shorten_names
)

# names(labels) <- rev(1:nrow(go_enrich_df))
labels <- as.factor(rev(KEGG_enrich_df$Description))
# DATA <- KEGG_enrich_df[1:20,]
## colors for bar // green, blue, orange
CPCOLS <- c("brown1")

p <- ggplot(data = KEGG_enrich_df, aes(x = ID, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.15) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS) +
  theme_bw() +
  scale_x_discrete(labels = labels) +
  xlab("") +
  theme(axis.text = element_text(face = "bold", color = "gray10")) +
  labs(title = "The Most Enriched KEGG Pathways")

p

pdf("KEGG_barplot.pdf", width = 10, height = 8)
p
dev.off()

########按照p值上色
k = result_KEGG
pdf(file="KEGG_barplot_p.pdf",width = 9,height = 6)
# png(file="KEGG_barplot.png",width = 900,height = 600)
barplot(k,showCategory = 15,font.size = 15,fill = pvalue)
dev.off()

##KEGG条形图
pdf(file="KEGG_bubble.pdf",width = 9,height = 6)
# png(file="KEGG_bubble.png",width = 900,height = 600)
enrichplot::dotplot(k,showCategory = 15,font.size = 15)
dev.off()


#KEGG八卦图

library(GOplot)
diffSig1 <- read.csv("diffSig(0.5).csv")
kkid_file <- read.csv('enrichKEGG.csv',check.names = F)
kkid <- data.frame(Category = 'All',
                   ID = kkid_file$ID,
                   Term = kkid_file$Description,
                   Genes = gsub('/',',',kkid_file$geneID),
                   adj_pval = kkid_file$p.adjust)
allentrez <- mapIds(org.Hs.eg.db,keys = diffSig1$X,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='first')
diffSig1 <- data.frame(ID = allentrez, logFC = diffSig1$logFC) 
cirkegg <- circle_dat(kkid,diffSig1) 

pdf('enrichKEGG.gossip.pdf',width = 10,height = 8)
GOCircle(cirkegg,rad1=2.5,rad2=3.5,label.size=3.5,nsub=15,zsc.col=c('blue', 'white','yellow'))  #nsub是富集出来的通路数目
dev.off()

#########################3. 中性粒细胞相关细胞焦亡基因筛选#######
setwd("")

a <- read.csv("PBMC33828074-细胞焦亡基因.csv",check.names = F)
b <- read.csv("差异关键模块基因.csv",check.names = F)


library(VennDiagram)
# summary(lasso_fea %in% top.features[1:which.min(errors), "FeatureName"]) 

pdf("DE-Pyroptosis.pdf", width = 6, height = 6)
grid.newpage()

venn.plot <- venn.diagram(list(Drug_targets = a$Genes,b = b$Symbol),NULL, 
                          fill = c("turquoise","coral1"), 
                          alpha = c(0.8,0.8), cex = 4, cat.fontface=3, 
                          category.names = c("DE-WGCNA","Pyroptosis-RGs"), 
                          main = "Overlap")
grid.draw(venn.plot)
dev.off()

d <- read.csv("候选基因.csv",check.names = F)
c <- read.csv("diffSig(0.5).csv",check.names = F)
cd = c[c$Symbol %in% d$Symbol,]
write.csv(cd,"DE-Pyroptosis.csv")


########################4.预后模型数据准备############
setwd("")

a <- read.csv("GSE65682_expr.csv",check.name = F)
b <- read.csv("候选基因.csv",check.name = F)

ab = a[a$ID %in% b$Symbol,]
ab = ab[!duplicated(ab$ID),]
ab = data.frame(t(ab))
write.csv(ab,"input1.csv")


a <- read.csv("input1.csv",check.name = F)
b <- read.csv("GSE65682-sample.csv",check.name = F)

ab = merge(a,b,by = "ID")
ab = ab[!duplicated(ab$ID),]
write.csv(ab,"input1.csv")


########################4. 通过循环筛选出具有生存意义的基因############
library(survival)
library(Hmisc)
rt<-read.csv("input-KM.csv",row.names = 1)

OSinput = rt
KM <- data.frame()

for(i in colnames(OSinput[,3:ncol(OSinput)])){
  
  # tryCatch({
    survival <- OSinput[,c('fustat','futime',i)]
    survival$exp <- ifelse(survival[,i]>quantile(survival[,i],0.5),
                           'High_Expression','Low_Expression')
    diff <- survdiff(Surv(futime,fustat) ~ exp, data = survival)
    pValue=format((1-pchisq(diff$chisq,df=1)),digits=5)
    KM <- rbind(KM,data.frame(Gene=i,pvalue=pValue))
  # },error=function(e){})
}
write.csv(KM,file = 'KM.csv')

########################4. 基因对样本的区分能力######


library(RColorBrewer)
library(pROC)

data <- read.csv('input1.csv',row.names = 1,check.names = F)


data$Type <- c(rep('1',760),rep('0',42))

genes <- c('CASP4','PLCG1','GSDMB','ELANE', 'CASP6','NLRP3')

for( gene in genes){
  
  # gene = "ELANE"
  hubexp = data.frame(exp = data[,gene],Type=data$Type)
  
  
  
  hubexp <- na.omit(hubexp)#删除缺失值
  
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp$exp,
                        Type = hubexp$Type) 
  #mycol <- brewer.pal(10,'Set3')
  #install.packages("RColorBrewer")
  
  mycol <- brewer.pal(8,'Set2')
  
  pdf(file = paste0(gene,'-区分样本ROC曲线.pdf'),width = 6,height = 6)
  
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="Sepsis vs. Healthy",
                # print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[3],#线的颜色
                lwd=3, #线的粗细
                legacy.axes=T)
  
  j=1
  auc.test <- paste0(gene,' AUC : ',format(as.numeric(x$auc),digits=3))
  
  legend(0.6,0.2, auc.test,lwd=2,bty="n",col=mycol[3],cex = 1.2)
  
  dev.off()
  
}
########################4. 基因的生存曲线#############

setwd("4. 中性粒细胞相关的细胞焦亡基因的生存")
library(survminer)
library(survival)
data <- read.csv("input-KM.csv",row.names = 1)



# 

genes <- c('CASP4','PLCG1','GSDMB','ELANE', 'CASP6','NLRP3')

for( gene in genes){
gene = "CASP4"
exp = data[,gene]
# group<- ifelse(exp > quantile(exp,0.5), 'high','low')

# ##确定最佳阈值
res.cut <- surv_cutpoint(data, time = "futime", event = "fustat",
                         variables = c(gene),minprop = 0.05,
                         progressbar = TRUE)
pdf(file = paste0(gene,"-最佳阈值.pdf"),width=6,height = 6)
p1 <-plot(res.cut,gene,palette = "npn")
print(p1)
dev.off()

xx = res.cut$cutpoint$cutpoint
group<- ifelse(data[,gene] > res.cut$cutpoint$cutpoint, 'high','low')

cox.data.plot = data.frame(futime= data$futime,fustat=data$fustat,exp,group)

# group <- ifelse(cox.data.plot$ELANE > quantile(cox.data.plot$CASP4,0.5), 'high','low')
# cox.data.plot$futime = cox.data.plot$futime/30
fit <- survfit(Surv(futime,fustat) ~ group, data = cox.data.plot)
survdf <- survdiff(Surv(futime,fustat) ~ group, data = cox.data.plot)
p.value <- 1 - pchisq(survdf$chisq, length(survdf$n) -1)



pdf(file=paste0(gene,"-28天生存曲线.pdf"),height = 9,width = 9)
# pdf(file="CGGA-生存曲线.pdf",height = 9,width = 9)

ggsurvplot(fit,
           # surv.median.line = "hv",
           conf.int=TRUE,
           pval=TRUE,
           risk.table=TRUE,
           #legend.labs=legend.labs,
           risk.table.title = "",
           xlab="Time (days)",
           legend.title=paste0(gene," Expression"),
           #linetype = lty,
           palette = c("coral1", "royalblue"),
           
           title="28 days survival in patients with sepsis",
           # title="CGGA dataset \r\n (Overall Survival)",
           risk.table.height=.15) %>% print()
dev.off()



}


########################5. 基因高低表达组的临床差异######


setwd("")

# a <- read.csv("file:///E:/22年5月项目/52. YQ345-4/input1.csv",check.names = F)
# b <- read.csv("input-KM.csv",check.names = F)
# ab = a[a$ID%in% b$ID,]
# write.csv(ab,"临床input.csv")

library(boot)
library(htmlTable)
library(Rcpp)
library(Gmisc)

data <- read.csv("临床input.csv",row.names = 1)



# 

genes <- c('CASP4','PLCG1','GSDMB','ELANE', 'CASP6','NLRP3')

for( gene in genes){
  gene = "NLRP3"
  exp = data[,gene]
  # group<- ifelse(exp > quantile(exp,0.5), 'high','low')
  
  # ##确定最佳阈值
  res.cut <- surv_cutpoint(data, time = "futime", event = "fustat",
                           variables = c(gene),minprop = 0.05,
                           progressbar = TRUE)

  
  xx = res.cut$cutpoint$cutpoint

annotation_col1=read.csv("临床input.csv",row.names=1,check.names = F)

annotation_col1$Risk = ifelse(exp > xx,"high",'low')
# View(annotation_col1)
annotation_col1$Risk <- factor(annotation_col1$Risk,
                               levels = c("high","low"),
                               labels = c("Risk-high","Risk-low"))

annotation_col1$gender <- factor(annotation_col1$gender)
annotation_col1$`pneumonia diagnoses` <- factor(annotation_col1$`pneumonia diagnoses`)
annotation_col1$endotype_class <- factor(annotation_col1$endotype_class)
annotation_col1$endotype_cohort <- factor(annotation_col1$endotype_cohort)

annotation_col1$diabetes_mellitus<- factor(annotation_col1$diabetes_mellitus)
annotation_col1$thrombocytopenia <- factor(annotation_col1$thrombocytopenia)
annotation_col1$icu_acquired_infection <- factor(annotation_col1$icu_acquired_infection)
annotation_col1$icu_acquired_infection_paired <- factor(annotation_col1$icu_acquired_infection_paired)



####自定义函数
getT1Stat <- function(varname, digits=1){
  getDescriptionStatsBy(annotation_col1[,varname],
                        annotation_col1$Risk,
                        add_total_col = TRUE,
                        show_all_values = TRUE,
                        hrzl_prop = FALSE,
                        useNA = "no",
                        statistics = TRUE,
                        statistics.sig_lim = 10^-3,
                        html=TRUE,
                        digits = digits
  )
}

table_data <- list()
table_data[["age(year)"]]<- getT1Stat("age")
table_data[["gender"]]<- getT1Stat("gender")
table_data[["pneumonia diagnoses"]]<- getT1Stat("pneumonia diagnoses")
table_data[["endotype_class"]]<- getT1Stat("endotype_class")
table_data[["endotype_cohort"]]<- getT1Stat("endotype_cohort")
table_data[["diabetes_mellitus"]]<- getT1Stat("diabetes_mellitus")
table_data[["thrombocytopenia"]] <- getT1Stat("thrombocytopenia")
table_data[["icu_acquired_infection"]]<- getT1Stat("icu_acquired_infection")

table_data[["icu_acquired_infection_paired"]]<- getT1Stat("icu_acquired_infection_paired")


rgroup <- c()
n.rgroup <- c()
output_data <- NULL

for(varlable in names(table_data)){
  output_data <- rbind(output_data,
                       table_data[[varlable]])
  
  rgroup <- c(rgroup,
              varlable)
  
  n.rgroup <- c(n.rgroup,
                nrow(table_data[[varlable]]))
}

cgroup <- c("","Risk","")
n.cgroup <- c(1,2,1)
colnames(output_data)<- gsub("[ ]*Risk","",colnames(output_data))


htmlTable(output_data,align = "rrrr",
          rgroup = rgroup,
          n.rgroup = n.rgroup,
          rgroupCSSseparator = "",
          n.cgroup = n.cgroup,
          rowlable = "",
          caption = paste0("Table 1 ",gene," and clinical data"),
          tfoot = "",
          ctable = TRUE
          
)
}

########################5. 基因表达与ICU住院率##########
library(ggplot2)
library(ggstatsplot)
library(ggpubr)

rt <- read.csv("临床input.csv",row.names = 1)
genes <- c('CASP4','PLCG1','GSDMB','ELANE', 'CASP6','NLRP3')

for( gene in genes){
  # gene = "NLRP3"
  exp = rt[,gene]
  # group<- ifelse(exp > quantile(exp,0.5), 'high','low')
  
  # ##确定最佳阈值
  res.cut <- surv_cutpoint(rt, time = "futime", event = "fustat",
                           variables = c(gene),minprop = 0.05,
                           progressbar = TRUE)
  
  
  xx = res.cut$cutpoint$cutpoint
  
  dataall=read.csv("临床input.csv",row.names=1,check.names = F)
  
  dataall$group = ifelse(exp > xx,"high",'low')
data = data.frame(icu_acquired_infection=dataall$icu_acquired_infection,group=dataall$group)

ggbarstats(data,icu_acquired_infection,group, palette = 'Set2',ylab = "",xlab =paste0(gene,"-Expression"))

ggsave(file=paste0(gene,"-高低风险组ICU住院率差异.pdf"),width = 6,height = 6)

}

#######################5. 基因表达与endotype class######

setwd("")


dbox1 <- read.csv("临床input.csv",check.names = F,row.names = 1)

table(dbox1$endotype_class)
compare_means(ELANE ~ endotype_class, data = dbox1)

my_comparition <- list(c("Mars1","Mars2"),c("Mars1","Mars3"),c("Mars1","Mars4")
                       ,c("Mars2","Mars3"),c("Mars2","Mars4"),c("Mars3","Mars4"))
mycol <- c("coral1","turquoise","cornflowerblue","lightpink")

pdf('endotype_class与基因表达.pdf',width = 9,height = 6)

# p <- ggboxplot(dbox1,x='immune',y='Score',fill ='Type',
#                ylab = 'Expression',
#                xlab ='',palette = mycol,
#                # add ='jitter',
#                size =0.4)+
#   rotate_x_text(45)+
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title = element_text(size=9),
#         legend.text = element_text(size=9),
#         legend.title = element_text(size=9),
#         axis.line = element_line(size=0.3))+
#   stat_compare_means(size = 3,aes(group=Type),
#                      label = 'p.signif',method = 'wilcox.test')
# 
# p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
#   font('x.text',face = 'bold')+font('y.text',face = 'bold')+
#   font('legend.title',face = 'bold')+font('legend.text',face = 'bold')

##小提琴图
# p <- dbox1 %>%
#   ggviolin(x= "group",y = c(colnames(dbox1)[1:9]),
#            fill = "group",combine = T,
#            palette = mycol,
#            ylab = "Nomolized Expression",
#            add = "boxplot",add.params = list(fill="white"))
# 
# p + stat_compare_means(method = "wilcox.test",
#                        label = "p.signif",comparisons = my_comparition)

##箱线图
# View(dbox1)
p <- dbox1 %>%
  ggboxplot(x= "endotype_class",y = c(colnames(dbox1)[3:8]),
            fill = "endotype_class",combine = T,
            palette = mycol,
            ylab = "Nomolized Expression",xlab = "endotype_class")

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition,size=4)

dev.off()




######################################6. 独立预后分析#######
setwd("")
library(survival)
library(survminer)
library(survivalROC)



##单因素
rt <- read.csv("临床input.csv",row.names=1,check.names = F)

# rt[,3:110] <- apply(rt[,3:110],2,function(x){log2(x+1)})###谨记！！一定要log

rt$futime=as.numeric(
  rt$futime)
rt$fustat=as.numeric(rt$fustat)


class(rt)
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.csv(outTab,file="uniCox_prognostic.csv")


##绘图


rt <- read.csv("uniCox_prognostic.csv",row.names = 1,check.names = F)
# rt <- read.csv("file:///E:/22年4月项目/YQ312-2/2. 风险模型构建/2. 预后模型构建及验证/uniCox.csv",row.names = 1)
# rt = rt[rt$pvalue < 0.05,]


gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
# pdf(file="单因素Cox回归分析森林图.pdf", width = 7,height = 4)
pdf(file="单因素独立预后森林图.pdf", width = 8,height = 6)
# 
# png(file="单因素独立预后森林图.png", width = 500,height = 300)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "coral1","seagreen1")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.0)
axis(1)
dev.off()


##多因素独立预后

rt1 <- read.csv('多因素input.csv',row.names = 1,check.names = F)


rt1$futime=as.numeric(
  rt1$futime)
rt1$fustat=as.numeric(rt1$fustat)


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

outTab=cbind(id=row.names(outTab),outTab)

write.csv(outTab,"mulCox_prognostic.csv")

rt <- read.csv("mulCox_prognostic.csv",row.names = 1,check.names = F)

gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="多因素独立预后森林图.pdf", width = 6,height = 4)

# png(file="多因素独立预后.png", width = 500,height = 300)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "coral1","seagreen1")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.0)
axis(1)
dev.off()

#####################################6. 列线图###########


library(rms)

rt <-read.csv("临床input.csv",row.names = 1,check.names = F)


pbc <-rt
dd <- datadist(pbc)
options(datadist="dd")
options(na.action="na.delete")
summary(pbc$futime)
coxpbc <- cph(formula = Surv(futime,fustat) ~  endotype_class+ELANE+age,data=pbc,x=T,y=T,surv = T,na.action=na.delete)


surv <- Survival(coxpbc) 
surv1 <- function(x) surv(7,x)
surv2 <- function(x) surv(14,x)
surv3 <- function(x) surv(28,x)

x <- nomogram(coxpbc,fun = list(surv1,surv2,surv3),lp=T,
              funlabel = c('7-days survival Probability','14-days survival Probability','28-days survival Probability'),
              maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))



pdf("预测脓毒症生存的列线图.pdf",width = 8,height = 8)

# png("TCGA_nomogram.png",width = 800,height = 600)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

#计算C指数
# set.seed(1)
v <- validate(coxpbc, dxy=TRUE, B=1000)
Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.corrected']
orig_Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.orig']
bias_corrected_c_index  <- abs(Dxy)/2+0.5  # 计算校正c-index
orig_c_index <- abs(orig_Dxy)/2+0.5  # 计算未校正c-index
bias_corrected_c_index
orig_c_index


#####################################6. 矫正曲线###########

pbcox=pbc


f1<- cph(formula = Surv(futime,fustat) ~   endotype_class+ELANE+age,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 7) 

#参数m=50表示每组50个样本进行重复计算
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=7,m=100,B=100) 

f2 <- cph(formula = Surv(futime,fustat) ~   endotype_class+ELANE+age,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 14) 
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=14,m=100,B=100)

f3 <- cph(formula = Surv(futime,fustat) ~   endotype_class+ELANE+age,data=pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 28) 
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=28,m=100,B=100)

pdf("7天生存校正曲线.pdf",width = 5,height = 5)
plot(cal1,lwd = 2,lty = 1,errbar.col = c("turquoise"),
     #bty = "l", #只画左边和下边框
     xlim = c(0.6,1),ylim= c(0.4,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("turquoise"),
     cex.lab=1.2,cex.axis=1, cex.main=1, cex.sub=0.5)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("turquoise"), pch = 10)
mtext("")
abline(0,1, lwd = 2, lty = 3, col = c("#224444"))
legend("bottomright", #图例的位置
       legend = c("7-days survival"), #图例文字
       col =c("turquoise"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "")#不显示图例边框
dev.off()

pdf("14天生存校正曲线.pdf",width = 5,height = 5)

plot(cal2,lwd = 3,lty = 1,errbar.col = c("cornflowerblue"),
     
     xlim = c(0.6,1),ylim= c(0.4,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("cornflowerblue"),
     cex.lab=1.2,cex.axis=1, cex.main=1, cex.sub=0.5)

lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("cornflowerblue"), pch = 16)
abline(0,1, lwd = 2, lty = 3, col = c("#224444"))
legend("bottomright", #图例的位置
       legend = c("14-days survival"), #图例文字
       col =c("cornflowerblue"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "")#不显示图例边框
dev.off()

pdf("28天生存校正曲线.pdf",width = 5,height = 5)
plot(cal3,lwd = 3,lty = 1,errbar.col = c("coral1"),
     
     xlim = c(0.4,1),ylim= c(0.4,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("coral1"),
     cex.lab=1.2,cex.axis=1, cex.main=1, cex.sub=0.5)


lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("coral1"), pch = 16)
abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("bottomright", #图例的位置
       legend = c("28-days survival"), #图例文字
       col =c("coral1"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "")#不显示图例边框
dev.off()


########   计算斜率
# library(stringr)
# caldat <- data.frame(summary(cal1))
# cal1rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -1) ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -1))[["coefficients"]][["(Intercept)"]]
# # summary(lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -1) ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -1)))
# caldat <- data.frame(summary(cal3))
# cal3rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -3)[1:6])[["coefficients"]][["(Intercept)"]]
# caldat <- data.frame(summary(cal2))
# cal2rate <- lm( str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -10, -3)[1:6])[["coefficients"]][["(Intercept)"]]



#######################7. ELANE高低表达组GSEA富集分析#######

setwd("")

matrix<-read.csv("GSE65682_expr.csv",check.names = F,row.names = 1)
# matrix<-read.csv("xCell.csv",check.names = F,row.names = 1)
# View(matrix)
clini<-read.csv("临床input.csv",check.names = F)##按照high-low排序

id <- clini$ID
ES.exp <- matrix[,id]
write.csv(cbind(Symbol=rownames(ES.exp),ES.exp),'Estimate.csv')


##差异表达分析
library(limma)
library(ggplot2)
library(pheatmap)
# exprSet <- read.csv("RPRD1B.csv",row.names = 1)

exprSet = ES.exp
# exprSet = apply(exprSet,2,function(x){log2(x+1)})
group_list=c(rep('high',376),rep('low',103))
group_list <- factor(group_list,levels = c("high","low"))

pvalue <-0.05
logFoldChange <-1
dat <- exprSet
design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dat)
design
# 构建差异比较矩阵
contrast.matrix <- makeContrasts(high-low, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.csv(cbind(Symbol=rownames(allDiff),allDiff),file="Risk_limmaOut.csv")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(RColorBrewer)

DEGs <- read.csv("Risk_limmaOut.csv",row.names = 1)

DEGs <- DEGs[order(DEGs$logFC,decreasing = T),]
entrezIDs <- mget(rownames(DEGs), org.Hs.egSYMBOL2EG, ifnotfound=NA)    
entrezIDs <- as.character(entrezIDs)
out=cbind(DEGs,entrezID=entrezIDs)
out = out[,c(1,7)]
out <- cbind(symbol = rownames(out),out)
write.table(out,file="GSEA_id.xls",sep="\t",quote=F,row.names = F)  

# GSEA_data <- read.table("GSEA_id.xls",sep="\t",header = T) 
# GSEA_data <- na.omit(GSEA_data)
# gene <- GSEA_data$logFC
# names(gene) <- GSEA_data$symbol

#gmt 文件--GO


GSEA_data <- read.table("GSEA_id.xls",sep="\t",header = T) 
GSEA_data <- na.omit(GSEA_data)
GSEA_gene <- GSEA_data$logFC
names(GSEA_gene) <- GSEA_data$entrezID

ggsea <-  gseGO(
  GSEA_gene,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  eps = 0,
  verbose = TRUE,
  seed = FALSE
  
)

ggsea.symbol <- setReadable(ggsea,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.csv(ggsea,'高低风险组GSEA.GO.csv')

pdf('高低风险组.GO.GSEA.pdf',width = 10,height = 10)
# png('高低风险组.GO.GSEA.png',width = 800,height = 750)
gseaplot2(ggsea, c('GO:0010880',
                   'GO:0010881',
                   'GO:1903514',
                   'GO:0014808',
                   'GO:0007158',
                   'GO:0086010',
                   'GO:0042340',
                   'GO:1905606',
                   'GO:0099174',
                   'GO:0090042'
                   
),
# rel_heights = c(1.5, 0.5, 0.5),
color= brewer.pal(10,'Paired'),pvalue_table = FALSE,rel_heights = c(1.5, 0.5, 0.5))
ggsave("高风险组TOP10-GO.pdf",width = 11,height = 10)
dev.off()

gseaplot2(ggsea,c('GO:0006271',
                  'GO:0006270',
                  'GO:0043038',
                  'GO:0000460',
                  'GO:0031055',
                  'GO:0000959',
                  'GO:0022616',
                  'GO:0030261',
                  'GO:0051383',
                  'GO:0034508'
),color= brewer.pal(10,'Paired'),pvalue_table = TRUE,rel_heights = c(1.5, 0.5, 0.5))

ggsave("低风险组TOP10-GO.pdf",width = 12,height = 12)
#KEGG--gmt 文件

# c2 <- clusterProfiler::read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
# gsea_KEGG <- GSEA(gene, TERM2GENE=c2, verbose=FALSE, pvalueCutoff = 0.05); head(gsea)
# write.csv(ggsea_KEGG,'高低风险组GSEA.KEGG.csv')

kegggse <-  gseKEGG(
  GSEA_gene,
  organism = "hsa",
  keyType = "kegg",
  
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = TRUE,
  seed = FALSE
)
kegggse.symbol <- setReadable(kegggse,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
write.csv(kegggse,'ELANE-GSEA.KEGG.csv')


gseaplot2(kegggse, c(1:10),pvalue_table = TRUE,
rel_heights = c(1.5, 0.5, 0.5),
color= brewer.pal(10,'Paired'))

ggsave("ELANE-GSEA-KEGG.pdf",width = 10,height = 11)
dev.off()

# gseaplot2(kegggse, c('hsa00020',
#                      'hsa03430',
#                      'hsa03420',
#                      'hsa00630',
#                      'hsa03050',
#                      'hsa03460',
#                      'hsa03410',
#                      'hsa03440',
#                      'hsa03030',
#                      'hsa00970'
#                      
#                      
# ),pvalue_table = TRUE,
# rel_heights = c(1.5, 0.5, 0.5),
# color= brewer.pal(10,'Paired'))
# 
# ggsave("低风险组TOP10-KEGG.pdf",width = 13,height = 14)

#######################7. ELANE-DEGs-KEGG富集分析#################

deg.aging.related <- read.csv('ELANE-DEGs.csv') %>% .$Symbol


target_gene_id<- bitr(deg.aging.related,fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db",drop = T)
write.csv(target_gene_id,"gene_entriz.csv")
result_KEGG <- enrichKEGG(target_gene_id$ENTREZID,organism = "hsa", pvalueCutoff =1)

write.csv(result_KEGG,"ELANE-enrichKEGG.csv")

###将KEGG富集文件的entrize改为Symbol

# enrichOutput<- read.csv("ELANE-enrichKEGG.csv",check.names = F)
# biotype <- read.csv("gene_entriz.csv",check.names = F)
# 
# enrichOutput$geneSymbol <- unlist(lapply(enrichOutput$geneID, function(v)
# 
#   paste(biotype$SYMBOL[match(strsplit(as.character(v), '/', fixed=TRUE)[[1]],
#                              biotype$ENTREZID)], collapse = '/')))
# write.csv(enrichOutput,"ELANE-enrichKEGG.csv")


##提取重叠通路对应基因

k = result_KEGG
keggenrich_out = k
keggenrich_sig <- dplyr::filter(keggenrich_out, pvalue < 0.05)
kegg_node <- keggenrich_sig[c(1,2,5,6), c(1,8)]
network_node <- tidyr::separate_rows(kegg_node, geneID,sep="/")
write.csv(network_node, "kegg_network_node.csv",quote = F,row.names = F)

# 属性文件
gene_type <- data.frame(terms = gene, type = rep("gene", length(gene)))
kegg_type <- data.frame(terms = keggenrich_sig$ID, type = rep("KEGG", length(keggenrich_sig$ID)))
type <- rbind(gene_type, kegg_type)
write.csv(type, "kegg_network_type.csv")

a <- read.csv("kegg_network_node.csv")
b <- read.csv("gene_entriz.csv")
ab = merge(a,b,by="geneID")

write.csv(ab,"kegg_network_node.csv")


##基因数目条形图
result_KEGG <- as.data.frame(result_KEGG)
result_KEGG = result_KEGG[1:15,]
result_KEGG$number <- factor(rev(1:nrow(result_KEGG)))

KEGG_enrich_df <- data.frame(
  ID = c(result_KEGG$ID),
  Description = c(result_KEGG$Description),
  GeneNumber = c(result_KEGG$Count),
  type = factor(c(
    rep("KEGG Pathway"))
  ), levels = c("KEGG Pathway"))


labels <- sapply(
  levels(KEGG_enrich_df$Description)[as.numeric(KEGG_enrich_df$Description)],
  shorten_names
)

# names(labels) <- rev(1:nrow(go_enrich_df))
labels <- as.factor(rev(KEGG_enrich_df$Description))
# DATA <- KEGG_enrich_df[1:20,]
## colors for bar // green, blue, orange
CPCOLS <- c("brown1")

p <- ggplot(data = KEGG_enrich_df, aes(x = ID, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.15) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS) +
  theme_bw() +
  scale_x_discrete(labels = labels) +
  xlab("") +
  theme(axis.text = element_text(face = "bold", color = "gray10")) +
  labs(title = "The Most Enriched KEGG Pathways")

p

pdf("KEGG_barplot.pdf", width = 10, height = 8)
p
dev.off()

########按照p值上色
k = result_KEGG
pdf(file="ELANE-KEGG_barplot_p.pdf",width = 9,height = 6)
# png(file="KEGG_barplot.png",width = 900,height = 600)
barplot(k,showCategory = 15,font.size = 15,fill = pvalue)
dev.off()

##KEGG条形图
pdf(file="KEGG_bubble.pdf",width = 9,height = 6)
# png(file="KEGG_bubble.png",width = 900,height = 600)
enrichplot::dotplot(k,showCategory = 15,font.size = 15)
dev.off()


#KEGG八卦图

library(GOplot)
diffSig1 <- read.csv("diffSig(0.5).csv")
kkid_file <- read.csv('enrichKEGG.csv',check.names = F)
kkid <- data.frame(Category = 'All',
                   ID = kkid_file$ID,
                   Term = kkid_file$Description,
                   Genes = gsub('/',',',kkid_file$geneID),
                   adj_pval = kkid_file$p.adjust)
allentrez <- mapIds(org.Hs.eg.db,keys = diffSig1$X,column = 'ENTREZID',keytype = 'SYMBOL',multiVals='first')
diffSig1 <- data.frame(ID = allentrez, logFC = diffSig1$logFC) 
cirkegg <- circle_dat(kkid,diffSig1) 

pdf('enrichKEGG.gossip.pdf',width = 10,height = 8)
GOCircle(cirkegg,rad1=2.5,rad2=3.5,label.size=3.5,nsub=15,zsc.col=c('blue', 'white','yellow'))  #nsub是富集出来的通路数目
dev.off()


#######################8. ELANE高低表达组细胞差异-mcp_counter#######

setwd("")


dbox1 <- read.csv("mcp_counter.csv",check.names = F,row.names = 1)

group <- c(rep("High",376),rep("Low",103))
length(group) = dim(dbox1)[[1]]
dbox1$group <- group
table(dbox1$group)
my_comparition <- list(c("High","Low"))
mycol <- c("coral1","turquoise")

pdf('mcp_counter细胞在ELANE高低表达组差异.pdf',width = 9,height = 8)

# p <- ggboxplot(dbox1,x='immune',y='Score',fill ='Type',
#                ylab = 'Expression',
#                xlab ='',palette = mycol,
#                # add ='jitter',
#                size =0.4)+
#   rotate_x_text(45)+
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title = element_text(size=9),
#         legend.text = element_text(size=9),
#         legend.title = element_text(size=9),
#         axis.line = element_line(size=0.3))+
#   stat_compare_means(size = 3,aes(group=Type),
#                      label = 'p.signif',method = 'wilcox.test')
# 
# p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
#   font('x.text',face = 'bold')+font('y.text',face = 'bold')+
#   font('legend.title',face = 'bold')+font('legend.text',face = 'bold')

##小提琴图
# p <- dbox1 %>%
#   ggviolin(x= "group",y = c(colnames(dbox1)[1:9]),
#            fill = "group",combine = T,
#            palette = mycol,
#            ylab = "Nomolized Expression",
#            add = "boxplot",add.params = list(fill="white"))
# 
# p + stat_compare_means(method = "wilcox.test",
#                        label = "p.signif",comparisons = my_comparition)

##箱线图
p <- dbox1 %>%
  ggboxplot(x= "group",y = c(colnames(dbox1)[1:11]),
           fill = "group",combine = T,
           palette = mycol,
           ylab = "Proportion",xlab = "Expression of ELANE")

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition,size=4)

dev.off()




#######################8. ELANE高低表达组细胞差异-quanTIseq#######

setwd("")


dbox1 <- read.csv("quanTIseq.csv",check.names = F,row.names = 1)

group <- c(rep("High",376),rep("Low",103))
length(group) = dim(dbox1)[[1]]
dbox1$group <- group
table(dbox1$group)
my_comparition <- list(c("High","Low"))
mycol <- c("coral1","turquoise")

pdf('quanTIseq细胞在ELANE高低表达组差异.pdf',width = 9,height = 8)

# p <- ggboxplot(dbox1,x='immune',y='Score',fill ='Type',
#                ylab = 'Expression',
#                xlab ='',palette = mycol,
#                # add ='jitter',
#                size =0.4)+
#   rotate_x_text(45)+
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title = element_text(size=9),
#         legend.text = element_text(size=9),
#         legend.title = element_text(size=9),
#         axis.line = element_line(size=0.3))+
#   stat_compare_means(size = 3,aes(group=Type),
#                      label = 'p.signif',method = 'wilcox.test')
# 
# p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
#   font('x.text',face = 'bold')+font('y.text',face = 'bold')+
#   font('legend.title',face = 'bold')+font('legend.text',face = 'bold')

##小提琴图
# p <- dbox1 %>%
#   ggviolin(x= "group",y = c(colnames(dbox1)[1:11]),
#            fill = "group",combine = T,
#            palette = mycol,
#            ylab = "Nomolized Expression",
#            add = "boxplot",add.params = list(fill="white"))
# 
# p + stat_compare_means(method = "wilcox.test",
#                        label = "p.signif",comparisons = my_comparition)

##箱线图
p <- dbox1 %>%
  ggboxplot(x= "group",y = c(colnames(dbox1)[1:11]),
            fill = "group",combine = T,
            palette = mycol,
            ylab = "Proportion",xlab = "Expression of ELANE")

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition,size=4)

dev.off()

#######################8. ELANE与quanTIseq细胞相关性分析#############

### 读入ssGSEA的结果

tcga_gsva<- read.csv("quanTIseq.csv",check.names = F,row.names = 1)

# View(tcga_gsva)
### 读入表达量数据
a <- read.csv("Estimate.csv",check.names = F,row.names = 1)
gene <- "ELANE"

ab = a[gene,]
tcga_expr = ab[,-1]
# View(tcga_expr)
# tcga_expr<- read.csv("m5C_expr_DE.csv",check.names = F,row.names = 1)##基因是一列
tcga_expr = as.matrix(tcga_expr)


### 读入感兴趣的基因
# genelist <- read.csv("m5C_expr_DE.csv",check.names = F)
genelist <-gene
##写函数
gene <- genelist
immuscore <- function(gene){
  y <- tcga_expr[gene,]
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
# immuscore("")
### 批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
# p = 0.01
# cor = 0.6
# data1 = data[(data$p.value < p & (data$cor > cor |  data$cor < -(cor))),]
write.csv(data, "ELANE与quanTIseq免疫细胞相关性.csv", quote = F, row.names = F)

data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

##画图
# library(ggplot2)
# 
# pdf('lncRNA_差异免疫细胞相关性热图.pdf',width = 7.5,height = 5.5 )
# ggplot(data, aes(immune_cells, gene)) + 
#   geom_tile(aes(fill = cor))+
#   scale_fill_gradient2(low = "#4b0082",mid = "white",high = "#ff4500")+
#   geom_text(aes(label=pstar),col ="#2f4f4f",size = 5)+
#   theme_minimal()+# 不要背景
#   theme(axis.title.x=element_blank(),#不要title
#         axis.ticks.x=element_blank(),#不要x轴
#         axis.title.y=element_blank(),#不要y轴
#         axis.text.x = element_text(angle = 45, hjust = 1,colour = "#000000"),# 调整x轴文字
#         axis.text.y = element_text(size = 8,colour = "#000000"))+#调整y轴文字
#   #调整legen
#   labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
# dev.off()


#######################8. ELANE与mcp_counter细胞相关性分析#############

### 读入ssGSEA的结果

tcga_gsva<- read.csv("mcp_counter.csv",check.names = F,row.names = 1)

# View(tcga_gsva)
### 读入表达量数据
a <- read.csv("Estimate.csv",check.names = F,row.names = 1)
gene <- "ELANE"

ab = a[gene,]
tcga_expr = ab[,-1]
# View(tcga_expr)
# tcga_expr<- read.csv("m5C_expr_DE.csv",check.names = F,row.names = 1)##基因是一列
tcga_expr = as.matrix(tcga_expr)


### 读入感兴趣的基因
# genelist <- read.csv("m5C_expr_DE.csv",check.names = F)
genelist <-gene
##写函数
gene <- genelist
immuscore <- function(gene){
  y <- tcga_expr[gene,]
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
# immuscore("")
### 批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
# p = 0.01
# cor = 0.6
# data1 = data[(data$p.value < p & (data$cor > cor |  data$cor < -(cor))),]
write.csv(data, "ELANE与mcp_counter免疫细胞相关性.csv", quote = F, row.names = F)

data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

##画图
# library(ggplot2)
# 
# pdf('lncRNA_差异免疫细胞相关性热图.pdf',width = 7.5,height = 5.5 )
# ggplot(data, aes(immune_cells, gene)) + 
#   geom_tile(aes(fill = cor))+
#   scale_fill_gradient2(low = "#4b0082",mid = "white",high = "#ff4500")+
#   geom_text(aes(label=pstar),col ="#2f4f4f",size = 5)+
#   theme_minimal()+# 不要背景
#   theme(axis.title.x=element_blank(),#不要title
#         axis.ticks.x=element_blank(),#不要x轴
#         axis.title.y=element_blank(),#不要y轴
#         axis.text.x = element_text(angle = 45, hjust = 1,colour = "#000000"),# 调整x轴文字
#         axis.text.y = element_text(size = 8,colour = "#000000"))+#调整y轴文字
#   #调整legen
#   labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
# dev.off()


#######################8. 相关性棒棒糖图#############

library(ggplot2)

Corr <- read.csv("ELANE与mcp_counter免疫细胞相关性.csv",check.names = F)

# p <- ggplot(Corr,aes(cor,immune_cells,fill=p.value)) +
#   geom_segment(aes(xend=0,yend=immune_cells))+
#   geom_point(shape=21,aes(size = cor,colour=p.value))+
#   scale_fill_gradient2(mid = "LightSeaGreen",high = "yellow")+
#   scale_color_gradient2(mid = "LightSeaGreen",high = "yellow") +
#   scale_size_continuous(range = c(3,8))+
#   ggtitle('FPI')+
#   theme(
#     axis.title=element_text(size=13,face="plain",color="black"),
#     axis.text = element_text(size=13,face="plain",color="black"),
#     legend.title=element_text(size=12,face="plain",color="black"),
#     legend.background = element_blank(),
#     legend.position = "left"
#   )+theme_bw()
# # theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))
# 
# print(p)
# dev.off()
# 
# 
# 
# View(data)

data = Corr
data <- data %>% 
  mutate(Correlation = ifelse(cor>0, "1-positive", "negative")) #设置分组

## fig 2 ：简单优化后的图形
fig4 <- ggplot(data, aes(x=reorder(immune_cells,cor), y=cor,fill = Correlation)) +
  geom_segment( aes(x=reorder(immune_cells,cor), xend=reorder(immune_cells,cor), y=0, yend=cor,color=Correlation),
                size=0.5,linetype=1)+#使用reorder()排序变量
  
  geom_point(shape=21,aes(size = -log10(p.value),colour=Correlation)) + 
  #在散点上显示数值并保留两位小数
  geom_text(aes(label =sprintf("%.2f",cor)), color = "black", size = 2.4)+ 
  xlab("Immunocyte") +
  scale_y_continuous("Correlation",breaks  = c(-0.5,-0.25,-0.1,0,0.1,0.25,0.5)) +
  coord_flip() +
  theme_minimal()+
  theme(
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=13,face="plain",color="black"),
        legend.title=element_text(size=12,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position = "left"
      )+theme_bw()
   

fig4
ggsave(file="ELANE与mcp_counter免疫细胞相关棒棒糖图.pdf",width = 6,height = 6)


###ELANE与quanTIseq免疫细胞相关性
Corr <- read.csv("风险评分与细胞相关性.csv",check.names = F)

data = Corr
data <- data %>% 
  mutate(Correlation = ifelse(cor>0, "1-positive", "negative")) #设置分组

## fig 2 ：简单优化后的图形
fig4 <- ggplot(data, aes(x=reorder(immune_cells,cor), y=cor,fill = Correlation)) +
  geom_segment( aes(x=reorder(immune_cells,cor), xend=reorder(immune_cells,cor), y=0, yend=cor,color=Correlation),
                size=0.5,linetype=1)+#使用reorder()排序变量
  
  geom_point(shape=21,aes(size = -log10(p.value),colour=Correlation)) + 
  #在散点上显示数值并保留两位小数
  geom_text(aes(label =sprintf("%.2f",cor)), color = "black", size = 2.3)+ 
  xlab("Immunocyte") +
  scale_y_continuous("Correlation",breaks  = c(-0.5,-0.25,-0.1,0,0.1,0.25,0.5)) +
  coord_flip() +
  theme_minimal()+
  theme(
    axis.title=element_text(size=13,face="plain",color="black"),
    axis.text = element_text(size=13,face="plain",color="black"),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background = element_blank(),
    legend.position = "left"
  )+theme_bw()


fig4
ggsave(file="风险评分与细胞相关棒棒糖图.pdf",width = 6,height = 6)

#######################9. ELANE的表达验证################
#######################8. ELANE高低表达组细胞差异-quanTIseq#######

library(ggplot2)
library(ggpubr)
dbox1 <- read.csv("GSE95233-ELANE.csv",check.names = F,row.names = 1)

group <- c(rep("Healthy",22),rep("Sepsis",102))
length(group) = dim(dbox1)[[1]]
dbox1$group <- group
table(dbox1$group)
my_comparition <- list(c("Sepsis","Healthy"))
mycol <- c("turquoise","coral1")

pdf('GSE95233-ELANE表达.pdf',width = 5,height = 5)

# p <- ggboxplot(dbox1,x='immune',y='Score',fill ='Type',
#                ylab = 'Expression',
#                xlab ='',palette = mycol,
#                # add ='jitter',
#                size =0.4)+
#   rotate_x_text(45)+
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title = element_text(size=9),
#         legend.text = element_text(size=9),
#         legend.title = element_text(size=9),
#         axis.line = element_line(size=0.3))+
#   stat_compare_means(size = 3,aes(group=Type),
#                      label = 'p.signif',method = 'wilcox.test')
# 
# p + font('xlab',face = 'bold')+font('ylab',face = 'bold')+
#   font('x.text',face = 'bold')+font('y.text',face = 'bold')+
#   font('legend.title',face = 'bold')+font('legend.text',face = 'bold')

##小提琴图
# p <- dbox1 %>%
#   ggviolin(x= "group",y = c(colnames(dbox1)[1:11]),
#            fill = "group",combine = T,
#            palette = mycol,
#            ylab = "Nomolized Expression",
#            add = "boxplot",add.params = list(fill="white"))
# 
# p + stat_compare_means(method = "wilcox.test",
#                        label = "p.signif",comparisons = my_comparition)

##箱线图
p <- dbox1 %>%
  ggboxplot(x= "group",y = c(colnames(dbox1)[1]),
            fill = "group",combine = T,
            palette = mycol,
            ylab = "ELANE Expression",xlab = "GSE95233")

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition,size=4)

dev.off()


###############################补充. 中性粒细胞的比较箱线图######

data <- read.csv("GSE65682_expr.csv",
                 check.names = F,row.names=1)


xx <- "ELANE"

dbox1 = data[xx,]
dbox1 = data.frame(t(dbox1))
group =c(rep('Healthy',42),rep('Sepsis',760))



length(group) = dim(dbox1)[[1]]
dbox1$group <- group
table(dbox1$group)
my_comparition <- list(c("Sepsis","Healthy"))

mycol <- c("turquoise","coral1")

p <- dbox1 %>%
  ggboxplot(x= "group",y = c(colnames(dbox1)[1]),
            fill = "group",combine = T,
            palette = mycol,
            ylab = "Normalized Expression")

p + stat_compare_means(method = "wilcox.test",
                       label = "p.signif",comparisons = my_comparition)
ggsave(file="GSE65682-ELANE的箱线图.pdf",width = 5,height = 5)



