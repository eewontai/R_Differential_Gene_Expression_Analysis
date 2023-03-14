#### install package 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("limma")
#BiocManager::install("GEOquery")
#install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz',repos=NULL)
#BiocManager::install("oligo")
#BiocManager::install("EnhancedVolcano")


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("EnhancedVolcano")

library(GEOquery)
library(limma)
library(affy)
library(oligo)

setwd("C:/Users/USER/Desktop/GSE22575_RAW")  
getwd()

# raw data download (.cel file)
getGEOSuppFiles("GSE22575")
untar("GSE22575/GSE22575_RAW.tar", exdir = "suppdata")
cels = list.files("suppdata/", pattern = "[gz]")
sapply(paste("suppdata", cels, sep = "/"), gunzip)
exp.cel_22575 <- read.celfiles(list.celfiles("suppdata/", full.names = T))
#exp.cel_22575 <- read.celfiles(list.celfiles("/home/user/athero_periodontitis/periodontitis/GSE22575", full.names = T))

# Feature and pheno data download  
GSE22575<- getGEO("GSE22575") 
GSE22575 <- GSE22575[[1]]

#data normalization 
exp.rma_22575 <- rma(exp.cel_22575)

#Quality assessment
#histogram
hist(exprs(GSE22575))
hist(exprs(exp.rma_22575),breaks =10)

#boxplot
boxplot(exprs(GSE22575))
boxplot(exprs(exp.rma_22575))

#MVA plot 
mva.pairs(exprs(GSE22575)[,1:4],log.it = TRUE)
mva.pairs(exprs(exp.rma_22575)[,1:4],log.it = TRUE)


#dataframe 

eset22575<-exprs(exp.rma_22575)

eset22575 <- data.frame(eset22575)

head(eset22575)

pset22575 <- pData(GSE22575)
fset22575 <- fData(GSE22575)

head(fset22575,3) # Probe ID, Gene symbol 
head(pset22575) # tissue:ch1 

fset22575_id <- fset22575[, c("ID", "Gene Symbol")]

eset22575$gene <- fset22575_id[,2]

eset22575 <- aggregate(. ~ gene, data = eset22575, mean)

eset22575[1:10,]
nrow(eset22575)

rownames(eset22575) <- eset22575$gene

eset22575 <- eset22575[,-1]
#???????????????????????????????????????????????????????????????????????????????????????????????????????????
grp_22575 <- pset22575$'source_name_ch1'
head(pset22575)
grp_22575
#####limma

design_22575 <- model.matrix(~0+ grp_22575)

colnames(design_22575)
head(design_22575)
#????????????????????????????????????????????????????????????????????????????????????????????????????????????
colnames(design_22575) <- c("ctrl","test")

fit_22575 <- lmFit(eset22575,design_22575)

design_22575
cont_22575 <- makeContrasts(test-ctrl,levels=design_22575)

fit.cont_22575 <- contrasts.fit(fit_22575,cont_22575)

fit.cont_22575 <- eBayes(fit.cont_22575)

res_22575 <- topTable(fit.cont_22575,number=Inf,lfc=1,adjust="BH")

head(res_22575)
# p < 0.05 & log2FC > abs(2) = DEGs

up_res22575 <-res_22575[res_22575$P.Value < 0.05 & res_22575$logFC > 1.5,]
down_res22575 <-res_22575[res_22575$P.Value < 0.05 & res_22575$logFC < -1.5,]

nrow(up_res22575)
nrow(down_res22575)

write.csv(up_res22575, file='C:/Users/USER/Desktop/GSE22575_RAW/up_res22575.csv')
write.csv(down_res22575, file='C:/Users/USER/Desktop/GSE22575_RAW/down_res22575.csv')

#devtools::install_github('kevinblighe/EnhancedVolcano')

#volcanoplot
#??????????????????????????????????????????????????????????????????????????????????????????????/
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res_22575,
                lab = rownames(res_22575),
                x = 'logFC',
                y = 'P.Value',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 2.0,
                labSize = 5.0,
                colAlpha = 1,
#                legend=c('NS','Log (base 2) fold-change','P value',
#                       'P value & Log (base 2) fold-change'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)

# downstream analysis_GO(Gene Ontology) analysis 


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
#install.packages("org.Mm.eg.db")
#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
#nrow(up_res22575)
#Input genes; convert to ENTREZID 
eg = bitr(rownames(up_res22575), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
eg1 = bitr(rownames(down_res22575), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
#head(eg)
write.csv(eg,file='C:/Users/USER/Desktop/GSE22575_RAW/up_id.csv')
write.csv(eg1,file='C:/Users/USER/Desktop/GSE22575_RAW/down_id.csv')

#GO-Enrichment analysis 
go <- enrichGO(eg$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all")
#keyType ='SYMBOL',pAdjustMethod = "BH",
#pvalueCutoff  = 0.01, qvalueCutoff  = 0.05
head(go)
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all")

#ego2 <- enrichGO(gene         = rownames(up_res22575),
#                 OrgDb         = org.Hs.eg.db,
#                 keyType       = 'SYMBOL',
#                 ont           = "all",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05)

dotplot(go, split="ONTOLOGY", showCategory=10) + facet_grid(ONTOLOGY~., scale="free")
dotplot(go1, split="ONTOLOGY", showCategory=10) + facet_grid(ONTOLOGY~., scale="free")

kk <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
dotplot(kk)










