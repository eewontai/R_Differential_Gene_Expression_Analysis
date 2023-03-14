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

setwd("C:/Users/USER/Desktop/GSE48397_RAW")  
getwd()

# raw data download (.cel file)
getGEOSuppFiles("GSE48397")
untar("GSE48397/GSE48397_RAW.tar", exdir = "suppdata")
cels = list.files("suppdata/", pattern = "[gz]")
sapply(paste("suppdata", cels, sep = "/"), gunzip)
exp.cel_48397 <- read.celfiles(list.celfiles("suppdata/", full.names = T))
#exp.cel_48397 <- read.celfiles(list.celfiles("/home/user/athero_periodontitis/periodontitis/GSE48397", full.names = T))

# Feature and pheno data download  
GSE48397<- getGEO("GSE48397") 
GSE48397 <- GSE48397[[1]]

#data normalization 
exp.rma_48397 <- rma(exp.cel_48397)

#Quality assessment
#histogram
hist(exprs(GSE48397))
hist(exprs(exp.rma_48397),breaks =10)

#boxplot
boxplot(exprs(GSE48397))
boxplot(exprs(exp.rma_48397))

#MVA plot 
mva.pairs(exprs(GSE48397)[,1:4],log.it = TRUE)
mva.pairs(exprs(exp.rma_48397)[,1:4],log.it = TRUE)


#dataframe 

eset48397<-exprs(exp.rma_48397)

eset48397 <- data.frame(eset48397)

head(eset48397)

pset48397 <- pData(GSE48397)
fset48397 <- fData(GSE48397)

head(fset48397,3) # Probe ID, Gene symbol 
head(pset48397) # tissue:ch1 

fset48397_id <- fset48397[, c("ID", "Gene Symbol")]

eset48397$gene <- fset48397_id[,2]

eset48397 <- aggregate(. ~ gene, data = eset48397, mean)

eset48397[1:10,]
nrow(eset48397)

rownames(eset48397) <- eset48397$gene

eset48397 <- eset48397[,-1]
#???????????????????????????????????????????????????????????????????????????????????????????????????????????
grp_48397 <- pset48397$'characteristics_ch1.1'
head(pset48397)
grp_48397
#####limma

design_48397 <- model.matrix(~0+ grp_48397)

colnames(design_48397)
head(design_48397)
#????????????????????????????????????????????????????????????????????????????????????????????????????????????
colnames(design_48397) <- c("ctrl","test")

fit_48397 <- lmFit(eset48397,design_48397)

design_48397


cont_48397 <- makeContrasts(test-ctrl,levels=design_48397)

fit.cont_48397 <- contrasts.fit(fit_48397,cont_48397)

fit.cont_48397 <- eBayes(fit.cont_48397)

res_48397 <- topTable(fit.cont_48397,number=Inf,lfc=1,adjust="BH")

head(res_48397)
# p < 0.05 & log2FC > abs(2) = DEGs

up_res48397 <-res_48397[res_48397$P.Value < 0.05 & res_48397$logFC > 1.5,]
down_res48397 <-res_48397[res_48397$P.Value < 0.05 & res_48397$logFC < -1.5,]

nrow(up_res48397)
nrow(down_res48397)

write.csv(up_res48397, file='C:/Users/USER/Desktop/GSE48397_RAW/up_res48397.csv')
write.csv(down_res48397, file='C:/Users/USER/Desktop/GSE48397_RAW/down_res48397.csv')

#devtools::install_github('kevinblighe/EnhancedVolcano')

#volcanoplot
#??????????????????????????????????????????????????????????????????????????????????????????????/
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res_48397,
                lab = rownames(res_48397),
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

BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
install.packages("org.Mm.eg.db")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
#nrow(up_res48397)
#Input genes; convert to ENTREZID 
eg = bitr(rownames(up_res48397), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
eg1 = bitr(rownames(down_res48397), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
#head(eg)
write.csv(eg,file='C:/Users/USER/Desktop/GSE48397_RAW/up_id.csv')
write.csv(eg1,file='C:/Users/USER/Desktop/GSE48397_RAW/down_id.csv')

#GO-Enrichment analysis 
go <- enrichGO(eg$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all")
#keyType ='SYMBOL',pAdjustMethod = "BH",
#pvalueCutoff  = 0.01, qvalueCutoff  = 0.05
head(go)
go1 <- enrichGO(eg1$ENTREZID, OrgDb = "org.Mm.eg.db", ont="all")

#ego2 <- enrichGO(gene         = rownames(up_res48397),
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
                 pvalueCutoff = 0.1)
dotplot(kk)
