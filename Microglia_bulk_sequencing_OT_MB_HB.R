# Load in required packaged: 
library("DESeq2")
library(org.Dr.eg.db)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
print(sessionInfo())

# 1. Generate PCA, get normalized counts and heatmaps:
# Set up Sample Table:
directory<- ("~/Count_files")
samplename<- c("OT_1", "OT_2", "OT_3", "MB_1", "MB_2", "MB_3", "HB_1", "HB_2", "HB_3")
samples <- c("OT_1_count.txt", "OT_2_count.txt", "OT_3_count.txt", "MB_1_count.txt", "MB_2_count.txt", "MB_3_count.txt", "HB_1_count.txt", "HB_2_count.txt", "HB_3_count.txt")
region <- c("OT", "OT", "OT", "MB", "MB", "MB", "HB", "HB", "HB")
date <- c("Day1", "Day2", "Day3", "Day1", "Day2", "Day3", "Day1", "Day2", "Day3")
sampleTable<-data.frame(sampleName=samplename, fileName=samples, region=region, date=date)
sampleTable

#Build DESeq2 object from count files:
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design= ~date + region)
ddsHTSeq

# Generate normalized counts:
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
sizeFactors(ddsHTSeq) # Look at the normalization factor applied to each sample
normalized_counts <- counts(ddsHTSeq, normalized=TRUE)

# Get gene names and export table:
ids=rownames(normalized_counts)
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Dr.eg.db
selRes <- AnnotationDbi::select(db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey))
x=selRes[ match( ids, selRes[,1] ), 2 ]
normalized_countsTable=data.frame(Ensembl_ID=rownames(normalized_counts),Gene_ID=x,normalized_counts)
write.table(normalized_countsTable, file="~/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Make a heatmap for purity check depicting a select set of genes (normalized counts + 1) with an average per region for each gene = 
dat=read.csv("~/normalized_counts+1_AVG.csv")
rownames(dat)=as.character(dat[,2])
dat=dat[,c(3:dim(dat)[2])]

# Generate heatmap: 
col_map = colorRamp2(c(0, 8, 15), c("#313695", "#676B93", "#D10909"))
column_order = c("OT", "MB", "HB")
Heatmap(dat, name = "Log2 (normalized counts + 1)", column_order = column_order, row_order = rownames(dat), col = col_map, column_title = "Purity", border = TRUE, width = unit(6, "cm"), height = unit(10, "cm"))

# Remove Batch effect of day of sample prep: 
vsd <- vst(ddsHTSeq)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$date)
plotPCA(vsd, "region")
p <- plotPCA(vsd, intgroup=c("region"))
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
print(p)

# 2. Compare OT to HB:
directory<- ("~/Count_files")
samplename<- c("OT_1", "OT_2", "OT_3", "HB_1", "HB_2", "HB_3")
samples <- c("OT_1_count.txt", "OT_2_count.txt", "OT_3_count.txt", "HB_1_count.txt", "HB_2_count.txt", "HB_3_count.txt")
region <- c("OT", "OT", "OT", "HB", "HB", "HB")
date <- c("Day1", "Day2", "Day3", "Day1", "Day2", "Day3")
sampleTableOvH<-data.frame(sampleName=samplename, fileName=samples, region=region, date=date)
sampleTableOvH

#BUILD A DESEQ2 OBJECT:  
ddsHTSeqOTvsHB<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTableOvH, directory=directory, design= ~date + region)
ddsHTSeqOTvsHB

# Differential expression
DEOTvsHB <- DESeq(ddsHTSeqOTvsHB)
resultsNames(DEOTvsHB)

# Do individual comparisons with p-value of 0.01: 
OTvsHB <- results(DEOTvsHB, contrast=c("region","OT","HB"), alpha=0.01) # if p-value of 0.05 set alpha=0.05
sum(OTvsHB$padj < 0.01, na.rm=TRUE)

ids=rownames(OTvsHB)
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Dr.eg.db
selRes <- AnnotationDbi::select(db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey))
x=selRes[ match( ids, selRes[,1] ), 2 ]

#Add gene symbols to differences
DEgenesOTvsHB=data.frame(Ensembl_ID=rownames(OTvsHB),Gene_ID=x,OTvsHB)
write.table(DEgenesOTvsHB, "~/DEgenesOTvsHB001.txt", sep="\t")

# 3. Compare OT to MB:
directory<- ("/Count_files")
samplename<- c("OT_1", "OT_2", "OT_3", "MB_1", "MB_2", "MB_3")
samples <- c("OT_1_count.txt", "OT_2_count.txt", "OT_3_count.txt", "MB_1_count.txt", "MB_2_count.txt", "MB_3_count.txt")
region <- c("OT", "OT", "OT", "MB", "MB", "MB")
date <- c("Day1", "Day2", "Day3", "Day1", "Day2", "Day3")
sampleTableOvM<-data.frame(sampleName=samplename, fileName=samples, region=region, date=date)
sampleTableOvM

#BUILD A DESEQ2 OBJECT:  
ddsHTSeqOTvsMB<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTableOvM, directory=directory, design= ~date + region)
ddsHTSeqOTvsMB

# Differential expression
DEOTvsMB <- DESeq(ddsHTSeqOTvsMB)
resultsNames(DEOTvsMB)

# Do individual comparisons with p-value of 0.01: 
OTvsMB <- results(DEOTvsMB, contrast=c("region","OT","MB"), alpha=0.01) # if p-value of 0.05 set alpha=0.05
sum(OTvsMB$padj < 0.01, na.rm=TRUE)

ids=rownames(OTvsMB)
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Dr.eg.db
selRes <- AnnotationDbi::select(db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey))
x=selRes[ match( ids, selRes[,1] ), 2 ]

#Add gene symbols to differences
DEgenesOTvsMB=data.frame(Ensembl_ID=rownames(OTvsMB),Gene_ID=x,OTvsMB)
write.table(DEgenesOTvsMB, "~/DEgenesOTvsMB001.txt", sep="\t")

# 4. Compare HB to MB:
directory<- ("~/Count_files")
samplename<- c("MB_1", "MB_2", "MB_3", "HB_1", "HB_2", "HB_3")
samples <- c("MB_1_count.txt", "MB_2_count.txt", "MB_3_count.txt", "HB_1_count.txt", "HB_2_count.txt", "HB_3_count.txt")
region <- c("MB", "MB", "MB", "HB", "HB", "HB")
date <- c("Day1", "Day2", "Day3", "Day1", "Day2", "Day3")
sampleTableHvM<-data.frame(sampleName=samplename, fileName=samples, region=region, date=date)
sampleTableHvM

#BUILD A DESEQ2 OBJECT:  
ddsHTSeqHBvsMB<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTableHvM, directory=directory, design= ~date + region)
ddsHTSeqHBvsMB

# Differential expression
DEHBvsMB <- DESeq(ddsHTSeqHBvsMB)
resultsNames(DEHBvsMB)

# Do individual comparisons with p-value of 0.01: 
HBvsMB <- results(DEHBvsMB, contrast=c("region","HB","MB"), alpha=0.01) # if p-value of 0.05 set alpha=0.05
sum(HBvsMB$padj < 0.01, na.rm=TRUE)

ids=rownames(HBvsMB)
fromKey="ENSEMBL"
toKey="SYMBOL"
db=org.Dr.eg.db
selRes <- AnnotationDbi::select(db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey))
x=selRes[ match( ids, selRes[,1] ), 2 ]

#Add gene symbols to differences
DEgenesHBvsMB=data.frame(Ensembl_ID=rownames(HBvsMB),Gene_ID=x,HBvsMB)
write.table(DEgenesHBvsMB, "~/DEgenesHBvsMB001.txt", sep="\t")

# 5. Make heatmap of OT, MB and HB using genes from OT vs HB (batch corrected counts with limma remove batch effects see section 1):
# filter for genes that are found in OT vs HB (p<0.05)
filter1=read.csv("~/OT_MD_HB/DEgenesOTvsHB005.csv")

# Filter in the dataset for only genes of interest:
dat=assay(vsd)
dat=dat[is.element(rownames(dat),as.character(filter1[,1])),]
dat=data.frame(Ensembl_ID=rownames(dat),dat)
dat=merge(filter1[,1:2],dat,by.x=1,by.y=1)
dat <- na.omit(dat) # remove the NA one
rownames(dat)=as.character(dat[,2])
dat=dat[,c(3:dim(dat)[2])]
dat=t(scale(t(dat)))

col_map = colorRamp2(c(-2, 1, 0, 1, 2), c("#313695", "#F88D51", "#E0F3F8", "74ADD1", "#BE1826"))
Heatmap(dat, name = "Expression", column_order = colnames(dat), col = col_map, column_title = "Up and downregulated", row_km = 2, border = TRUE, width = unit(8, "cm"), height = unit(8, "cm"))

# End




