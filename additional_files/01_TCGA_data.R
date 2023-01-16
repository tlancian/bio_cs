##Data from TCGAbiolinks, mRNA FPKM
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

##function to change gene ids, keeping the probe with the highest expression for each gene
changenames<-function(data, anno){
  annotation_sel=anno[match( rownames(data), anno[,1]),2]
  
  if(length(which(annotation_sel==""))>0){
    data<-data[-which(annotation_sel==""),]
    annotation_sel<-annotation_sel[-which(annotation_sel=="")]
  }
  
  a<-which(duplicated(annotation_sel))
  while(length(a)>0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
        m=which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
        data=data[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
        annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
      }
    }
    
    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }
  
  rownames(data)=annotation_sel
  return(data)
}

# mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
query.exp.hg38 <- GDCquery(project = "TCGA-BRCA", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "exp.rda")
save(expdat, file="TCGA_BRCA_FPKM.RData")

load("D:/Drive/PC lab/Documents/R analyses/TCGA biolinks/exp.rda")

TCGA<-assay(expdat)

##convert ENSEMBL ID in GENE SYMBOL
library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
genes_id<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                  values = rownames(TCGA), 
                  mart = ensembl)

##Probes merged by taking the one with the highest average expression for each gene symbol
TCGA<-changenames(TCGA, anno=genes_id)

##only tumor, normal, metastatic
TCGA_tumor<-TCGA[,expdat$sample_type_id=="01"]
TCGA_norm<-TCGA[,expdat$sample_type_id=="11"]
TCGA_meta<-TCGA[,expdat$sample_type_id=="06"]
colnames(TCGA_tumor)<-gsub('([^-]*)-([^-]*)-([^-]*)-([^-]*)-([^-]*)-([^-]*)-.*', '\\1-\\2-\\3-\\4', colnames(TCGA_tumor))
colnames(TCGA_norm)<-gsub('([^-]*)-([^-]*)-([^-]*)-([^-]*)-([^-]*)-([^-]*)-.*', '\\1-\\2-\\3-\\4', colnames(TCGA_norm))
colnames(TCGA_meta)<-gsub('([^-]*)-([^-]*)-([^-]*)-([^-]*)-([^-]*)-([^-]*)-.*', '\\1-\\2-\\3-\\4', colnames(TCGA_meta))

save(TCGA_tumor, file="TCGA_BRCA_tumor.RData")
save(TCGA_norm, file="TCGA_BRCA_norm.RData")
save(TCGA_meta, file="TCGA_BRCA_meta.RData")
