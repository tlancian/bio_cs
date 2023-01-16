library(openxlsx)
library(TCGAbiolinks)

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

######################################
#####CPTAC and TCGA data filtering
######################################

#supplementary table 1 from Mertins et al. doi: 10.1038/nature18003
meta<-read.xlsx("path_to_mertin_supptable1",1)
rownames(meta)<-meta[,2]

#supplementary table 3 from Mertins et al. doi: 10.1038/nature18003
CPTAC_BRCA<-read.csv("path_to_mertin_supptable3", sep=";", dec=",")
rownames(CPTAC_BRCA)<-CPTAC_BRCA[,1]
##remove replicates
CPTAC_BRCA<-CPTAC_BRCA[,-c(93:95)]
#change colnames
colnames(CPTAC_BRCA)[13:92]<-gsub("[.]", "-", colnames(CPTAC_BRCA)[13:92])
colnames(CPTAC_BRCA)[13:92]<-gsub("TCGA", "", colnames(CPTAC_BRCA)[13:92])
colnames(CPTAC_BRCA)[13:92]<-gsub('([^-]*)-([^-]*)-.*', '\\1-\\2', colnames(CPTAC_BRCA)[13:92])

###TCGA
load("TCGA_BRCA_tumor.RData")
colnames(TCGA_tumor)<-gsub('([^-]*)-([^-]*)-([^-]*)-.*', '\\1-\\2-\\3', colnames(TCGA_tumor))
subtype<-TCGAquery_subtype(tumor = "brca")
colnames(TCGA_tumor)<-gsub("TCGA-", "", colnames(TCGA_tumor))

##gene filtering
TCGA_tumor<-TCGA_tumor[which(rowSums(TCGA_tumor<1)<=1040),]
#log transformation
TCGA_tumor<-log2(TCGA_tumor+1)

##kept only samples with both proteomic and transcriptomic profiling
inboth<-intersect(colnames(CPTAC_BRCA), colnames(TCGA_tumor))
CPTAC_BRCA_sel<-cbind.data.frame(CPTAC_BRCA[,11],CPTAC_BRCA[,inboth])
TCGA_sel<-TCGA_tumor[,inboth]


CPTAC_BRCA_sel<-changenames(CPTAC_BRCA_sel[,-1], anno=cbind(rownames(CPTAC_BRCA_sel), CPTAC_BRCA_sel[,1]))

##keep transcripts with less than 20 NAs
CPTAC_BRCA_sel<-CPTAC_BRCA_sel[which(rowSums(is.na(CPTAC_BRCA_sel))<20),]

######################################
##### compute correlations
#######################################

cor_TCGA<-cor(t(TCGA_sel), use="pairwise.complete.obs")
write.csv(cor_TCGA, file="cor_TCGA.csv")

cor_CPTAC<-cor(t(CPTAC_BRCA_sel), use="pairwise.complete.obs")
write.csv(cor_CPTAC, file="cor_CPTAC.csv")

