
###METABRIC networks
##metabric data can be accessed from the Synapse database (syn2160410), upon request
metabric<-read.csv("path_to_metabric_data")
meta<-read.csv("path_to_metabric_metadata")

metabric_Basal<-metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]
metabric_LumA<-metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumA"]

cor_basal<-cor(t(metabric_Basal))
cor_lumA<-cor(t(metabric_LumA))

write.csv(cor_basal, "cor_basal_METABRIC.csv")
write.csv(cor_lumA, "cor_lumA_METABRIC.csv")


###TCGA networks
library(TCGAbiolinks)
load("TCGA_BRCA_tumor.RData")
colnames(TCGA_tumor)<-gsub('([^-]*)-([^-]*)-([^-]*)-.*', '\\1-\\2-\\3', colnames(TCGA_tumor))
subtype<-TCGAquery_subtype(tumor = "brca")

##check that TCGA id is in subtype
TCGA_tumor<-TCGA_tumor[,which(colnames(TCGA_tumor) %in% c(subtype[,1])$patient)]

##gene filtering
TCGA_tumor<-TCGA_tumor[which(rowSums(TCGA_tumor<1)<=1040),]
##log transforming
TCGA_tumor<-log2(TCGA_tumor+1)

TCGA_basal<-TCGA_tumor[,subtype$BRCA_Subtype_PAM50[match(colnames(TCGA_tumor), c(subtype[,1])$patient)]=="Basal"]
TCGA_lumA<-TCGA_tumor[,subtype$BRCA_Subtype_PAM50[match(colnames(TCGA_tumor), c(subtype[,1])$patient)]=="LumA"]

cor_basal_TCGA<-cor(t(TCGA_basal))
cor_lumA_TCGA<-cor(t(TCGA_lumA))

write.csv(cor_basal_TCGA, "cor_basal_TCGA.csv")
write.csv(cor_lumA_TCGA, "cor_lumA_TCGA.csv")