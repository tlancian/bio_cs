library(clusterProfiler)
library(openxlsx)
library(enrichplot)
require(DOSE)

##########
###GO enrichment of differential subraphs comparing proteome and transcriptome
#########

setwd("path_to_subgraphs")
files<-list.files(recursive = T)
files<-files[-c(1:2)]
background<-read.csv("transcr_atlas.txt", sep="\t")[,1]

module<-read.csv(files[1], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_CPTAC_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

module<-read.csv(files[11], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_TCGA_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

for(i in 1:length(files)){
  module<-read.csv(files[i], sep="\t")[,1]
  ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)
  write.xlsx(ego, gsub(".txt", "_GO.xlsx",files[i]))
}

##########
###GO enrichment of differential subraphs comparing cell lines
#########

setwd("path_to_subgraphs")
files<-list.files(recursive = T)
files<-files[-c(1,62)]
background<-read.csv("properseq_atlas.txt", sep="\t")[,1]

for(i in 1:length(files)){
  tryCatch({
  module<-read.csv(files[i], sep="\t")[,1]
  ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)
  write.xlsx(ego, gsub(".txt", "_GO.xlsx",files[i]), overwrite = T)
  }, error=function(e){print("no lines")})
}


module<-read.csv(files[2], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_HEK_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

##########
###GO enrichment of differential subraphs comparing breast cancer subtypes, TCGA
#########

setwd("path_to_subgraphs")
files<-list.files(recursive = T)
files<-files[-c(1,2)]
background<-read.csv("tcga_atlas.txt", sep="\t")[,1]


module<-read.csv(files[1], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_tcga_basal_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

module<-read.csv(files[11], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_tcga_lumA_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()


##########
###GO enrichment of differential subraphs comparing breast cancer subtypes, METABRIC
#########

setwd("path_to_subgraphs")

files<-list.files(recursive = T)
files<-files[-c(1)]
background<-read.csv("metabric_atlas.txt", sep="\t")[,1]


module<-read.csv(files[1], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_metabric_basal_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

module<-read.csv(files[11], sep="\t")[,1]
ego<-enrichGO(module,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_metabric_lumA_1.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()


###########
##Compare correlation and Proportionality
###########

#####TCGA

path_to_prop_subgraphs<-"../prop_tcga"
path_to_corr_subgraphs<-"../tcga"

background<-read.csv(paste(path_to_prop_subgraphs,"tcga_atlas.txt", sep="/"), sep="\t")[,1]
TCGA_basal_orig<-read.csv(paste(path_to_corr_subgraphs, "tcga_basal_full_10/sol_1.txt", sep="/"), sep="\t")[,1]
TCGA_basal_prop<-read.csv(paste(path_to_prop_subgraphs, "prop_tcga_basal_full_10/sol_1.txt", sep="/"), sep="\t")[,1]

a<-length(intersect(TCGA_basal_orig,TCGA_basal_prop))
b<-length(setdiff(TCGA_basal_orig,TCGA_basal_prop))
c<-length(setdiff(TCGA_basal_prop, TCGA_basal_orig))
d<-length(setdiff(background, union(TCGA_basal_orig,TCGA_basal_prop)))

a/length(union(TCGA_basal_orig,TCGA_basal_prop)) #Jaccard
fisher.test(matrix(c(a,b,c,d), nrow = 2)) #Fisher test
mat<-matrix(c(a,b,c,d), nrow = 2, dimnames=list(c("In prop", "Not in prop"), c("In corr", "Not in corr"))) #contingency matrix
write.csv(mat, "contingency_TCGA_basal.csv")

ego<-enrichGO(TCGA_basal_prop,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_TCGA_basal_prop.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

background<-read.csv(paste(path_to_prop_subgraphs,"tcga_atlas.txt", sep="/"), sep="\t")[,1]
TCGA_lumA_orig<-read.csv(paste(path_to_corr_subgraphs, "tcga_lumA_full_10/sol_1.txt", sep="/"), sep="\t")[,1]
TCGA_lumA_prop<-read.csv(paste(path_to_prop_subgraphs, "prop_tcga_lumA_full_10/sol_1.txt", sep="/"), sep="\t")[,1]

a<-length(intersect(TCGA_lumA_orig,TCGA_lumA_prop))
b<-length(setdiff(TCGA_lumA_orig,TCGA_lumA_prop))
c<-length(setdiff(TCGA_lumA_prop, TCGA_lumA_orig))
d<-length(setdiff(background, union(TCGA_lumA_orig,TCGA_lumA_prop)))

a/length(union(TCGA_lumA_orig,TCGA_lumA_prop))
fisher.test(matrix(c(a,b,c,d), nrow = 2))

mat<-matrix(c(a,b,c,d), nrow = 2, dimnames=list(c("In prop", "Not in prop"), c("In corr", "Not in corr")))
write.csv(mat, "../contingency_TCGA_lumA.csv")

ego<-enrichGO(TCGA_lumA_prop,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_TCGA_lumA_prop.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

### METABRIC

path_to_prop_subgraphs<-"../prop_metabric"
path_to_corr_subgraphs<-"../metabric"

background<-read.csv(paste(path_to_prop_subgraphs,"metabric_atlas.txt", sep="/"), sep="\t")[,1]
METABRIC_basal_orig<-read.csv(paste(path_to_corr_subgraphs, "metabric_basal_full_10/sol_1.txt", sep="/"), sep="\t")[,1]
METABRIC_basal_prop<-read.csv(paste(path_to_prop_subgraphs, "prop_metabric_basal_full_10/sol_1.txt", sep="/"), sep="\t")[,1]

a<-length(intersect(METABRIC_basal_orig,METABRIC_basal_prop))
b<-length(setdiff(METABRIC_basal_orig,METABRIC_basal_prop))
c<-length(setdiff(METABRIC_basal_prop, METABRIC_basal_orig))
d<-length(setdiff(background, union(METABRIC_basal_orig,METABRIC_basal_prop)))

a/length(union(METABRIC_basal_orig,METABRIC_basal_prop))
fisher.test(matrix(c(a,b,c,d), nrow = 2))

mat<-matrix(c(a,b,c,d), nrow = 2, dimnames=list(c("In prop", "Not in prop"), c("In corr", "Not in corr")))
write.csv(mat, "contingency_Metabric_basal.csv")

ego<-enrichGO(METABRIC_basal_prop,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_Metabric_basal_prop.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

background<-read.csv(paste(path_to_prop_subgraphs,"metabric_atlas.txt", sep="/"), sep="\t")[,1]
METABRIC_lumA_orig<-read.csv(paste(path_to_corr_subgraphs, "metabric_lumA_full_10/sol_1.txt", sep="/"), sep="\t")[,1]
METABRIC_lumA_prop<-read.csv(paste(path_to_prop_subgraphs, "prop_metabric_lumA_full_10/sol_1.txt", sep="/"), sep="\t")[,1]

a<-length(intersect(METABRIC_lumA_orig,METABRIC_lumA_prop))
b<-length(setdiff(METABRIC_lumA_orig,METABRIC_lumA_prop))
c<-length(setdiff(METABRIC_lumA_prop, METABRIC_lumA_orig))
d<-length(setdiff(background, union(METABRIC_lumA_orig,METABRIC_lumA_prop)))

a/length(union(METABRIC_lumA_orig,METABRIC_lumA_prop))
fisher.test(matrix(c(a,b,c,d), nrow = 2))

mat<-matrix(c(a,b,c,d), nrow = 2, dimnames=list(c("In prop", "Not in prop"), c("In corr", "Not in corr")))
write.csv(mat, "../contingency_Metabric_lumA.csv")

ego<-enrichGO(METABRIC_lumA_prop,OrgDb = 'org.Hs.eg.db',keyType = "SYMBOL", ont="BP", universe = background)

pdf("GO_Metabric_lumA_prop.pdf", 7 ,7)
dotplot(ego, showCategory=10)
dev.off()

