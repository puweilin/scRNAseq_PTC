library(plyr)
library(dplyr)
library(sctransform)
library(Seurat)
library(ggplot2)
library(ggsci)
library(readr)
library(readxl)


#Identify the top variable genes in the thyroid cells 
Thyroid = FindVariableFeatures(Thyroid, nfeatures = 2000)
remove.ribo = Thyroid@assays$RNA@var.features[grep("^(RPL|RPS)+", Thyroid@assays$RNA@var.features)]
remove.mito = Thyroid@assays$RNA@var.features[grep("(MT-)+", Thyroid@assays$RNA@var.features)]
Topfeatures = setdiff(Thyroid@assays$RNA@var.features, c(remove.ribo, remove.mito))

#Read the bulk-rna sequencing dataset and its clinical information 
Bulk_RNA = read_rds("../Bulk_RNA/Bulk.rds")
Clinical = read_excel("../Bulk_RNA/Samplesheet.xlsx")
Clinical$Type = mapvalues(Clinical$Type, from = c("FA", "FA-N", "FTC", "FTC-N"), 
                          to = c("Tumor", "Control", "Tumor", "Control"))

#Subset the Bulk datasets with the top variable genes in the thyroid cells; 
shared.features = intersect(Topfeatures, rownames(Bulk_RNA))
Bulk_RNA_sub = Bulk_RNA[match(Topfeatures, rownames(Bulk_RNA)), ]

Korean = Bulk_RNA_sub[, which(Clinical$Source == "Korean")]
TCGA = Bulk_RNA_sub[, which(Clinical$Source == "TCGA")]
TCGA_type = Clinical$Type[match(colnames(TCGA), Clinical$NewID)]

#Subset the scRNA datasets;
dat = as.data.frame(Thyroid@assays$RNA@data[match(shared.features, rownames(Thyroid)), ])


#KNN method was implemented (K = 9) 
#To classify the cells according to the status of its most correlated sample in Bulk RNA-seq dataset;

classifier = function(cell, Bulk, identity){
  cor_value = apply(Bulk, 2, function(x) cor(cell, x, method = "spearman"))
  max9 = identity[order(cor_value, decreasing = T)[1:9]]
  return(names(sort(table(max9), decreasing = T)[1]))
}

#Verify the classifier using TCGA as Bulk and Korean samples as test
Korean_pred = apply(Korean, 2, function(x) classifier(x, TCGA, identity = Clinical$Type[match(colnames(TCGA), Clinical$NewID)]))
Korean_pred = as.data.frame.array(table(Korean_pred, Clinical$Type[match(colnames(Korean), Clinical$NewID)]))


#predict the malignant/non-malignant cells of scRNA data
library(doParallel)
library(snow)
library(pbapply)


n = dim(dat)[2]
f <- function(){
  pb <- txtProgressBar(min=1, max=n-1,style=3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count)
    Sys.sleep(0.01)
    flush.console()
    c(...)
  }
}

cl2 <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl2)
pred = foreach(i = 1:dim(dat)[2], .combine = f()) %dopar% {
  classifier(dat[,i], TCGA, identity = TCGA_type)
}
scRNA_pred = as.data.frame.array(prop.table(table(pred, Thyroid$Type), margin = 2))
