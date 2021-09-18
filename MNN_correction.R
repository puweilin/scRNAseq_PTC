library(plyr)
library(dplyr)
library(sctransform)
library(Seurat)
library(readr)
library(readxl)

source("d:/ScientificProject/THCA_scSeq/ReferenceData/Immune_cell_identification.R")
source("d:/ScientificProject/THCA_scSeq/NewAnalysis/Code/FastFindMarkers.R")
source("d:/ScientificProject/THCA_scSeq/NewAnalysis/Code/Embedding_plot.R")
source("d:/ScientificProject/THCA_scSeq/NewAnalysis/Code/scEnrichment.R")
source("d:/ScientificProject/THCA_scSeq/NewAnalysis/Code/batch_effect_removal.R")

#Load the corrected dataset using MNN to the Seurat object; 
load_corrected_csv = function(seurat_object, object_dir, batch_order = c()){
  corrected_dat = read_csv(paste(object_dir, "X.csv", sep = "/"), col_names = F)
  corrected_dat = as.matrix(t(corrected_dat))
  cellnames = c()
  scaled.data = seurat_object@assays$RNA@scale.data
  if(is.vector(batch_order)){
    for(i in 1:length(batch_order)){
      batchname = strsplit(batch_order[i], split = "_before")[[1]][1]
      cellnames = append(cellnames, colnames(scaled.data)[which(seurat_object$orig.ident == batchname)])
    }
  }else{
    batchname = as.character(sort(unique(origial_idents)))
    for(i in 1:length(as.character(sort(unique(origial_idents))))){
      cellnames = append(cellnames, colnames(scaled.data)[which(seurat_object$orig.ident == batchname[i])])
    }
  }
  colnames(corrected_dat) = cellnames
  rownames(corrected_dat) = rownames(scaled.data)
  
  #reformat
  corrected_dat = corrected_dat[, match(colnames(scaled.data), colnames(corrected_dat))]
  return(corrected_dat)
}

batch_order = c("LN10r_beforeCorrection.csv", "LN11r_beforeCorrection.csv", "LN2l_beforeCorrection.csv",
                "LN3l_beforeCorrection.csv", "LN3r_beforeCorrection.csv", "LN5r_beforeCorrection.csv",
                "LN6r_beforeCorrection.csv", "LN7r_beforeCorrection.csv", "P1_beforeCorrection.csv", 
                "P2_beforeCorrection.csv", "P3_beforeCorrection.csv", "P5_beforeCorrection.csv", 
                "P8_beforeCorrection.csv", "P9_beforeCorrection.csv", "SC11_beforeCorrection.csv", 
                "SC4_beforeCorrection.csv", "T10_beforeCorrection.csv", "T1_beforeCorrection.csv",
                "T2_beforeCorrection.csv", "T3_beforeCorrection.csv", "T5_beforeCorrection.csv", 
                "T8_beforeCorrection.csv", "T9_beforeCorrection.csv")

corrected_dat = load_corrected_csv(THCA, 
                                   object_dir = "d:/ScientificProject/THCA_scSeq/NewAnalysis/Sample_Integration/AfterMNN/output/", 
                                   batch_order = batch_order)


#Replace the scaled.data slot for further analysis;
THCA@assays$RNA@scale.data = corrected_dat
THCA = RunPCA(THCA)
THCA = RunUMAP(THCA, reduction = "pca", dims = 1:20)
THCA = FindNeighbors(THCA, reduction = "pca", dims = 1:20)
THCA = FindClusters(THCA, resolution = 1)


