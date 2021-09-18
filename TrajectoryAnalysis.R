library(plyr)
library(dplyr)
library(sctransform)
library(Seurat)
library(ggplot2)
library(ggsci)
library(readr)
library(readxl)
library(monocle)


#Trajectory analysis using all thyroid cells; 

  #Construct the monocle object
  dat = GetAssayData(THCA, slot = "counts", assay = "RNA")
  dat = as.matrix(dat)
  gene_annotation = data.frame(gene_short_name = rownames(dat))
  rownames(gene_annotation) = gene_annotation$gene_short_name
  phenotype = THCA@meta.data
  sample_sheet = phenotype[, c("orig.ident", "RNA_snn_res.0.6","Origin")]
  gd = new("AnnotatedDataFrame", data = gene_annotation)
  pd = new("AnnotatedDataFrame", data = sample_sheet)
  THCA_monocle = newCellDataSet(as(dat, "sparseMatrix"), phenoData = pd, 
                                  featureData = gd, lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())

  
  #Estimate size factors and dispersions 
  THCA_monocle = estimateSizeFactors(THCA_monocle)
  options(DelayedArray.auto.block.size = 1000e6)
  THCA_monocle = estimateDispersions(THCA_monocle)
  
  
  #constructing single cell trajectories
  ordering_genes = subset(disp_table, mean_expression >= 0.5 &
                            dispersion_empirical >=1 * dispersion_fit)$gene_id
  remove.idx1 = grep("^(RPL|RPS)+", ordering_genes)
  remove.idx2 = grep("(MT-)+", ordering_genes)
  ordering_genes = ordering_genes[-c(remove.idx1, remove.idx2)]
  my_object = setOrderingFilter(my_object, ordering_genes = ordering_genes)
  my_object = reduceDimension(my_object, method = "DDRTree", max_components = 2, verbose = T)
  print("Finished Dimension Reduction")
  my_object = orderCells(my_object)


  
  
  
# Trajectory analysis using two types of normal thyroid cells; 
  dat = GetAssayData(THCA_Control, slot = "counts", assay = "RNA")
  dat = as.matrix(dat)
  gene_annotation = data.frame(gene_short_name = rownames(dat))
  rownames(gene_annotation) = gene_annotation$gene_short_name
  phenotype = THCA_Control@meta.data
  sample_sheet = phenotype[, c("orig.ident", "RNA_snn_res.0.1","Origin")]
  gd = new("AnnotatedDataFrame", data = gene_annotation)
  pd = new("AnnotatedDataFrame", data = sample_sheet)
  monocle_con = newCellDataSet(as(dat, "sparseMatrix"), phenoData = pd, 
                                featureData = gd, lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  
  monocle_con = estimateSizeFactors(monocle_con)
  options(DelayedArray.auto.block.size = 1000e6)
  monocle_con = estimateDispersions(monocle_con)
  
  disp_table = dispersionTable(monocle_con)
  ordering_genes = subset(disp_table, mean_expression >= 1 &
                            dispersion_empirical >=1 * dispersion_fit)$gene_id
  monocle_con = setOrderingFilter(monocle_con, ordering_genes = ordering_genes)
  monocle_con = reduceDimension(monocle_con, reduction_method = "DDRTree",verbose = T, max_components = 2)
  print("Finished Dimension Reduction")
  monocle_con = orderCells(monocle_con)