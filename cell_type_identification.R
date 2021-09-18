library(plyr)
library(dplyr)
library(sctransform)
library(Seurat)
library(ggplot2)
library(ggsci)
library(readr)
library(readxl)
library(viridis)

q =  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(legend.position=c(0.15,0.95))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 16, face = "bold"))+
  theme(axis.text.x=element_text(size=25, face = "bold"))+
  theme(axis.title.x=element_text(size=25, face = "bold"))+
  theme(axis.text.y=element_text(size=20, face = "bold"))+
  theme(axis.title.y=element_text(size=25, face = "bold"))+
  theme(plot.title = element_text(hjust = 0.5, vjust=0))

#Load the batch effect corrected data to the slot "data" of THCA_MNN object
remove.ribo = THCA_MNN@assays$RNA@var.features[grep("^(RPL|RPS)+", THCA_MNN@assays$RNA@var.features)]
remove.mito = THCA_MNN@assays$RNA@var.features[grep("(MT-)+", THCA_MNN@assays$RNA@var.features)]
selected.features = setdiff(rownames(THCA_MNN@assays$RNA@scale.data), c(remove.ribo, remove.mito))

THCA_MNN = RunPCA(THCA_MNN, npcs = 50, features = selected.features )
THCA_MNN = RunUMAP(THCA_MNN, reduction = "pca", dims = 1:15)

THCA_MNN = FindNeighbors(THCA_MNN, reduction = "pca", dims = 1:15)
THCA_MNN = FindClusters(THCA_MNN, resolution = 1)

Diff_genes = FindAllMarkers(THCA_MNN, logfc.threshold = 0.25)

## Identify the cell type based on the DEGs of each cluster
# Cluster 0: IL7R, CCR7, CD3D, CD3E, CD3G -- T cells 
# Cluster 1: CD79A, CD74, IGHM, MS4A1, CD19 -- B cells
# Cluster 2: GZMK, GZMA, CD8A, NKG7 --NKT cells
# Cluster 3: TIGIT, TNFRSF4, CXCL13, IL32, TNFRSF18, FOXP3, CTLA4 -- Treg cells
# Cluster 4: GZMK, GZMA, CD8A, NKG7 --NKT cells
# Cluster 5: FN1, KRT19 -- Epithelial cells (thyroid cells)
# Cluster 6: GZMK, NKG7 -- NK cells 
# Cluster 7: TG, CITED1 -- Epithelial cells (thyroid cells)
# Cluster 8: CD14 -- Macrophages
# Cluster 9: IL7R, CCR7, CD3D, CD3E, CD3G -- T cells 
# Cluster 10: IL7R, CCR7, CD3D, CD3E, CD3G -- T cells 
# Cluster 11: TG, CITED1, -- Epithelial cells (thyroid cells)
# Cluster 12: TG, TSHR -- Epithelial cells (thyroid cells)
# Cluster 13: FN1, KRT19 -- Thyroid cells 
# Cluster 14: CD14 -- Macrophages 
# Cluster 15: TAGLN, ACTA2 -- Fibroblasts
# Cluster 16: FN1, KRT19 -- Thyroid cells 
# Cluster 17: FN1, KRT19 -- Thyroid cells 
# Cluster 18: SPARCL1, PLVAP, VWF -- Endothelial cells 
# Cluster 19: FN1, KRT19 -- Thyroid cells 
# Cluster 20: TG, TSHR -- Epitheial cells  
# Cluster 21: MMRN1, LYVE1 -- Lymphatic Endothelial cells 
# Cluster 22: TG, TSHR -- Epitheial cells 
# Cluster 23: DCN,SPARC -- Fibroblasts
# Cluster 24: MS4A1, CD79A, JCHAIN -- B cells 
# Cluster 25：CITED1, KRT19 -- thyroid cells



new_clusters = c("T & NK cells", "B cells", "T & NK cells", "T & NK cells", "T & NK cells", "Thyroid cells", "T & NK cells", "Thyroid cells", "Myeloid cells", "T & NK cells", "T & NK cells",  "Thyroid cells", "Thyroid cells","Thyroid cells", "Myeloid cells", "Fibroblasts", "Thyroid cells", "Thyroid cells", "Endothelial cells", "Thyroid cells", "Thyroid cells",
                 "Endothelial cells",  "Thyroid cells", "Fibroblasts", "B cells", "Thyroid cells")


THCA_MNN@active.ident = mapvalues(THCA_MNN@active.ident, from = seq(0, 25, 1), to = new_clusters)
THCA_MNN = AddMetaData(THCA_MNN, metadata = Idents(THCA_MNN), col.name = "manual_identity")

#UMAP visualization
p = Embedding_plot(THCA_MNN,  vis_type = "umap", group = "manual_identity", sort_type = FALSE, legend_distance = 3, text_size = 6, user_cols = My_colors[c(1:9, 11)], point_size = 0.6, point_alpha = 0.9)

