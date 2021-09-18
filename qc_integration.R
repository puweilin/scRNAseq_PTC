library(plyr)
library(dplyr)
library(sctransform)
library(Seurat)
library(ggplot2)
library(ggsci)
library(readr)
library(readxl)
library(DoubletFinder)

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


#T1
T1 = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "T1" )
T1[["percent.mt"]] = PercentageFeatureSet(object = T1, pattern = "^MT-")
VlnPlot(T1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T1 = subset(T1, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

T1 = NormalizeData(T1)
T1 = ScaleData(T1)
T1 = FindVariableFeatures(T1, selection.method = "vst", nfeatures = 10000)
T1 = RunPCA(T1, features = T1@assays$RNA@var.features)
T1 =  FindNeighbors(T1, reduction = "pca", dims = 1:50)
T1 = FindClusters(T1, resolution = 0.6)

sweep.data = paramSweep_v3(T1, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T1@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T1$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T1=doubletFinder_v3(T1, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T1@meta.data[,8]
T1 = SubsetData(T1, cells = colnames(T1)[which(doubletsID == "Singlet")])
T1 = NormalizeData(T1)
T1 = ScaleData(T1)
T1 = FindVariableFeatures(T1, selection.method = "vst", nfeatures = 10000)
T1 = RunPCA(T1, features = T1@assays$RNA@var.features)
T1 =  FindNeighbors(T1, reduction = "pca", dims = 1:50)
T1 = FindClusters(T1, resolution = 0.6)



#P1
P1 = CreateSeuratObject(counts = P1.data, min.cells = 3, min.features = 200, project = "P1" )
P1[["percent.mt"]] = PercentageFeatureSet(object = P1, pattern = "^MT-")
VlnPlot(P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P1 = subset(P1, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
P1 = NormalizeData(P1)
P1 = ScaleData(P1)
P1 = FindVariableFeatures(P1, selection.method = "vst", nfeatures = 10000)
P1 = RunPCA(P1, features = P1@assays$RNA@var.features)
P1 =  FindNeighbors(P1, reduction = "pca", dims = 1:50)
P1 = FindClusters(P1, resolution = 0.6)

sweep.data = paramSweep_v3(P1, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(P1@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(P1$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
P1=doubletFinder_v3(P1, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = P1@meta.data[,8]
P1 = SubsetData(P1, cells = colnames(P1)[which(doubletsID == "Singlet")])

P1 = NormalizeData(P1)
P1 = ScaleData(P1)
P1 = FindVariableFeatures(P1, selection.method = "vst", nfeatures = 10000)
P1 = RunPCA(P1, features = P1@assays$RNA@var.features)
P1 =  FindNeighbors(P1, reduction = "pca", dims = 1:50)
P1 = FindClusters(P1, resolution = 0.6)


#T2
T2 = CreateSeuratObject(counts = T2.data, min.cells = 3, min.features = 200, project = "T2" )
T2[["percent.mt"]] = PercentageFeatureSet(object = T2, pattern = "^MT-")
VlnPlot(T2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T2 = subset(T2, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
T2 = NormalizeData(T2)
T2 = ScaleData(T2)
T2 = FindVariableFeatures(T2, selection.method = "vst", nfeatures = 10000)
T2 = RunPCA(T2, features = T2@assays$RNA@var.features)
T2 =  FindNeighbors(T2, reduction = "pca", dims = 1:50)
T2 = FindClusters(T2, resolution = 0.6)

sweep.data = paramSweep_v3(T2, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T2@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T2$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T2=doubletFinder_v3(T2, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T2@meta.data[,8]
T2 = SubsetData(T2, cells = colnames(T2)[which(doubletsID == "Singlet")])

T2 = NormalizeData(T2)
T2 = ScaleData(T2)
T2 = FindVariableFeatures(T2, selection.method = "vst", nfeatures = 10000)
T2 = RunPCA(T2, features = T2@assays$RNA@var.features)
T2 =  FindNeighbors(T2, reduction = "pca", dims = 1:50)
T2 = FindClusters(T2, resolution = 0.6)


#P2
P2 = CreateSeuratObject(counts = P2.data, min.cells = 3, min.features = 200, project = "P2" )
P2[["percent.mt"]] = PercentageFeatureSet(object = P2, pattern = "^MT-")
VlnPlot(P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P2 = subset(P2, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
P2 = NormalizeData(P2)
P2 = ScaleData(P2)
P2 = FindVariableFeatures(P2, selection.method = "vst", nfeatures = 10000)
P2 = RunPCA(P2, features = P2@assays$RNA@var.features)
P2 =  FindNeighbors(P2, reduction = "pca", dims = 1:50)
P2 = FindClusters(P2, resolution = 0.6)

sweep.data = paramSweep_v3(P2, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(P2@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(P2$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
P2=doubletFinder_v3(P2, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = P2@meta.data[,8]
P2 = SubsetData(P2, cells = colnames(P2)[which(doubletsID == "Singlet")])
P2 = NormalizeData(P2)
P2 = ScaleData(P2)
P2 = FindVariableFeatures(P2, selection.method = "vst", nfeatures = 10000)
P2 = RunPCA(P2, features = P2@assays$RNA@var.features)
P2 =  FindNeighbors(P2, reduction = "pca", dims = 1:50)
P2 = FindClusters(P2, resolution = 0.6)



#LN2l
LN2l = CreateSeuratObject(counts = LN2l.data, min.cells = 3, min.features = 200, project = "LN2l" )
LN2l[["percent.mt"]] = PercentageFeatureSet(object = LN2l, pattern = "^MT-")
VlnPlot(LN2l, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN2l = subset(LN2l, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN2l = NormalizeData(LN2l)
LN2l = ScaleData(LN2l)
LN2l = FindVariableFeatures(LN2l, selection.method = "vst", nfeatures = 10000)
LN2l = RunPCA(LN2l, features = LN2l@assays$RNA@var.features)
LN2l =  FindNeighbors(LN2l, reduction = "pca", dims = 1:50)
LN2l = FindClusters(LN2l, resolution = 0.6)

sweep.data = paramSweep_v3(LN2l, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN2l@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN2l$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN2l=doubletFinder_v3(LN2l, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN2l@meta.data[,8]
LN2l = SubsetData(LN2l, cells = colnames(LN2l)[which(doubletsID == "Singlet")])

LN2l = NormalizeData(LN2l)
LN2l = ScaleData(LN2l)
LN2l = FindVariableFeatures(LN2l, selection.method = "vst", nfeatures = 10000)
LN2l = RunPCA(LN2l, features = LN2l@assays$RNA@var.features)
LN2l =  FindNeighbors(LN2l, reduction = "pca", dims = 1:50)
LN2l = FindClusters(LN2l, resolution = 0.6)




#T3
T3 = CreateSeuratObject(counts = T3.data, min.cells = 3, min.features = 200, project = "T3" )
T3[["percent.mt"]] = PercentageFeatureSet(object = T3, pattern = "^MT-")
VlnPlot(T3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T3 = subset(T3, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
T3 = NormalizeData(T3)
T3 = ScaleData(T3)
T3 = FindVariableFeatures(T3, selection.method = "vst", nfeatures = 10000)
T3 = RunPCA(T3, features = T3@assays$RNA@var.features)
T3 =  FindNeighbors(T3, reduction = "pca", dims = 1:50)
T3 = FindClusters(T3, resolution = 0.6)

sweep.data = paramSweep_v3(T3, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T3@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T3$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T3=doubletFinder_v3(T3, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T3@meta.data[,8]
T3 = SubsetData(T3, cells = colnames(T3)[which(doubletsID == "Singlet")])

T3 = NormalizeData(T3)
T3 = ScaleData(T3)
T3 = FindVariableFeatures(T3, selection.method = "vst", nfeatures = 10000)
T3 = RunPCA(T3, features = T3@assays$RNA@var.features)
T3 =  FindNeighbors(T3, reduction = "pca", dims = 1:50)
T3 = FindClusters(T3, resolution = 0.6)



#P3
P3 = CreateSeuratObject(counts = P3.data, min.cells = 3, min.features = 200, project = "P3" )
P3[["percent.mt"]] = PercentageFeatureSet(object = P3, pattern = "^MT-")
VlnPlot(P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P3 = subset(P3, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
P3 = NormalizeData(P3)
P3 = ScaleData(P3)
P3 = FindVariableFeatures(P3, selection.method = "vst", nfeatures = 10000)
P3 = RunPCA(P3, features = P3@assays$RNA@var.features)
P3 =  FindNeighbors(P3, reduction = "pca", dims = 1:50)
P3 = FindClusters(P3, resolution = 0.6)

sweep.data = paramSweep_v3(P3, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(P3@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(P3$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
P3=doubletFinder_v3(P3, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = P3@meta.data[,8]
P3 = SubsetData(P3, cells = colnames(P3)[which(doubletsID == "Singlet")])

P3 = NormalizeData(P3)
P3 = ScaleData(P3)
P3 = FindVariableFeatures(P3, selection.method = "vst", nfeatures = 10000)
P3 = RunPCA(P3, features = P3@assays$RNA@var.features)
P3 =  FindNeighbors(P3, reduction = "pca", dims = 1:50)
P3 = FindClusters(P3, resolution = 0.6)


#LN3l
LN3l = CreateSeuratObject(counts = LN3l.data, min.cells = 3, min.features = 200, project = "LN3l" )
LN3l[["percent.mt"]] = PercentageFeatureSet(object = LN3l, pattern = "^MT-")
VlnPlot(LN3l, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN3l = subset(LN3l, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN3l = NormalizeData(LN3l)
LN3l = ScaleData(LN3l)
LN3l = FindVariableFeatures(LN3l, selection.method = "vst", nfeatures = 10000)
LN3l = RunPCA(LN3l, features = LN3l@assays$RNA@var.features)
LN3l =  FindNeighbors(LN3l, reduction = "pca", dims = 1:50)
LN3l = FindClusters(LN3l, resolution = 0.6)

sweep.data = paramSweep_v3(LN3l, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN3l@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN3l$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN3l=doubletFinder_v3(LN3l, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN3l@meta.data[,8]
LN3l = SubsetData(LN3l, cells = colnames(LN3l)[which(doubletsID == "Singlet")])

LN3l = NormalizeData(LN3l)
LN3l = ScaleData(LN3l)
LN3l = FindVariableFeatures(LN3l, selection.method = "vst", nfeatures = 10000)
LN3l = RunPCA(LN3l, features = LN3l@assays$RNA@var.features)
LN3l =  FindNeighbors(LN3l, reduction = "pca", dims = 1:50)
LN3l = FindClusters(LN3l, resolution = 0.6)

#LN3r
LN3r = CreateSeuratObject(counts = LN3r.data, min.cells = 3, min.features = 200, project = "LN3r" )
LN3r[["percent.mt"]] = PercentageFeatureSet(object = LN3r, pattern = "^MT-")
VlnPlot(LN3r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN3r = subset(LN3r, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN3r = NormalizeData(LN3r)
LN3r = ScaleData(LN3r)
LN3r = FindVariableFeatures(LN3r, selection.method = "vst", nfeatures = 10000)
LN3r = RunPCA(LN3r, features = LN3r@assays$RNA@var.features)
LN3r =  FindNeighbors(LN3r, reduction = "pca", dims = 1:50)
LN3r = FindClusters(LN3r, resolution = 0.6)

sweep.data = paramSweep_v3(LN3r, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN3r@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN3r$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN3r=doubletFinder_v3(LN3r, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN3r@meta.data[,8]
LN3r = SubsetData(LN3r, cells = colnames(LN3r)[which(doubletsID == "Singlet")])

LN3r = NormalizeData(LN3r)
LN3r = ScaleData(LN3r)
LN3r = FindVariableFeatures(LN3r, selection.method = "vst", nfeatures = 10000)
LN3r = RunPCA(LN3r, features = LN3r@assays$RNA@var.features)
LN3r =  FindNeighbors(LN3r, reduction = "pca", dims = 1:50)
LN3r = FindClusters(LN3r, resolution = 0.6)



#SC4
SC4 = CreateSeuratObject(counts = SC4.data, min.cells = 3, min.features = 200, project = "SC4" )
SC4[["percent.mt"]] = PercentageFeatureSet(object = SC4, pattern = "^MT-")
VlnPlot(SC4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SC4 = subset(SC4, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
SC4 = NormalizeData(SC4)
SC4 = ScaleData(SC4)
SC4 = FindVariableFeatures(SC4, selection.method = "vst", nfeatures = 10000)
SC4 = RunPCA(SC4, features = SC4@assays$RNA@var.features)
SC4 =  FindNeighbors(SC4, reduction = "pca", dims = 1:50)
SC4 = FindClusters(SC4, resolution = 0.6)

sweep.data = paramSweep_v3(SC4, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(SC4@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(SC4$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
SC4=doubletFinder_v3(SC4, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                     nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = SC4@meta.data[,8]
SC4 = SubsetData(SC4, cells = colnames(SC4)[which(doubletsID == "Singlet")])
SC4 = NormalizeData(SC4)
SC4 = ScaleData(SC4)
SC4 = FindVariableFeatures(SC4, selection.method = "vst", nfeatures = 10000)
SC4 = RunPCA(SC4, features = SC4@assays$RNA@var.features)
SC4 =  FindNeighbors(SC4, reduction = "pca", dims = 1:50)
SC4 = FindClusters(SC4, resolution = 0.6)



#T5
T5 = CreateSeuratObject(counts = T5.data, min.cells = 3, min.features = 200, project = "T5" )
T5[["percent.mt"]] = PercentageFeatureSet(object = T5, pattern = "^MT-")
VlnPlot(T5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T5 = subset(T5, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
T5 = NormalizeData(T5)
T5 = ScaleData(T5)
T5 = FindVariableFeatures(T5, selection.method = "vst", nfeatures = 10000)
T5 = RunPCA(T5, features = T5@assays$RNA@var.features)
T5 =  FindNeighbors(T5, reduction = "pca", dims = 1:50)
T5 = FindClusters(T5, resolution = 0.6)

sweep.data = paramSweep_v3(T5, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T5@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T5$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T5=doubletFinder_v3(T5, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T5@meta.data[,8]
T5 = SubsetData(T5, cells = colnames(T5)[which(doubletsID == "Singlet")])

T5 = NormalizeData(T5)
T5 = ScaleData(T5)
T5 = FindVariableFeatures(T5, selection.method = "vst", nfeatures = 10000)
T5 = RunPCA(T5, features = T5@assays$RNA@var.features)
T5 =  FindNeighbors(T5, reduction = "pca", dims = 1:50)
T5 = FindClusters(T5, resolution = 0.6)





#P5
P5 = CreateSeuratObject(counts = P5.data, min.cells = 3, min.features = 200, project = "P5" )
P5[["percent.mt"]] = PercentageFeatureSet(object = P5, pattern = "^MT-")
VlnPlot(P5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P5 = subset(P5, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
P5 = NormalizeData(P5)
P5 = ScaleData(P5)
P5 = FindVariableFeatures(P5, selection.method = "vst", nfeatures = 10000)
P5 = RunPCA(P5, features = P5@assays$RNA@var.features)
P5 =  FindNeighbors(P5, reduction = "pca", dims = 1:50)
P5 = FindClusters(P5, resolution = 0.6)

sweep.data = paramSweep_v3(P5, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(P5@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(P5$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
P5=doubletFinder_v3(P5, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = P5@meta.data[,8]
P5 = SubsetData(P5, cells = colnames(P5)[which(doubletsID == "Singlet")])


P5 = NormalizeData(P5)
P5 = ScaleData(P5)
P5 = FindVariableFeatures(P5, selection.method = "vst", nfeatures = 10000)
P5 = RunPCA(P5, features = P5@assays$RNA@var.features)
P5 =  FindNeighbors(P5, reduction = "pca", dims = 1:50)
P5 = FindClusters(P5, resolution = 0.6)


#LN5r
LN5r = CreateSeuratObject(counts = LN5r.data, min.cells = 3, min.features = 200, project = "LN5r" )
LN5r[["percent.mt"]] = PercentageFeatureSet(object = LN5r, pattern = "^MT-")
VlnPlot(LN5r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN5r = subset(LN5r, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN5r = NormalizeData(LN5r)
LN5r = ScaleData(LN5r)
LN5r = FindVariableFeatures(LN5r, selection.method = "vst", nfeatures = 10000)
LN5r = RunPCA(LN5r, features = LN5r@assays$RNA@var.features)
LN5r =  FindNeighbors(LN5r, reduction = "pca", dims = 1:50)
LN5r = FindClusters(LN5r, resolution = 0.6)

sweep.data = paramSweep_v3(LN5r, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN5r@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN5r$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN5r=doubletFinder_v3(LN5r, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN5r@meta.data[,8]
LN5r = SubsetData(LN5r, cells = colnames(LN5r)[which(doubletsID == "Singlet")])

LN5r = NormalizeData(LN5r)
LN5r = ScaleData(LN5r)
LN5r = FindVariableFeatures(LN5r, selection.method = "vst", nfeatures = 10000)
LN5r = RunPCA(LN5r, features = LN5r@assays$RNA@var.features)
LN5r =  FindNeighbors(LN5r, reduction = "pca", dims = 1:50)
LN5r = FindClusters(LN5r, resolution = 0.6)



#LN6r
LN6r = CreateSeuratObject(counts = LN6r.data, min.cells = 3, min.features = 200, project = "LN6r" )
LN6r[["percent.mt"]] = PercentageFeatureSet(object = LN6r, pattern = "^MT-")
VLNPlot(LN6r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN6r = subset(LN6r, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN6r = NormalizeData(LN6r)
LN6r = ScaleData(LN6r)
LN6r = FindVariableFeatures(LN6r, selection.method = "vst", nfeatures = 10000)
LN6r = RunPCA(LN6r, features = LN6r@assays$RNA@var.features)
LN6r =  FindNeighbors(LN6r, reduction = "pca", dims = 1:50)
LN6r = FindClusters(LN6r, resolution = 0.6)

sweep.data = paramSweep_v3(LN6r, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN6r@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN6r$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN6r=doubletFinder_v3(LN6r, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN6r@meta.data[,8]
LN6r = SubsetData(LN6r, cells = colnames(LN6r)[which(doubletsID == "Singlet")])

LN6r = NormalizeData(LN6r)
LN6r = ScaleData(LN6r)
LN6r = FindVariableFeatures(LN6r, selection.method = "vst", nfeatures = 10000)
LN6r = RunPCA(LN6r, features = LN6r@assays$RNA@var.features)
LN6r =  FindNeighbors(LN6r, reduction = "pca", dims = 1:50)
LN6r = FindClusters(LN6r, resolution = 0.6)



#LN7r
LN7r = CreateSeuratObject(counts = LN7r.data, min.cells = 3, min.features = 200, project = "LN7r" )
LN7r[["percent.mt"]] = PercentageFeatureSet(object = LN7r, pattern = "^MT-")
VlnPlot(LN7r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN7r = subset(LN7r, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN7r = NormalizeData(LN7r)
LN7r = ScaleData(LN7r)
LN7r = FindVariableFeatures(LN7r, selection.method = "vst", nfeatures = 10000)
LN7r = RunPCA(LN7r, features = LN7r@assays$RNA@var.features)
LN7r =  FindNeighbors(LN7r, reduction = "pca", dims = 1:50)
LN7r = FindClusters(LN7r, resolution = 0.6)

sweep.data = paramSweep_v3(LN7r, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN7r@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN7r$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN7r=doubletFinder_v3(LN7r, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN7r@meta.data[,8]
LN7r = SubsetData(LN7r, cells = colnames(LN7r)[which(doubletsID == "Singlet")])

LN7r = NormalizeData(LN7r)
LN7r = ScaleData(LN7r)
LN7r = FindVariableFeatures(LN7r, selection.method = "vst", nfeatures = 10000)
LN7r = RunPCA(LN7r, features = LN7r@assays$RNA@var.features)
LN7r =  FindNeighbors(LN7r, reduction = "pca", dims = 1:50)
LN7r = FindClusters(LN7r, resolution = 0.6)




#T8
T8 = CreateSeuratObject(counts = T8.data, min.cells = 3, min.features = 200, project = "T8" )
T8[["percent.mt"]] = PercentageFeatureSet(object = T8, pattern = "^MT-")
VlnPlot(T8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T8 = subset(T8, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
T8 = NormalizeData(T8)
T8 = ScaleData(T8)
T8 = FindVariableFeatures(T8, selection.method = "vst", nfeatures = 10000)
T8 = RunPCA(T8, features = T8@assays$RNA@var.features)
T8 =  FindNeighbors(T8, reduction = "pca", dims = 1:50)
T8 = FindClusters(T8, resolution = 0.6)

sweep.data = paramSweep_v3(T8, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T8@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T8$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T8=doubletFinder_v3(T8, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T8@meta.data[,8]
T8 = SubsetData(T8, cells = colnames(T8)[which(doubletsID == "Singlet")])

T8 = NormalizeData(T8)
T8 = ScaleData(T8)
T8 = FindVariableFeatures(T8, selection.method = "vst", nfeatures = 10000)
T8 = RunPCA(T8, features = T8@assays$RNA@var.features)
T8 =  FindNeighbors(T8, reduction = "pca", dims = 1:50)
T8 = FindClusters(T8, resolution = 0.6)


#P8
P8 = CreateSeuratObject(counts = P8.data, min.cells = 3, min.features = 200, project = "P8" )
P8[["percent.mt"]] = PercentageFeatureSet(object = P8, pattern = "^MT-")
VlnPlot(P8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P8 = subset(P8, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
P8 = NormalizeData(P8)
P8 = ScaleData(P8)
P8 = FindVariableFeatures(P8, selection.method = "vst", nfeatures = 10000)
P8 = RunPCA(P8, features = P8@assays$RNA@var.features)
P8 =  FindNeighbors(P8, reduction = "pca", dims = 1:50)
P8 = FindClusters(P8, resolution = 0.6)

sweep.data = paramSweep_v3(P8, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(P8@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(P8$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
P8=doubletFinder_v3(P8, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = P8@meta.data[,8]
P8 = SubsetData(P8, cells = colnames(P8)[which(doubletsID == "Singlet")])

P8 = NormalizeData(P8)
P8 = ScaleData(P8)
P8 = FindVariableFeatures(P8, selection.method = "vst", nfeatures = 10000)
P8 = RunPCA(P8, features = P8@assays$RNA@var.features)
P8 =  FindNeighbors(P8, reduction = "pca", dims = 1:50)
P8 = FindClusters(P8, resolution = 0.6)



#T9
T9 = CreateSeuratObject(counts = T9.data, min.cells = 3, min.features = 200, project = "T9" )
T9[["percent.mt"]] = PercentageFeatureSet(object = T9, pattern = "^MT-")
VlnPlot(T9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T9 = subset(T9, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
T9 = NormalizeData(T9)
T9 = ScaleData(T9)
T9 = FindVariableFeatures(T9, selection.method = "vst", nfeatures = 10000)
T9 = RunPCA(T9, features = T9@assays$RNA@var.features)
T9 =  FindNeighbors(T9, reduction = "pca", dims = 1:50)
T9 = FindClusters(T9, resolution = 0.6)

sweep.data = paramSweep_v3(T9, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T9@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T9$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T9=doubletFinder_v3(T9, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T9@meta.data[,8]
T9 = SubsetData(T9, cells = colnames(T9)[which(doubletsID == "Singlet")])

T9 = NormalizeData(T9)
T9 = ScaleData(T9)
T9 = FindVariableFeatures(T9, selection.method = "vst", nfeatures = 10000)
T9 = RunPCA(T9, features = T9@assays$RNA@var.features)
T9 =  FindNeighbors(T9, reduction = "pca", dims = 1:50)
T9 = FindClusters(T9, resolution = 0.6)

#P9
P9 = CreateSeuratObject(counts = P9.data, min.cells = 3, min.features = 200, project = "P9" )
P9[["percent.mt"]] = PercentageFeatureSet(object = P9, pattern = "^MT-")
VlnPlot(P9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P9 = subset(P9, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
P9 = NormalizeData(P9)
P9 = ScaleData(P9)
P9 = FindVariableFeatures(P9, selection.method = "vst", nfeatures = 10000)
P9 = RunPCA(P9, features = P9@assays$RNA@var.features)
P9 =  FindNeighbors(P9, reduction = "pca", dims = 1:50)
P9 = FindClusters(P9, resolution = 0.6)

sweep.data = paramSweep_v3(P9, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(P9@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(P9$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
P9=doubletFinder_v3(P9, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = P9@meta.data[,8]
P9 = SubsetData(P9, cells = colnames(P9)[which(doubletsID == "Singlet")])

P9 = NormalizeData(P9)
P9 = ScaleData(P9)
P9 = FindVariableFeatures(P9, selection.method = "vst", nfeatures = 10000)
P9 = RunPCA(P9, features = P9@assays$RNA@var.features)
P9 =  FindNeighbors(P9, reduction = "pca", dims = 1:50)
P9 = FindClusters(P9, resolution = 0.6)




#T10
T10 = CreateSeuratObject(counts = T10.data, min.cells = 3, min.features = 200, project = "T10" )
T10[["percent.mt"]] = PercentageFeatureSet(object = T10, pattern = "^MT-")
VlnPlot(T10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
T10 = subset(T10, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
T10 = NormalizeData(T10)
T10 = ScaleData(T10)
T10 = FindVariableFeatures(T10, selection.method = "vst", nfeatures = 2000)
T10 = RunPCA(T10, features = T10@assays$RNA@var.features)
T10 =  FindNeighbors(T10, reduction = "pca", dims = 1:50)
T10 = FindClusters(T10, resolution = 0.6)
T10 = SubsetData(T10, subset.name = "IGKC", high.threshold = 1)
sweep.data = paramSweep_v3(T10, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T10@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.2*length(T10$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T10=doubletFinder_v3(T10, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                     nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T10@meta.data[,8]
T10 = SubsetData(T10, cells = colnames(T10)[which(doubletsID == "Singlet")])

T10 = NormalizeData(T10)
T10 = ScaleData(T10)
T10 = FindVariableFeatures(T10, selection.method = "vst", nfeatures = 10000)
T10 = RunPCA(T10, features = T10@assays$RNA@var.features)
T10 =  FindNeighbors(T10, reduction = "pca", dims = 1:50)
T10 = FindClusters(T10, resolution = 0.6)



#LN10r
LN10r = CreateSeuratObject(counts = LN10r.data, min.cells = 3, min.features = 200, project = "LN10r" )
LN10r[["percent.mt"]] = PercentageFeatureSet(object = LN10r, pattern = "^MT-")
VlnPlot(LN10r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN10r = subset(LN10r, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN10r = NormalizeData(LN10r)
LN10r = ScaleData(LN10r)
LN10r = FindVariableFeatures(LN10r, selection.method = "vst", nfeatures = 10000)
LN10r = RunPCA(LN10r, features = LN10r@assays$RNA@var.features)
LN10r =  FindNeighbors(LN10r, reduction = "pca", dims = 1:50)
LN10r = FindClusters(LN10r, resolution = 0.6)

sweep.data = paramSweep_v3(LN10r, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN10r@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN10r$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN10r=doubletFinder_v3(LN10r, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                       nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN10r@meta.data[,8]
LN10r = SubsetData(LN10r, cells = colnames(LN10r)[which(doubletsID == "Singlet")])

LN10r = NormalizeData(LN10r)
LN10r = ScaleData(LN10r)
LN10r = FindVariableFeatures(LN10r, selection.method = "vst", nfeatures = 10000)
LN10r = RunPCA(LN10r, features = LN10r@assays$RNA@var.features)
LN10r =  FindNeighbors(LN10r, reduction = "pca", dims = 1:50)
LN10r = FindClusters(LN10r, resolution = 0.6)



#SC11
SC11 = CreateSeuratObject(counts = SC11.data, min.cells = 3, min.features = 200, project = "SC11" )
SC11[["percent.mt"]] = PercentageFeatureSet(object = SC11, pattern = "^MT-")
VlnPlot(SC11, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SC11 = subset(SC11, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
SC11 = NormalizeData(SC11)
SC11 = ScaleData(SC11)
SC11 = FindVariableFeatures(SC11, selection.method = "vst", nfeatures = 10000)
SC11 = RunPCA(SC11, features = SC11@assays$RNA@var.features)
SC11 =  FindNeighbors(SC11, reduction = "pca", dims = 1:50)
SC11 = FindClusters(SC11, resolution = 0.6)

sweep.data = paramSweep_v3(SC11, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(SC11@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(SC11$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
SC11=doubletFinder_v3(SC11, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                      nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = SC11@meta.data[,8]
SC11 = SubsetData(SC11, cells = colnames(SC11)[which(doubletsID == "Singlet")])

SC11 = NormalizeData(SC11)
SC11 = ScaleData(SC11)
SC11 = FindVariableFeatures(SC11, selection.method = "vst", nfeatures = 10000)
SC11 = RunPCA(SC11, features = SC11@assays$RNA@var.features)
SC11 =  FindNeighbors(SC11, reduction = "pca", dims = 1:50)
SC11 = FindClusters(SC11, resolution = 0.6)


#LN11r
LN11r = CreateSeuratObject(counts = LN11r.data, min.cells = 3, min.features = 200, project = "LN11r" )
LN11r[["percent.mt"]] = PercentageFeatureSet(object = LN11r, pattern = "^MT-")
VlnPlot(LN11r, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LN11r = subset(LN11r, subset =  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10 )
LN11r = NormalizeData(LN11r)
LN11r = ScaleData(LN11r)
LN11r = FindVariableFeatures(LN11r, selection.method = "vst", nfeatures = 10000)
LN11r = RunPCA(LN11r, features = LN11r@assays$RNA@var.features)
LN11r =  FindNeighbors(LN11r, reduction = "pca", dims = 1:50)
LN11r = FindClusters(LN11r, resolution = 0.6)

sweep.data = paramSweep_v3(LN11r, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(LN11r@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(LN11r$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
LN11r=doubletFinder_v3(LN11r, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                       nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = LN11r@meta.data[,8]
LN11r = SubsetData(LN11r, cells = colnames(LN11r)[which(doubletsID == "Singlet")])

LN11r = NormalizeData(LN11r)
LN11r = ScaleData(LN11r)
LN11r = FindVariableFeatures(LN11r, selection.method = "vst", nfeatures = 10000)
LN11r = RunPCA(LN11r, features = LN11r@assays$RNA@var.features)
LN11r =  FindNeighbors(LN11r, reduction = "pca", dims = 1:50)
LN11r = FindClusters(LN11r, resolution = 0.6)





T1 = RenameCells(T1, add.cell.id = "T1")
P1 = RenameCells(P1, add.cell.id = "P1")
T2 = RenameCells(T2, add.cell.id = "T2")
P2 = RenameCells(P2, add.cell.id = "P2")
LN2l = RenameCells(LN2l, add.cell.id = "LN2l")
T3 = RenameCells(T3, add.cell.id = "T3")
P3 = RenameCells(P3, add.cell.id = "P3")
LN3l = RenameCells(LN3l, add.cell.id = "LN3l")
LN3r = RenameCells(LN3r, add.cell.id = "LN3r")
SC4 = RenameCells(SC4, add.cell.id = "SC4")
T5 = RenameCells(T5, add.cell.id = "T5")
P5 = RenameCells(P5, add.cell.id = "P5")
LN5r = RenameCells(LN5r, add.cell.id = "LN5r")
LN6r = RenameCells(LN6r, add.cell.id = "LN6r")
LN7r = RenameCells(LN7r, add.cell.id = "LN7r")
T8 = RenameCells(T8, add.cell.id = "T8")
P8 = RenameCells(P8, add.cell.id = "P8")
T9 = RenameCells(T9, add.cell.id = "T9")
P9 = RenameCells(P9, add.cell.id = "P9")
T10 = RenameCells(T10, add.cell.id = "T10")
LN10r = RenameCells(LN10r, add.cell.id = "LN10r")
SC11 = RenameCells(SC11, add.cell.id = "SC11")
LN11r = RenameCells(LN11r, add.cell.id = "LN11r")

THCA = merge(T1, c(P1, T2, P2, LN2l, T3, P3, LN3l, LN3r, SC4, 
                   T5, P5, LN5r, LN6r, LN7r, T8, P8, T9,
                   P9, T10, LN10r, SC11, LN11r), project = "THCA")

THCA = NormalizeData(THCA)
THCA = FindVariableFeatures(THCA, selection.method = "vst", nfeatures = 5000)

THCA = ScaleData(THCA)

THCA = RunPCA(THCA, npcs = 50)
THCA = RunUMAP(THCA, reduction = "pca", dims = 1:20)

THCA = FindNeighbors(THCA, reduction = "pca", dims = 1:20)
THCA = FindClusters(THCA, resolution = 1)

#UMAP visualization
DimPlot(THCA, reduction = "umap",  pt.size = 0.6, label =TRUE, repel =  TRUE, label.size = 5,
        cols = c(pal_npg()(10), pal_aaas()(7), pal_jama()(7), pal_lancet()(7))) + NoLegend() 








