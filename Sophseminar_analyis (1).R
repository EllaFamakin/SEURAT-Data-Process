#FIRST THING TO DO WAS INSTALLL THE SEURAT TOOLKIT
install.packages("Seurat")
library(dplyr)
library(Seurat)

library(patchwork)

#Load the dataset

temporal_data <- load("~/Downloads/.RData")

#Quality control
#filtering off cells with less than 200 genes for Quality control and cells that have below the 99.7th percentile for nCount and NFeature


# get quartiles for nCount
quantile(integratedSeurat@meta.data[["nCount_RNA"]], prob=c(.25,.5,.75, .997))
# get quartiles for nFeature
quantile(integratedSeurat@meta.data[["nFeature_RNA"]], prob=c(.25,.5,.75, .997))

#Filtering out cells with less than 200 genes, and also cells that have more than the 99.7th percentile nFeature and nCount

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
filtered_temporal <-  subset(integratedSeurat, subset = nFeature_RNA > 200 & nCount_RNA < 44445.75  & nFeature_RNA < 8347.424  )


#Data normalization

temporal_data <- NormalizeData(temporal_data)


#Data scaling
#all.genes <- rownames(temporal_data)
temporal_data <- ScaleData(temporal_data, features = all.genes)


#PCA/ Feature reduction: principal component reduction 
temporal_data <- RunPCA(temporal_data, features = VariableFeatures(object = temporal_data))






#clustering
#This step takes all similar components and groups them together.


temporal_data <- FindNeighbors(temporal_data, dims = 1:10)
temporal_data <- FindClusters(temporal_data, resolution = 0.5)

#then we genarte umaps

temporal_data <- RunUMAP(temporal_data, dims = 1:10)



#Cluster diagram of cell types

temporal_data <- RunUMAP(temporal_data, dims = 1:10)

DimPlot(temporal_data, reduction = "umap")


#umaps of the cells: Here, I subset each cell type based on their numeric ID and derived images that represent their disease and control status(shows how much disease and control presence there is in a cel type)

cell0 <- subset(filtered_temporal, subset = seurat_clusters == 0 )

cell0 <- FindNeighbors(cell0, dims = 1:10)
cell0 <- FindClusters(cell0, resolution = 0.5)

cell0 <- RunUMAP(cell0, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(cell0, reduction = "umap", group.by ='Disease')

cell1 <- subset(filtered_temporal, subset = seurat_clusters == 1 )
cell1 <- FindNeighbors(cell1, dims = 1:10)
cell1 <- FindClusters(cell1, resolution = 0.5)
DimPlot(cell1, reduction = "umap", group.by ='Disease')


cell2 <- subset(filtered_temporal, subset = seurat_clusters == 2 )
cell2 <- FindNeighbors(cell2, dims = 1:10)
cell2 <- FindClusters(cell2, resolution = 0.5)
DimPlot(cell2, reduction = "umap", group.by ='Disease')


cell3 <- subset(filtered_temporal, subset = seurat_clusters == 3 )
cell3 <- FindNeighbors(cell3, dims = 1:10)
cell3 <- FindClusters(cell3, resolution = 0.5)
DimPlot(cell3, reduction = "umap", group.by ='Disease')


cell4 <- subset(filtered_temporal, subset = seurat_clusters == 4 )
cell4 <- FindNeighbors(cell4, dims = 1:10)
cell4 <- FindClusters(cell4, resolution = 0.5)
DimPlot(cell4, reduction = "umap", group.by ='Disease')


cell5 <- subset(filtered_temporal, subset = seurat_clusters == 5 )
cell5 <- FindNeighbors(cell5, dims = 1:10)
cell5 <- FindClusters(cell5, resolution = 0.5)
DimPlot(cell5, reduction = "umap", group.by ='Disease')


cell6 <- subset(filtered_temporal, subset = seurat_clusters == 6 )
cell6 <- FindNeighbors(cell6, dims = 1:10)
cell6 <- FindClusters(cell6, resolution = 0.5)
DimPlot(cell6, reduction = "umap", group.by ='Disease')


cell7 <- subset(filtered_temporal, subset = seurat_clusters == 7 )
cell7 <- FindNeighbors(cell7, dims = 1:10)
cell7 <- FindClusters(cell7, resolution = 0.5)
DimPlot(cell7, reduction = "umap", group.by ='Disease')


cell8 <- subset(filtered_temporal, subset = seurat_clusters == 8 )
cell8 <- FindNeighbors(cell8, dims = 1:10)
cell8 <- FindClusters(cell8, resolution = 0.5)
DimPlot(cell8, reduction = "umap", group.by ='Disease')


cell9 <- subset(filtered_temporal, subset = seurat_clusters == 9)
cell9 <- FindNeighbors(cell9, dims = 1:10)
cell9 <- FindClusters(cell9, resolution = 0.5)
DimPlot(cell9, reduction = "umap", group.by ='Disease')



