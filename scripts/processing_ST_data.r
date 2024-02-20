library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratObject)

### old ###
ST <- readRDS("./spatial_result/rawData/sample_ST_bin100_seurat.rds") # 此处加载不同的空转rds对象

sample_label <- "old" # 此处更改样本标签

# data processing
plot1 <- VlnPlot(ST, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(ST, features = "nCount_Spatial",stroke = NA,pt.size.factor = 4) + theme(legend.position = "right")
plot3 <- VlnPlot(ST, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot4 <- SpatialFeaturePlot(ST, features = "nFeature_Spatial", stroke = NA,pt.size.factor = 4) + theme(legend.position = "right")
ggsave(paste0("./spatial_result/Plots/spatial_",sample_label,"_data_processing_nCounts-1.pdf"),plot1)
ggsave(paste0("./spatial_result/Plots/spatial_",sample_label,"_data_processing_nCounts-2.pdf"),plot2)
ggsave(paste0("./spatial_result/Plots/spatial_",sample_label,"_data_processing_nFeatures-1.pdf"),plot3)
ggsave(paste0("./spatial_result/Plots/spatial_",sample_label,"_data_processing_nFeatures-2.pdf"),plot4)

ST <- SCTransform(ST, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction, clustering, and visualization
ST <- RunPCA(ST, assay = "SCT", verbose = FALSE)
ST <- FindNeighbors(ST, reduction = "pca", dims = 1:30)
ST <- FindClusters(ST, verbose = FALSE,resolution=0.4)
ST <- RunUMAP(ST, reduction = "pca", dims = 1:30)

subset.ST <- subset(ST, idents = c(4,5,6))

p1 <- DimPlot(subset.ST, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(subset.ST, label = TRUE, label.size = 2,stroke = NA,pt.size.factor = 5)
ggsave(paste0("./spatial_result/Plots/spatial_",sample_label,"_UMAP_and_STDimplot.pdf"),p1+p2,width=10,height=6)

# Identification of Spatially Variable Features
# ST <- FindSpatiallyVariableFeatures(ST, assay = "SCT", features = VariableFeatures(ST)[1:1000],
    # selection.method = "moransi")

# Integration with single-cell data
macaque <- readRDS("./cell_anno_result/macaque_scRNA_harmony.rds") # 这里更改为新的转录组rds对象

SC <- macaque[,macaque@meta.data$Retina.Samples %in% "E-OS-C"] # 选择与空转样本condition一致的单细胞转录组样本
SC <- SCTransform(SC, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
# DimPlot(SC, group.by = "celltype", label = TRUE)

anchors <- FindTransferAnchors(reference = SC, query = subset.ST, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = SC$celltype, prediction.assay = TRUE,
    weight.reduction = subset.ST[["pca"]], dims = 1:30)
subset.ST[["predictions"]] <- predictions.assay

DefaultAssay(subset.ST) <- "predictions"

p <- SpatialFeaturePlot(subset.ST, features = unique(SC$celltype), ncol = 4, crop = FALSE,stroke = NA,pt.size.factor = 5)
ggsave(paste0("./spatial_result/Plots/spatial_",sample_label,"_labeled_Featureplot.pdf"),p,width=20,height=20)

saveRDS(subset.ST,paste0("./spatial_result/result/",sample_label,"_subset.ST.rds"))
# save.image("./spatial_result/result/old_subset.ST-part1.Rdata")
