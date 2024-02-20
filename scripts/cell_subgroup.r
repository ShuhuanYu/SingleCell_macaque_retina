
library(Seurat)
library(ggplot2)
library(gplots)
library(patchwork)
library(dplyr)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(ArchR)
library(harmony)

# load data
macaque <- readRDS("./cell_anno_result/macaque_scRNA_harmony.rds")

# do subcluster
celltype <- "Müller"

celltype_cells <- macaque[, Idents(macaque) %in% celltype] 

celltype_cells <- FindVariableFeatures(celltype_cells, selection.method = 'vst', nfeatures = 2000)
celltype_cells <- ScaleData(celltype_cells,features=rownames(celltype_cells))
celltype_cells <- RunPCA(celltype_cells, features = VariableFeatures(object = celltype_cells)) 
ElbowPlot(celltype_cells, ndims=20, reduction="pca")

celltype_cells <- RunHarmony(celltype_cells, group.by.vars="Retina.Samples", plot_convergence=TRUE)

celltype_cells <- FindNeighbors(celltype_cells, dims = 1:15)
celltype_cells <- FindClusters(celltype_cells, resolution = 0.05 )
table(celltype_cells$seurat_clusters) 
  #  0    1    2 
# 4239  546   25
celltype_cells <- RenameIdents(celltype_cells, 
                        '0' = 'Müller1',
                        '1' = 'Müller2',
                        '2' = 'Müller3'
                        # '3' = 'Müller4',
                        # '4' = 'Müller5',
                        # '5' = 'Müller6'
                        )
celltype_cells$celltype <- Idents(celltype_cells)

celltype_cells <- RunUMAP(celltype_cells, dims = 1:15, reduction="harmony")

# Müller_subset <- celltype_cells[, Idents(celltype_cells) %in% c("Müller1","Müller2","Müller3","Müller4")] # remove 4,5
# Müller_subset$celltype <- Idents(Müller_subset)
p1 <- DimPlot(celltype_cells, reduction ="umap", label = TRUE)
p2 <- DimPlot(celltype_cells, reduction="umap", group.by="Retina.Samples")
ggsave("./cell_subsubgroup/Müller/Müller_subsubgroup_dimplot.pdf", p1+p2,width=8,height=4)

# find all markers
diff.wilcox = FindAllMarkers(celltype_cells) # Finds markers (differentially expressed genes) for each of the identity classes in a dataset
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "./cell_subsubgroup/Müller/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "./cell_subsubgroup/Müller/top10_diff_genes_wilcox.csv", row.names = F)

sample.distribution.in.subclusters <- as.matrix(table(celltype_cells$celltype,celltype_cells$Retina.Samples))
write.csv(sample.distribution.in.subclusters,"./cell_subsubgroup/Müller/sample.distribution.in.subclusters.csv",row.names=T)

# head(celltype_cells@assays[["RNA"]]@data)[,1:4]

# enrichment analysis of all markers genes(多个基因集富集分析)
ids=bitr(all.markers$gene,'SYMBOL','ENTREZID','org.Mmu.eg.db') ## 将SYMBOL转成ENTREZID
all.markers=merge(all.markers,ids,by.x='gene',by.y='SYMBOL')

# kegg
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(all.markers$ENTREZID, all.markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
  fun = "enrichKEGG",
  organism = "mcc", pvalueCutoff = 0.05
)
table(xx@compareClusterResult$Cluster)
head(as.data.frame(xx))
p1 <- dotplot(xx)


## GO
xx <- compareCluster(gcSample,
  fun = "enrichGO",
  OrgDb = "org.Mmu.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)
table(xx@compareClusterResult$Cluster)
p2 <- dotplot(xx)

ggsave("./cell_subsubgroup/Müller/Müller_subsubgroup_kegg_go.pdf",p1+p2,width=15,height=10)

# ATAC location
proj <- loadArchRProject(path = "./ArchRProject_Merged/")
cells_barocodes <- getCellNames(proj)[which(proj$predictedGroup=="Müller")]
Müller_proj = subsetArchRProject(
        ArchRProj = proj,
        cells = cells_barocodes,
        outputDirectory = "ArchRSubset",
        dropCells = TRUE,
        logFile = NULL,
        threads = getArchRThreads(),
        force = TRUE
    )
# unconstrained integration
Müller_proj <- addGeneIntegrationMatrix(
    ArchRProj = Müller_proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = celltype_cells,
    addToArrow = TRUE,
    groupRNA = "celltype",
    nameCell = "predictedCell2",
    nameGroup = "predictedGroup2",
    nameScore = "predictedScore2", 
    force = TRUE
)
table(Müller_proj$predictedGroup2)
# Müller1 Müller2 
#    4292      45 
# UMAP with cell type reannotated
p1 <- plotEmbedding(
    Müller_proj, 
    embedding = "UMAP",
    colorBy = "cellColData", 
    name = "predictedGroup2"
)
plotPDF(p1, name = "scATAC-Plot-UMAP-RNA-Integration.pdf", ArchRProj = Müller_proj, addDOC = FALSE, width = 5, height = 5)


# enrichment of each subtype-specific top50 marker genes
diff.roc <- FindAllMarkers(celltype_cells, test.use = "roc", only.pos = TRUE)
top50 = diff.roc %>% dplyr::select(gene, everything()) %>% subset(power > 0.4) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) 
write.csv(diff.roc, "./cell_subsubgroup/Müller/diff_genes_roc.csv",row.names=TRUE,quote=F)
write.csv(top50, "./cell_subsubgroup/Müller/top50_diff_genes_roc.csv", row.names=F,quote=F)
# top50 <- top50[which(grepl("ENSMMUG",top50$gene)=="FALSE"),]
p <- DoHeatmap(celltype_cells, features = top50$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))
ggsave("./cell_subsubgroup/Müller/subtype-specific_top50markers_expression_heatmap.pdf",p)


saveRDS(celltype_cells,"./cell_subsubgroup/Müller/celltype_cells.rds")
# save.image("./cell_subsubgroup/Müller/Müller-part1.Rdata")


####################################### compare expression between young and old ######################################
library(DESeq2)


celltype_cells <- readRDS("./cell_subsubgroup/Müller/celltype_cells.rds")

Müller1 <- celltype_cells[,Idents(celltype_cells)=="Müller1"]
Müller2 <- celltype_cells[,Idents(celltype_cells)=="Müller2"]
Müller3 <- celltype_cells[,Idents(celltype_cells)=="Müller3"]

Idents(Müller1) <- Müller1$Retina.Samples
age.related.degs.Müller1.central <- FindMarkers(Müller1,test.use="wilcox", ident.1="E-OS-C",ident.2="H-OS-C") %>% subset(p_val_adj<0.05 & abs(avg_log2FC)>0.5)
age.related.degs.Müller1.periphery <- FindMarkers(Müller1,test.use="wilcox", ident.1="D-OS-P",ident.2="H-OS-P") %>% subset(p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(age.related.degs.Müller1.central,"./cell_subsubgroup/Müller/age.related.degs.Müller1.central.csv",row.names=T,quote=F)
write.csv(age.related.degs.Müller1.periphery,"./cell_subsubgroup/Müller/age.related.degs.Müller1.periphery.csv",row.names=T,quote=F)

Idents(Müller2) <- Müller2$Retina.Samples
age.related.degs.Müller2.central <- FindMarkers(Müller2,test.use="wilcox", ident.1="E-OS-C",ident.2="H-OS-C") %>% subset(p_val_adj<0.05 & abs(avg_log2FC)>0.5)
age.related.degs.Müller2.periphery <- FindMarkers(Müller2,test.use="wilcox", ident.1="D-OS-P",ident.2="H-OS-P") %>% subset(p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(age.related.degs.Müller2.central,"./cell_subsubgroup/Müller/age.related.degs.Müller2.central.csv",row.names=T,quote=F)
write.csv(age.related.degs.Müller2.periphery,"./cell_subsubgroup/Müller/age.related.degs.Müller2.periphery.csv",row.names=T,quote=F)

Idents(Müller3) <- Müller3$Retina.Samples
age.related.degs.Müller3.central <- FindMarkers(Müller3,test.use="wilcox", ident.1="E-OS-C",ident.2="H-OS-C") %>% subset(p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(age.related.degs.Müller3.central,"./cell_subsubgroup/Müller/age.related.degs.Müller3.central.csv",row.names=T,quote=F)
# age.related.degs.Müller3.periphery <- FindMarkers(Müller3,test.use="wilcox", ident.1="D-OS-P",ident.2="H-OS-P")
# Error in ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2,  : 
#   Cell group 2 has fewer than 3 cells

save.image("./cell_subsubgroup/Müller/Müller-part2.Rdata")

############################################################################################################

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

youngST <- readRDS("./spatial_result/rawData/sample_youngST_bin100_seurat.rds")
youngST <- SCTransform(youngST, assay = "Spatial", verbose = FALSE)

# Dimensionality reduction, clustering, and visualization
youngST <- RunPCA(youngST, assay = "SCT", verbose = FALSE)
youngST <- FindNeighbors(youngST, reduction = "pca", dims = 1:30)
youngST <- FindClusters(youngST, verbose = FALSE,resolution=0.4)
youngST <- RunUMAP(youngST, reduction = "pca", dims = 1:30)

youngST <- subset(youngST,idents=c(3,4,7,8,9))

macaque <- readRDS("./cell_anno_result/macaque_scRNA_harmony.rds")
youngSC <- macaque[,macaque@meta.data$Retina.Samples %in% "E-OS-C"]
youngSC <- SCTransform(youngSC, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
DimPlot(youngSC, group.by = "celltype", label = TRUE)

celltype_cells <- readRDS("./cell_subsubgroup/Müller/celltype_cells.rds")

subclass <- as.data.frame(as.character(youngSC$celltype))
rownames(subclass) <- rownames(youngSC@meta.data)
Müller.young.subclass <- as.data.frame(as.character(celltype_cells$celltype))
rownames(Müller.young.subclass) <- rownames(celltype_cells@meta.data)
subclass[intersect(rownames(subclass),rownames(Müller.young.subclass)),1] <- Müller.young.subclass[intersect(rownames(subclass),rownames(Müller.young.subclass)),1]

youngSC$subclass <- subclass[,1]

anchors <- FindTransferAnchors(reference = youngSC, query = youngST, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = youngSC$subclass, prediction.assay = TRUE,
    weight.reduction = youngST[["pca"]], dims = 1:30)
youngST[["predictions_subclass"]] <- predictions.assay

DefaultAssay(youngST) <- "predictions_subclass"

p <- SpatialFeaturePlot(youngST, features = c("Müller1","Müller2","Müller3"), ncol = 3, crop = FALSE,stroke = NA,pt.size.factor = 4)
ggsave("./spatial_result/Plots/spatial_Müller_subtype_Featureplot.pdf",p,width=15,height=6)
