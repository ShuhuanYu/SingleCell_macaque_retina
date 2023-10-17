## Author: Yu Shuhuan
## Date: 2023-10-17
## Brief description: preprocessing macaque single cell data

```{r}
setwd("/opt/yushuhuan/macaque")

library(Seurat)
options(Seurat.object.assay.version = "v5")
library(EnsDb.Mmulatta.v109)
library(Signac)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(DoubletFinder)
library(patchwork)
library(hdf5r)
library(harmony)
library(cowplot)
library(SeuratWrappers)

### pre-processing data  

for(i in 1:6){
  samples=c("E-OS-C","A-OS-C","H-OS-C","D-OS-P","A-OS-P","H-OS-P")
  sample <- samples[i]
  print(paste0("Now processing: ",sample))

  h5_file <- paste0("./cellrangerARC_output/",sample,"/outs/filtered_feature_bc_matrix.h5")
  inputdata.10x <- Read10X_h5(h5_file)
  rna_counts <- inputdata.10x$`Gene Expression`
  macaque <- CreateSeuratObject(counts = rna_counts)
  macaque[["percent.mt"]] <- PercentageFeatureSet(macaque, pattern = "^MT-")
  raw_cell_number <- dim(macaque)[2]

  #quality control
  macaque <- subset(macaque, subset = nFeature_RNA > 200 & nFeature_RNA < 25000 & percent.mt < 5)
  cell_numbers_after_qc <- dim(macaque)[2]

  macaque <- NormalizeData(macaque)
  macaque <- FindVariableFeatures(macaque, selection.method = "vst", nfeatures = 2000)
  macaque <- ScaleData(macaque)
  macaque <- RunPCA(macaque)
  macaque <- RunUMAP(macaque, dims = 1:20)

  # cluster cells
  macaque <- FindNeighbors(macaque, dims = 1:20)
  macaque <- FindClusters(macaque,resolution=0.2)
  
  # pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep_v3(macaque, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # remofve doublets by DoubletFinder
  macaque$celltype <- Idents(macaque)
  DoubletRate=ncol(macaque)*8*1e-6
  homotypic.prop=modelHomotypic(macaque@meta.data$celltype)
  nExp_poi=round(DoubletRate*length(macaque$celltype))
  nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
  macaque=doubletFinder_v3(macaque, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  DF.clarrsications=names(macaque@meta.data)[grep("DF.classifications",names(macaque@meta.data))]
  
  # Plot results
  p1 <- DimPlot(macaque, reduction = "umap", group.by = DF.clarrsications) + ggtitle("Doublets")
  ggsave(paste0("./cell_anno_result/",sample,"_DoubletFinder.pdf"),p1,width = 10,height = 10)

  # filter doublets
  Idents(macaque ) <- DF.clarrsications
  macaque<-subset(x = macaque, idents="Singlet")   #过滤细胞，只保留经过检测的Singlet
  filtered_cell_number <- dim(macaque)[2]

  print(paste0("raw cell numbers: ",raw_cell_number))
  print(paste0("cell numbers after quality control: ", cell_numbers_after_qc))
  print(paste0("filtered cell numbers: ",filtered_cell_number))

  # save filtered rna_counts 
  filtered_barcodes <- colnames(macaque[["RNA"]]@counts)
  filtered.rna.counts = rna_counts[,filtered_barcodes]
  colnames(filtered.rna.counts) = paste0(sample,".",colnames(filtered.rna.counts))
  write.table(filtered.rna.counts, paste0("./cell_anno_result/",sample,"_filtered.rna.counts.txt"),row.names = T,sep = "\t",quote = F)

  rm(list = ls())
  gc()
}

### Aggregate all samples

E_OS_C.counts <- read.table("./cell_anno_result/E-OS-C_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
A_OS_C.counts <- read.table("./cell_anno_result/A-OS-C_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
H_OS_C.counts <- read.table("./cell_anno_result/H-OS-C_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
D_OS_P.counts <- read.table("./cell_anno_result/D-OS-P_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
A_OS_P.counts <- read.table("./cell_anno_result/A-OS-P_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
H_OS_P.counts <- read.table("./cell_anno_result/H-OS-P_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)

all.samples.counts <- cbind(E_OS_C.counts, A_OS_C.counts, H_OS_C.counts, D_OS_P.counts, A_OS_P.counts, H_OS_P.counts)
rm("E_OS_C.counts", "A_OS_C.counts", "H_OS_C.counts", "D_OS_P.counts", "A_OS_P.counts", "H_OS_P.counts")

macaque <- CreateSeuratObject(counts = all.samples.counts, project = "macaque", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = macaque@var.genes, npcs = 20, verbose = FALSE)

macaque@meta.data$Retina.Samples <- c(rep("E-OS-C", 12270), rep("A-OS-C", 9555), rep("H-OS-C", 12092), rep("D-OS-P", 9301), rep("A-OS-P", 8857), rep("H-OS-P", 11748))

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = macaque, reduction = "pca", pt.size = .1, group.by = "Retina.Samples")
p2 <- VlnPlot(object = macaque, features = "PC_1", group.by = "Retina.Samples", pt.size = .1)
ggsave("./cell_anno_result/plot_uncorrected_PCs_samples.pdf", plot_grid(p1,p2), width = 12, height = 6)

# run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
macaque <- macaque %>% 
    RunHarmony("Retina.Samples", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(macaque, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = macaque, reduction = "harmony", pt.size = .1, group.by = "Retina.Samples")
p2 <- VlnPlot(object = macaque, features = "harmony_1", group.by = "Retina.Samples", pt.size = .1)
ggsave("./cell_anno_result/plot_corrected_PCs_samples.pdf", plot_grid(p1,p2), width = 12, height = 6)

# UMAP analysis
macaque <- RunUMAP(macaque, reduction = "harmony", dims = 1:20) 
macaque <- FindNeighbors(macaque, reduction = "harmony", dims = 1:20)
macaque <- FindClusters(macaque, resolution = 0.2)

p1 <- DimPlot(macaque, reduction = "umap", group.by = "Retina.Samples", pt.size = .1, split.by = 'Retina.Samples')
ggsave("./cell_anno_result/plot_unlabel_Harmony_dimplot.pdf", p1, width = 30,height = 7)

p2 <- DimPlot(macaque, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("./cell_anno_result/plot_unlabel_Harmony_umap.pdf", p2, width = 7,height = 7)

# annotation clusters
markers=c(Rod=c("PDE6A","PDE6B"),
          Cone=c("ARR3"), #"KCNB2"
          Horizontal_cell=c("THSD7B","RET"),
          Amacrine_cell=c("SLC17A8","GRIP1","GAD1","GAD2","PTPRK","SLC6A9","CALB2"),
          RGC=c("SLC17A6","GALNTL6","IL1RAPL2"),
          Microglia=c("C1QB","C1QC","CX3CR1","APBB1IP"),
          Müller_cell=c("RLBP1","HKDC1","WIF1"),
          Bipolar_cell=c("GRIK1","VSX2","TRPM1","CCDC136")
          )
cd_genes=unique(markers)

p1 <- DotPlot(macaque,features = cd_genes)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("./cell_anno_result/number_anno_dotplot.pdf", p1, width = 15, height = 7)

macaque <- RenameIdents(macaque, 
                        '0' = 'Rod',
                        '1' = 'Rod',
                        '2' = 'Müller',
                        '3' = 'Cone',
                        '4' = 'GABA-amacrine',
                        '5' = 'ON-cone and rod bipolar',
                        '6' = 'OFF-cone and rod bipolar',
                        '7' = 'OFF-cone and rod bipolar',
                        '8' = 'ON-cone and rod bipolar',
                        '9' = 'Gly-amacrine',
                        '10' = 'RGC',
                        '11' = 'OFF-cone and rod bipolar',
                        '12' = 'Horizontal',
                        '13' = 'ON-cone and rod bipolar',
                        '14' = 'Gly-amacrine',
                        '15' = 'ON-cone and rod bipolar',
                        '16' = 'Amacrine',
                        '17' = 'OFF-cone and rod bipolar',
                        '18' = 'GABA-amacrine',
                        '19' = 'Bipolar',
                        '20' = 'RGC',
                        '21' = 'Microglia',
                        '22' = 'Cone'
                        )
# find all markers od cluster "XXX"
# XXX.markers <- FindMarkers(macaque, ident.1 = "19", min.pct = 0.25)
# XXX.markers[1:10,]

macaque$celltype <- Idents(macaque)
p1 <- DimPlot(macaque, reduction = "umap", group.by = "celltype", label = FALSE, label.size = 2.5, repel = TRUE, split.by = "Retina.Samples")
ggsave("./cell_anno_result/all_samples_anno_label_dimplot.pdf", p1, width = 30,height = 7)

levels(macaque)=c("Bipolar","ON-cone and rod bipolar","OFF-cone and rod bipolar","Müller",
                "Microglia","RGC","Gly-amacrine","GABA-amacrine","Amacrine","Horizontal",
                "Cone","Rod")
p2 <- DotPlot(macaque,features = cd_genes)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("./cell_anno_result/celltype_anno_dotplot.pdf", p2, width = 15, height = 7)

p3 <- DimPlot(macaque, reduction = "umap", group.by = "celltype", pt.size = .1)
ggsave("./cell_anno_result/plot_label_Harmony_dimplot.pdf", p3, width = 9,height = 7)

saveRDS(macaque,"./cell_anno_result/macaque_scRNA_harmony.rds")

### calculate cells percentage
data <- as.data.frame(table(macaque$Retina.Samples,macaque$celltype))
colnames(data) <- c("Retina.Samples","CellType","Freq")
df <- data %>% 
  group_by(Retina.Samples) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total)

df$CellType  <- factor(df$CellType,levels = unique(df$CellType))

write.csv(df,file = "./cell_anno_result/allsamples_cell_percent.csv",row.names = F,quote = F)

# save RNA counts and metadata
rna_count <- macaque[["RNA"]]@counts
metadata <- macaque@meta.data

saveRDS(rna_count, "./cell_anno_result/scRNA_count.rds")
saveRDS(metadata, "./cell_anno_result/scRNA_meta.data.rds")

```
