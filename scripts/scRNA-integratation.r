
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



### Integration 10 retina samples

E_OS_C.counts <- read.table("./result/1.remove_doublets-scRNA/E-OS-C/E-OS-C_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
A_OS_C.counts <- read.table("./result/1.remove_doublets-scRNA/A-OS-C/A-OS-C_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
H_OS_C.counts <- read.table("./result/1.remove_doublets-scRNA/H-OS-C/H-OS-C_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
D_OS_P.counts <- read.table("./result/1.remove_doublets-scRNA/D-OS-P/D-OS-P_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
A_OS_P.counts <- read.table("./result/1.remove_doublets-scRNA/A-OS-P/A-OS-P_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
H_OS_P.counts <- read.table("./result/1.remove_doublets-scRNA/H-OS-P/H-OS-P_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
COS_H.counts <- read.table("./result/1.remove_doublets-scRNA/COS-H/COS-H_filtered.rna.counts.txt", header = T, sep = "\t", row.names = 1)
EOS_Z.counts <- read.table("./result/1.remove_doublets-scRNA/EOS-Z/EOS-Z_filtered.rna.counts.txt", header = T, sep = "\t", row.names = 1)
IOS_H.counts <- read.table("./result/1.remove_doublets-scRNA/IOS-H/IOS-H_filtered.rna.counts.txt", header = T, sep = "\t", row.names = 1)
IOS_Z.counts <- read.table("./result/1.remove_doublets-scRNA/IOS-Z/IOS-Z_filtered.rna.counts.txt", header = T, sep = "\t", row.names = 1)

dim(E_OS_C.counts)
dim(A_OS_C.counts)
dim(H_OS_C.counts)
dim(D_OS_P.counts)
dim(A_OS_P.counts)
dim(H_OS_P.counts)
dim(COS_H.counts)
dim(EOS_Z.counts)
dim(IOS_H.counts)
dim(IOS_Z.counts)

all.samples.counts <- list(E_OS_C.counts,A_OS_C.counts,H_OS_C.counts,D_OS_P.counts,A_OS_P.counts,H_OS_P.counts,COS_H.counts,E_OS_C.counts,IOS_H.counts,IOS_Z.counts)
names(all.samples.counts) <- c("E-OS-C","A-OS-C","H-OS-C","D-OS-P","A-OS-P","H-OS-P","COS-H","EOS-Z","IOS-H","IOS-Z")

# 不同样本创建不同Seurat对象，需按照高变基因，取高变基因并集作为最终合并分析
samples <- names(all.samples.counts)
all.hvf <- vector()
macaque <- list()
for(i in 1:10){
  sample.seurat.object = CreateSeuratObject(counts = all.samples.counts[[samples[i]]], project = samples[i], min.cells = 3, min.features = 200) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  macaque = c(macaque, list(sample.seurat.object))
  top2000 = VariableFeatures(sample.seurat.object)
  all.hvf = c(all.hvf, top2000)
}
all.hvf <- unique(all.hvf)  # 7139个高变基因

rm(list = c("E_OS_C.counts","A_OS_C.counts","H_OS_C.counts","D_OS_P.counts","A_OS_P.counts","H_OS_P.counts","COS_H.counts","E_OS_C.counts","IOS_H.counts","IOS_Z.counts"))
# 帮助确定min.cells和min.features
# library(stringr)  
# fivenum(apply(all.samples.counts,1,function(x) sum(x>0) ))  # by row 
# boxplot(apply(all.samples.counts,1,function(x) sum(x>0) ))
# fivenum(apply(all.samples.counts,2,function(x) sum(x>0) ))  # by column
# hist(apply(all.samples.counts,2,function(x) sum(x>0) ))

macaque.combined <- SeuratObject::merge(macaque[[1]],y=c(macaque[[2]],macaque[[3]],macaque[[4]],macaque[[5]],
                        macaque[[6]],macaque[[7]],macaque[[8]],macaque[[9]],macaque[[10]]
                        ),
                        add.cells.ids = samples,
                        project = "macaque.retina.combined"
                        )
table(macaque.combined$orig.ident)
macaque.combined <- ScaleData(macaque.combined, features = all.hvf)
macaque.combined <- RunPCA(macaque.combined, features = all.hvf, npcs = 50, verbose = FALSE)
macaque.combined
# An object of class Seurat 
# 22639 features across 96674 samples within 1 assay 
# Active assay: RNA (22639 features, 0 variable features)
#  1 dimensional reduction calculated: pca

macaque.combined$Samples <- macaque.combined$orig.ident

p1 <- DimPlot(object = macaque.combined, reduction = "pca", pt.size = .1, group.by = "Samples")
p2 <- VlnPlot(object = macaque.combined, features = "PC_1", group.by = "Samples", pt.size = .1)
ggsave("./result/2.samples_integration/scRNA_integration//plot_uncorrected_PCs_samples.pdf", plot_grid(p1,p2), width = 12, height = 6)

# run Harmony
macaque.combined <- macaque.combined %>% 
    RunHarmony("Samples", plot_convergence = TRUE, max_iter = 20)

harmony_embeddings <- Embeddings(macaque.combined, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = macaque.combined, reduction = "harmony", pt.size = .1, group.by = "Samples")
p2 <- VlnPlot(object = macaque.combined, features = "harmony_1", group.by = "Samples", pt.size = .1)
ggsave("./result/2.samples_integration/scRNA_integration//plot_corrected_PCs_samples.pdf", plot_grid(p1,p2), width = 12, height = 6)

# UMAP analysis
ElbowPlot(macaque.combined,ndims = 50)
macaque.combined <- RunUMAP(macaque.combined, reduction = "harmony", dims = 1:25) 
macaque.combined <- FindNeighbors(macaque.combined, reduction = "harmony", dims = 1:25)
macaque.combined <- FindClusters(macaque.combined, resolution = 0.6)
#     0     1     2     3     4     5     6     7     8     9    10    11    12 
# 23476 18899 10818  7634  6494  3756  2312  1982  1809  1773  1592  1509  1383 
#    13    14    15    16    17    18    19    20    21    22    23    24    25 
#  1240  1222   936   921   835   785   784   755   715   683   507   485   395 
#    26    27    28    29    30    31    32    33    34    35    36    37 
#   355   352   350   313   295   287   279   274   210   129   128     2 

p1 <- DimPlot(macaque.combined, reduction = "umap", group.by = "Samples", pt.size = .1, split.by = 'Samples')
ggsave("./result/2.samples_integration/scRNA_integration//plot_unlabel_Harmony_dimplot.pdf", p1, width = 30,height = 7)

p2 <- DimPlot(macaque.combined, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("./result/2.samples_integration/scRNA_integration/plot_unlabel_Harmony_umap.pdf", p2, width = 7,height = 7)

macaque.combined <- subset(macaque.combined, idents = c())
# annotation clusters
markers=c(Rod=c("PDE6A","PDE6B"),
          Cone=c("ARR3"), #"KCNB2"
          Horizontal_cell=c("THSD7B","RET"),
          Amacrine_cell=c("SLC17A8","ERBB4","GRIP1","GAD1","GAD2","PTPRK","SLC6A9","CALB2"),
          RGC=c("SLC17A6","GALNTL6","IL1RAPL2"),
          Microglia=c("C1QB","C1QC","CX3CR1","APBB1IP"),
          Müller_cell=c("RLBP1","HKDC1","WIF1"),
          Bipolar_cell=c("CNTNAP5","GRIK1","VSX2","TRPM1","CCDC136")
          )
cd_genes=unique(markers)

p1 <- DotPlot(macaque.combined,features = cd_genes)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("./result/2.samples_integration/scRNA_integration/number_anno_dotplot.pdf", p1, width = 15, height = 10)
 
p2 <- DimPlot(macaque.combined,reduction = "umap",group.by = "Samples", label = TRUE, pt.size = .1)
ggsave("./result/2.samples_integration/scRNA_integration/sample_anno_dimplot.pdf", p2, width = 15, height = 10)

saveRDS(macaque.combined,"./result/2.samples_integration/scRNA_integration/macaque.combined.tmp.rds")
# t=FindMarkers(macaque,ident.1="24",ident.2="7")

macaque.combined <- readRDS("./result/2.samples_integration/scRNA_integration/macaque.combined.tmp.rds")
macaque.combined <- RenameIdents(macaque.combined, 
                        '0' = 'Rod',
                        '1' = 'Rod',
                        '2' = 'Rod',
                        '3' = 'Rod',
                        '4' = 'Müller glia',
                        '5' = 'Cone',
                        '6' = 'ON-cone and rod bipolar',
                        '7' = 'GABA-amacrine',
                        '8' = 'Gly-amacrine',
                        '9' = 'OFF-cone and rod bipolar',
                        '10' = 'Rod',
                        '11' = 'ON-cone and rod bipolar',
                        '12' = 'Retinal ganglion cell',
                        '13' = 'Horizontal',
                        '14' = 'OFF-cone and rod bipolar',
                        '15' = 'Gly-amacrine',
                        '16' = 'Bipolar',
                        '17' = 'ON-cone and rod bipolar',
                        '18' = 'Müller glia',
                        '19' = 'OFF-cone and rod bipolar',
                        '20' = 'Amacrine',
                        '21' = 'ON-cone and rod bipolar',
                        '22' = 'doublets',
                        '23' = 'GABA-amacrine',
                        '24' = 'doublets',
                        '25' = 'Müller glia',
                        '26' = 'GABA-amacrine',
                        '27' = 'Müller glia',
                        '28' = 'Microglia',
                        '29' = 'ON-cone and rod bipolar',
                        '30' = 'OFF-cone and rod bipolar',
                        '31' = 'ON-cone and rod bipolar',
                        '32' = 'ON-cone and rod bipolar',
                        '33' = 'OFF-cone and rod bipolar',
                        '34' = 'doublets',
                        '35' = 'doublets',
                        '36' = 'doublets',
                        '37' = 'doublets'
                        )
# find all markers od cluster "XXX"
# XXX.markers <- FindMarkers(macaque.combined, ident.1 = "24",ident.2 = "20")
# View(XXX.markers)

macaque.combined$celltype <- Idents(macaque.combined)
macaque.combined <- subset(macaque.combined, idents=c("Bipolar","ON-cone and rod bipolar","OFF-cone and rod bipolar","Müller glia",
                "Microglia","Retinal ganglion cell","Gly-amacrine","GABA-amacrine","Amacrine","Horizontal",
                "Cone","Rod"))

p1 <- DimPlot(macaque.combined, reduction = "umap", group.by = "celltype", label = FALSE, label.size = 2.5, repel = TRUE, split.by = "Samples")
ggsave("./result/2.samples_integration/scRNA_integration/all_samples_anno_label_dimplot.pdf", p1, width = 30,height = 7)

levels(macaque.combined)=c("Bipolar","ON-cone and rod bipolar","OFF-cone and rod bipolar","Müller glia",
                "Microglia","Retinal ganglion cell","Gly-amacrine","GABA-amacrine","Amacrine","Horizontal",
                "Cone","Rod")
p2 <- DotPlot(macaque.combined,features = cd_genes)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave("./result/2.samples_integration/scRNA_integration/celltype_anno_dotplot.pdf", p2, width = 15, height = 7)

p3 <- DimPlot(macaque.combined, reduction = "umap", group.by = "celltype", pt.size = .1, label = TRUE)
ggsave("./result/2.samples_integration/scRNA_integration/plot_label_Harmony_dimplot.pdf", p3, width = 9,height = 7)

saveRDS(macaque.combined,"./result/2.samples_integration/scRNA_integration/scRNA_integrated.rds")

# save RNA counts and metadata
rna_count <- macaque.combined[["RNA"]]@counts
metadata <- macaque.combined@meta.data

saveRDS(rna_count, "./result/2.samples_integration/scRNA_integration/scRNA_integrated_counts.rds")
saveRDS(metadata, "./result/2.samples_integration/scRNA_integration/scRNA_integrated_metadata.rds")
