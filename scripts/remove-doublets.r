### pre-processing data  
removeDoublet.cellNumber <- matrix(rep(0,12),10,3)
rownames(removeDoublet.cellNumber) <- c("E-OS-C","A-OS-C","H-OS-C","D-OS-P","A-OS-P","H-OS-P","COS-H",,"COS-H","EOS-Z","IOS-H","IOS-Z")
colnames(removeDoublet.cellNumber) <- c("raw_cell_numbers","cell_numbers_after_qc","cell_numbers_after_removedoublets")
for(i in 1:10){
  # samples=c("E-OS-C","A-OS-C","H-OS-C","D-OS-P","A-OS-P","H-OS-P","COS-H","EOS-Z","IOS-H","IOS-Z")
  samples=c("EOS-Z","IOS-H","IOS-Z")
  sample <- samples[i]
  print(paste0("Now processing: ",sample))

  # h5_file <- paste0("./cellrangerARC_output/",sample,"/outs/filtered_feature_bc_matrix.h5")
  h5_file <- paste0("./result/0.cellranger_count/",sample,"/filtered_feature_bc_matrix.h5")
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
  ggsave(paste0("./result/1.remove_doublets/",sample,"_DoubletFinder.pdf"),p1,width = 10,height = 10)

  # filter doublets
  Idents(macaque ) <- DF.clarrsications
  macaque<-subset(x = macaque, idents="Singlet")   #过滤细胞，只保留经过检测的Singlet
  filtered_cell_number <- dim(macaque)[2]

  removeDoublet.cellNumber[i,1] = raw_cell_number
  removeDoublet.cellNumber[i,2] = cell_numbers_after_qc
  removeDoublet.cellNumber[i,3] = filtered_cell_number

  print(paste0("raw cell numbers: ",raw_cell_number))
  print(paste0("cell numbers after quality control: ", cell_numbers_after_qc))
  print(paste0("filtered cell numbers: ",filtered_cell_number))

  # save filtered rna_counts 
  filtered_barcodes <- colnames(macaque[["RNA"]]@counts)
  filtered.rna.counts = rna_counts[,filtered_barcodes]
  colnames(filtered.rna.counts) = paste0(sample,".",colnames(filtered.rna.counts))
  write.table(filtered.rna.counts, paste0("./result/1.remove_doublets/",sample,"/",sample,"_filtered.rna.counts.txt"),row.names = T,sep = "\t",quote = F)

  # rm(list = ls())
  gc()
}
write.table(removeDoublet.cellNumber,"./result/1.remove_doublets/removeDoublet.cellNumber.txt",sep="\t",quote=F,row.names = T)
