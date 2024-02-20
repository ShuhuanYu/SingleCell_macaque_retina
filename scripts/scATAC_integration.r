setwd("./result/2.samples_integration/scATAC_integration")

library(ArchR)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(TxDb.Mmulatta.UCSC.rheMac10.refGene)
library(org.Mmu.eg.db)
library(magrittr)
library(tidyverse)
library(SummarizedExperiment)
library(Seurat)
library(SeuratObject)
set.seed(1)

inputFiles <- c("../../0.cellranger_count/A-OS-C/outs/atac_fragments.tsv.gz",
                "../../0.cellranger_count/A-OS-P/outs/atac_fragments.tsv.gz",
                "../../0.cellranger_count/COS-H/atac_fragments.tsv.gz",
                "../../0.cellranger_count/D-OS-P/outs/atac_fragments.tsv.gz",
                "../../0.cellranger_count/E-OS-C/outs/atac_fragments.tsv.gz",
                "../../0.cellranger_count/EOS-Z/atac_fragments.tsv.gz",
                "../../0.cellranger_count/H-OS-C/outs/atac_fragments.tsv.gz",
                "../../0.cellranger_count/H-OS-P/outs/atac_fragments.tsv.gz",
                "../../0.cellranger_count/IOS-H/atac_fragments.tsv.gz",
                "../../0.cellranger_count/IOS-Z/atac_fragments.tsv.gz"
                )
names(inputFiles) <- c("A_OS_C", "A_OS_P", "COS_H", "D_OS_P", "E_OS_C", "EOS_Z", "H_OS_C", "H_OS_P", "IOS_H", "IOS_Z")

# addArchRThreads(threads = 16) 

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmulatta.UCSC.rheMac10)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Mmulatta.UCSC.rheMac10.refGene, OrgDb = org.Mmu.eg.db)
loci <- grep("NA", geneAnnotation$genes$symbol)
gid <- geneAnnotation$genes$gene_id[-loci]
df <- AnnotationDbi::select(TxDb.Mmulatta.UCSC.rheMac10.refGene, keys = gid, columns="TXNAME", keytype="GENEID")
genes <- geneAnnotation$genes[-loci]
exons <- geneAnnotation$exons[-grep("NA", geneAnnotation$exons$symbol)]
tss <- geneAnnotation$TSS[which(geneAnnotation$TSS$tx_name %in% df$TXNAME)]
geneAnnotation <- createGeneAnnotation(genes = genes, 
                                             exons = exons, 
                                             TSS = tss)
# save(genomeAnnotation, geneAnnotation, file="./Macaca_mulatta_genomeAnnotation_geneAnnotationSubset.RData")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation,
  force = FALSE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
    force = FALSE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "macaqueRetinaArchR",
  copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation
)

proj <- filterDoublets(proj)
# Filtering 9713 cells from ArchRProject!
#         H_OS_P : 1422 of 12888 (11%)
#         IOS_H : 528 of 7270 (7.3%)
#         E_OS_C : 1690 of 13000 (13%)
#         H_OS_C : 1733 of 13167 (13.2%)
#         EOS_Z : 864 of 9297 (9.3%)
#         IOS_Z : 522 of 7230 (7.2%)
#         A_OS_C : 1023 of 10118 (10.1%)
#         COS_H : 79 of 2827 (2.8%)
#         A_OS_P : 894 of 9458 (9.5%)
#         D_OS_P : 958 of 9790 (9.8%)
saveArchRProject(ArchRProj = proj, outputDirectory = "macaqueRetinaArchR", load = TRUE)

## clusring by tile matrix
proj <- addIterativeLSI(
    # first round iterative LSI based on tilematrix, with default para, it will carry out estimated LSI
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    force = TRUE
  )
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
  )
proj <- addClusters(
    # add cluster based on LSI using seurat
    input = proj,
    reducedDims = "IterativeLSI",
    name = "Clusters",
    force = TRUE
  )
proj <- addUMAP(
    # add embedding
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    force = TRUE
  )
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAPHarmony",
  nNeighbors = 30
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = "macaqueRetinaArchR", load = TRUE)

## Cross-platform linkage of scATAC-seq cells with scRNA-seq cells
proj <- loadArchRProject(path = "./macaqueRetinaArchR")
getAvailableMatrices(proj)

counts <- readRDS("../scRNA_integration/scRNA_integrated_counts.rds")
meta <- readRDS("../scRNA_integration/scRNA_integrated_metadata.rds")
scRNA <- SummarizedExperiment(assays = list(counts = counts), colData = meta)

# Unconstrained Integration
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = scRNA,
    addToArrow = FALSE,
    groupRNA = "celltype",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
pal <- paletteDiscrete(values = colData(scRNA)$celltype)
p <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup", 
    embedding = "UMAPHarmony",
    pal = pal
)
plotPDF(p, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = proj, outputDirectory = "macaqueRetinaArchR", load = TRUE)

## marker Genes Imputation with MAGIC
proj <- addImputeWeights(proj)

################################## 细胞注释完成 #####################################

## Making Pseudo-bulk Replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "predictedGroup",
  minCells = 15,
  maxCells = 700, 
  force = TRUE
)
## call peaks
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "predictedGroup",
  genomeSize = 2.7e09,
  pathToMacs2 = "/opt/softwares/anaconda3/bin/macs2",
  force = TRUE
)
## calculate peak matrix
proj <- addPeakMatrix(proj, binarize = TRUE, force = TRUE)

## Identifying Marker Peaks with ArchR
table(proj$predictedGroup)
#                 Amacrine                  Bipolar                     Cone 
#                       30                       15                     2337 
#            GABA-amacrine             Gly-amacrine               Horizontal 
#                     3133                     1961                      992 
#                Microglia              Müller glia OFF-cone and rod bipolar 
#                       88                     6596                     4449 
#  ON-cone and rod bipolar    Retinal ganglion cell                      Rod 
#                     5416                     1530                    58785

saveArchRProject(ArchRProj = proj, outputDirectory = "macaqueRetinaArchR", load = TRUE)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

# Marker Peak Heatmaps
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
# draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

## Motif and Feature Enrichment with ArchR
# Motif Enrichment in Marker Peaks
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", species = "homo sapiens")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj, outputDirectory = "macaqueRetinaArchR", load = TRUE)
