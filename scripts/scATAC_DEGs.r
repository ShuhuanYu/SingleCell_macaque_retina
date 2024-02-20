library(ArchR)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(TxDb.Mmulatta.UCSC.rheMac10.refGene)
library(org.Mmu.eg.db)
library(xlsx)
library(dplyr)
library(chromVAR)
library(ComplexHeatmap)
library(circlize)
library(colorRampPalette)
library(corrplot)
library(Hmisc)


proj <- loadArchRProject(path = "/opt/yushuhuan/macaque/result/2.samples_integration/scATAC_integration/macaqueRetinaArchR")
getAvailableMatrices(proj)

state <- proj$Sample
state <- gsub("E_OS_C", "youngC", state)
state <- gsub("COS_H", "youngC", state)
state <- gsub("D_OS_P", "youngP", state)
state <- gsub("EOS_Z", "youngP", state)
state <- gsub("H_OS_C", "oldC", state)
state <- gsub("IOS_H", "oldC", state)
state <- gsub("H_OS_P", "oldP", state)
state <- gsub("IOS_Z", "oldP", state)
state <- gsub("A_OS_C", "middleC", state)
state <- gsub("A_OS_P", "middleP", state)
table(state)
# middleC middleP    oldC    oldP  youngC  youngP 
#    9095    8564   18176   18022   14058   17265
proj$state <- state
head(getCellColData(proj))

# Pairwise Testing Between Groups
compare1 <- c("oldC", "youngC")
compare2 <- c("oldP", "youngP")

# compare <- compare1
compare <- compare2

cell_clusters <- unique(getCellColData(proj)$predictedGroup)

for(celltype in cell_clusters){
    idxSample <- BiocGenerics::which(proj$predictedGroup %in% celltype)
    cellsSample <- proj$cellNames[idxSample]
    proj_celltype <- proj[cellsSample, ]

    markerTest <- getMarkerFeatures(
                                    ArchRProj = proj_celltype, 
                                    useMatrix = "PeakMatrix",
                                    groupBy = "state",
                                    testMethod = "wilcoxon",
                                    bias = c("TSSEnrichment", "log10(nFrags)"),
                                    useGroups = compare[1],
                                    bgdGroups = compare[2]
                                    )
    
    # visualize differential peaks
    pma <- plotMarkers(seMarker = markerTest, name = compare[1], cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")
    plotPDF(pma, name = paste0(compare[1], "-vs-", compare[2], "-", celltype, "-Markers-MA-plot"), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

    # Motif Enrichment in Differential Peaks
    # UP differeatial peaks 
    motifsUp <- peakAnnoEnrichment(
                                    seMarker = markerTest,
                                    ArchRProj = proj_celltype,
                                    peakAnnotation = "Motif",
                                    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
                                    )
    df.up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df.up <- df.up[order(df.up$mlog10Padj, decreasing = TRUE),]
    df.up$rank <- seq_len(nrow(df.up))
    ggUp <- ggplot(df.up, aes(rank, mlog10Padj, color = mlog10Padj)) + 
                    geom_point(size = 1) +
                    ggrepel::geom_label_repel(
                            data = df.up[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
                            size = 1.5,
                            nudge_x = 2,
                            color = "black"
                    ) + theme_ArchR() + 
                    ylab("-log10(P-adj) Motif Enrichment") + 
                    xlab("Rank Sorted TFs Enriched") +
                    scale_color_gradientn(colors = paletteContinuous(set = "comet")) + 
                    ggtitle("up peaks motif enrichment result")

    # Down differential peaks
    motifsDo <- peakAnnoEnrichment(
                                    seMarker = markerTest,
                                    ArchRProj = proj_celltype,
                                    peakAnnotation = "Motif",
                                    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
                                    )
    df.do <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
    df.do <- df.do[order(df.do$mlog10Padj, decreasing = TRUE),]
    df.do$rank <- seq_len(nrow(df.do))
    ggDo <- ggplot(df.do, aes(rank, mlog10Padj, color = mlog10Padj)) + 
                    geom_point(size = 1) +
                    ggrepel::geom_label_repel(
                            data = df.do[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
                            size = 1.5,
                            nudge_x = 2,
                            color = "black"
                    ) + theme_ArchR() + 
                    ylab("-log10(FDR) Motif Enrichment") +
                    xlab("Rank Sorted TFs Enriched") +
                    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
                    ggtitle("down peaks motif enrichment result")

    plotPDF(ggUp, ggDo, name = paste0(compare[1], "-vs-", compare[2], "-", celltype, "-Markers-Motifs-Enriched"), width = 5, height = 5, ArchRProj = proj_celltype, addDOC = FALSE)

    # output
    peaks.up <- getMarkers(markerTest, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
    peaks.do <- getMarkers(markerTest, cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
    peaks.celltype <- as.data.frame(rbind(peaks.up[[1]], peaks.do[[1]]))
    peaks.celltype$change <- rep(c("UP", "DOWN"), c(lengths(peaks.up), lengths(peaks.do)))
    write.xlsx(peaks.celltype, paste0("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/all.celltypes.", 
                                compare[1], "_", compare[2],
                                "_differential_peaks.xlsx"), sheetName = celltype, row.names = F, append = TRUE)

    motifs.up.enrichment <- matrix(as.data.frame(assays(motifsUp))$value, dim(motifsUp)[1], 10)
    colnames(motifs.up.enrichment) <- unique(as.data.frame(assays(motifsUp))$group_name)
    rownames(motifs.up.enrichment) <- rownames(motifsUp)
    motifs.do.enrichment <- matrix(as.data.frame(assays(motifsDo))$value, dim(motifsDo)[1], 10)
    colnames(motifs.do.enrichment) <- unique(as.data.frame(assays(motifsDo))$group_name)
    rownames(motifs.do.enrichment) <- rownames(motifsDo)
    motif.enrichment.re <- as.data.frame(rbind(motifs.up.enrichment, motifs.do.enrichment))
    motif.enrichment.re$change <- rep(c("UP", "DOWN"), c(dim(motifs.up.enrichment)[1], dim(motifs.do.enrichment)[1]))
    write.xlsx(motif.enrichment.re, paste0("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/all.celltypes.",
                                        compare[1], "_", compare[2],
                                        "_differential_peaks_motif_enrichment.xlsx"),
                sheetName = celltype, row.names = TRUE, append = TRUE)

    motif.enrichment.re.sig <- motif.enrichment.re %>% filter(mlog10Padj > 1) %>% arrange(desc(change))
    write.xlsx(motif.enrichment.re.sig, paste0("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/all.celltypes.",
                                        compare[1], "_", compare[2],
                                        "_differential_peaks_motif_enrichment_sig.xlsx"),
                                        sheetName = celltype, row.names = FALSE, append = TRUE)

}

### chromVAR deviation
proj1 <- addBgdPeaks(proj)
proj1 <- addDeviationsMatrix(
                            ArchRProj = proj1, 
                            peakAnnotation = "Motif",
                            force = FALSE,
                            binarize = TRUE
                            )
saveArchRProject(ArchRProj = proj1, outputDirectory = "./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif", load = TRUE)

proj1 <- loadArchRProject("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif")

motifMatrix <- getMatrixFromProject(ArchRProj = proj1,  # extract chromVAR result
                                        useMatrix = "MotifMatrix"
                                        )
saveRDS(motifMatrix, "./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/motifMatrix.rds")

zscoreMatrix <- assay(motifMatrix, "z") # get motif accessibility zscore
rownames(zscoreMatrix) <- apply(data.frame(rownames(zscoreMatrix)), 1, function(x){
                motifName = strsplit(x[1], "_")[[1]][1]
        })

scRNA <- readRDS("./result/2.samples_integration/scRNA_integration/scRNA_integrated.rds")
zscoreRNA <- scRNA[["RNA"]]@data # get RNA z-score

all.motif.names <- apply(data.frame(rownames(zscoreMatrix)), 1, function(x){
                motifName = strsplit(x[1], "_")[[1]][1]
        }) %>% intersect(rownames(zscoreRNA))

zscoreMatrix <- zscoreMatrix[all.motif.names, ]
zscoreRNA <- zscoreRNA[all.motif.names, ]

scRNA.meta <- readRDS("./result/2.samples_integration/scRNA_integration/scRNA_integrated_metadata.rds")
state <- scRNA.meta$Samples
state <- gsub("E-OS-C", "youngC", state)
state <- gsub("COS-H", "youngC", state)
state <- gsub("D-OS-P", "youngP", state)
state <- gsub("EOS-Z", "youngP", state)
state <- gsub("H-OS-C", "oldC", state)
state <- gsub("IOS-H", "oldC", state)
state <- gsub("H-OS-P", "oldP", state)
state <- gsub("IOS-Z", "oldP", state)
state <- gsub("A-OS-C", "middleC", state)
state <- gsub("A-OS-P", "middleP", state)
table(state)
scRNA.meta$state <- state

scATAC.meta <- as.data.frame(getCellColData(proj1))

table(scRNA.meta$state, scRNA.meta$celltype)
table(scATAC.meta$state, scATAC.meta$predictedGroup)

for(cell in cell_clusters){
        celltype.zscore.RNA <- zscoreRNA[, rownames(scRNA.meta %>% filter(state %in% c(compare[1], compare[2])) %>% filter(celltype == cell))]
        celltype.zscore.accessibility <- zscoreMatrix[, rownames(scATAC.meta %>% filter(state %in% c(compare[1], compare[2])) %>% filter(predictedGroup == cell))]

        celltype.motif.sig <- read.xlsx(paste0("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/all.celltypes.",
                                        compare[1], "_", compare[2],
                                        "_differential_peaks_motif_enrichment_sig.xlsx"), sheetName = cell)
        if(dim(celltype.motif.sig)[1] > 1){
                motif.names <- apply(celltype.motif.sig, 1, function(x){
                        motifName = strsplit(x[10], "_")[[1]][1]
                }) %>% intersect(rownames(zscoreRNA))

                # heatmap
                ATAC.top_anno <- HeatmapAnnotation(condition = scATAC.meta[colnames(celltype.zscore.accessibility), 20],
                                                        col = list(condition = c("oldP" = "#d76d43", "youngP" = "#2e99c7")))
                column_order.ATAC <- rownames(scATAC.meta[colnames(celltype.zscore.accessibility), ] %>% arrange(state))
                hm.ATAC <- Heatmap(as.matrix(celltype.zscore.accessibility[motif.names, ]), 
                                        top_annotation = ATAC.top_anno,
                                        width = unit(8, "cm"), height = unit(10, "cm"),
                                        cluster_columns = FALSE,
                                        cluster_rows = TRUE, 
                                        show_row_dend = FALSE,
                                        show_column_dend = FALSE,
                                        show_column_names = FALSE,
                                        col = colorRamp2(c(0, 0.1, 1), c("blue", "white", "red")),
                                        # col = colorRampPalette(colors = c("blue","white","red"))(5),
                                        heatmap_legend_param = list(
                                                title = "Motif accessibility z-score", at = c(-2, -1, 0, 1, 2)
                                        ),
                                        column_order = column_order.ATAC,
                                        row_names_gp = gpar(fontsize = 4)
                                        )
                RNA.top_anno <- HeatmapAnnotation(condition = scRNA.meta[colnames(celltype.zscore.RNA), 8], 
                                                        col = list(condition = c("oldP" = "#d76d43", "youngP" = "#2e99c7")))
                column_order.RNA <- rownames(scRNA.meta[colnames(celltype.zscore.RNA), ] %>% arrange(state))
                hm.RNA <- Heatmap(as.matrix(celltype.zscore.RNA[motif.names, ]), 
                                        top_annotation = RNA.top_anno,
                                        width = unit(8, "cm"), height = unit(10, "cm"),
                                        cluster_columns = FALSE,
                                        cluster_rows = TRUE,
                                        show_row_dend = FALSE,
                                        show_column_dend = FALSE,
                                        show_column_names = FALSE,
                                        col = colorRamp2(c(0, 0.1, 1), c("#ed4dd0", "black", "#e8e85f")),
                                        heatmap_legend_param = list(
                                                title = "RNA z-score", at = c(-2, -1, 0, 1, 2)
                                        ),
                                        column_order = column_order.RNA,
                                        row_names_gp = gpar(fontsize = 4)
                                        )
                hmlist <- hm.ATAC + hm.RNA
                pdf(paste0("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/", compare[1],"_vs_", compare[2], "_motif_heatmap/",
                                        cell, "_motif_accessibility_RNA_zscore_heatmap.pdf"),
                                        width = 15, height = 8
                                        )
                draw(hmlist, ht_gap = unit(3, "cm"))
                dev.off()

                # spearman and pearson
                celltype.motifs.RNA.cor <- data.frame(spearman.cor = rep(0, length(motif.names)), spearman.pvalue = rep(0, length(motif.names)),
                                                        pearson.cor = rep(0, length(motif.names)), pearson.pvalue = rep(0, length(motif.names))
                                                        )
                rownames(celltype.motifs.RNA.cor) <- motif.names  
                for(motif in motif.names){
                        atac = celltype.zscore.accessibility[motif, ]
                        rna = celltype.zscore.RNA[motif, ]

                        re.spearman = rcorr(t(rbind(as.double(atac), as.double(rna))), type = "spearman")
                        re.pearson = rcorr(t(rbind(as.double(atac), as.double(rna))), type = "pearson")
                        
                        celltype.motifs.RNA.cor[motif, 1] = re.spearman$r[1,2]
                        celltype.motifs.RNA.cor[motif, 2] = re.spearman$P[1,2]
                        celltype.motifs.RNA.cor[motif, 3] = re.pearson$r[1,2]
                        celltype.motifs.RNA.cor[motif, 4] = re.pearson$P[1,2]
                } 
                write.xlsx(celltype.motifs.RNA.cor, paste0("./result/3.senescence_related_DEGs/scATAC_DE_peaks.motif/",
                                compare[1], "_vs_", compare[2], "_celltype_motif_RNA_cor.xlsx"),
                                sheetName = cell, row.names = TRUE, append = TRUE
                                )   
        }    
}







###############################################################################################################################################
library(ChIPseeker)
