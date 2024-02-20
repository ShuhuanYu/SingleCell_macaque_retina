library(GSVA)
library(readxl)
library(limma)

scRNA <- readRDS("./result/2.samples_integration/scRNA_integration/scRNA_integrated.rds")
expr <- as.matrix(scRNA@assays[["RNA"]]@data)

geneSets <- read.table("./result/4.aging_pathways_analysis/gene_sets_interested.txt", header = T, sep = "\t", stringsAsFactors = F)
gs <- list()
for(i in 1:dim(geneSets)[1]){
    gs[[i]] = unlist(strsplit(geneSets[i,2], ","))
}
names(gs) <- geneSets[,1]

gsvaScore <- gsva(expr = expr[intersect(rownames(expr), unlist(gs)),], gset.idx.list = gs, 
                    method = "gsva", 
                    kcdf = "Gaussian",
                    min.sz = 5,
                    max.sz = 500,
                    mx.diff = TRUE,
                    abs.ranking = FALSE
                    )

saveRDS(gsvaScore, "./result/4.aging_pathways_analysis/aging_and_mitochondrial_gsvaScore.rds")
write.table(gsvaScore, "./result/4.aging_pathways_analysis/aging_and_mitochondrial_gsvaScore.txt")
