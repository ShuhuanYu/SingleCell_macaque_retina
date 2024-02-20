library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(org.Mmu.eg.db)
library(clusterProfiler)
library(xlsx)
library(enrichplot)
library(scRNAtoolVis)
library(RColorBrewer)
library(patchwork)
library(scales)

## functions ##
GOall_enrichment <- function(degs_result, change=c("UP", "DOWN")){
  DGEs.change <- degs_result %>% filter(Change == change)
  eg <- bitr(rownames(DGEs.change), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mmu.eg.db)

  # GO enrichment analysis
  go.all <-enrichGO(gene = eg$ENTREZID, 
            keyType = 'ENTREZID', 
            OrgDb = org.Mmu.eg.db, 
            ont = 'ALL', 
            pAdjustMethod = 'BH', 
            pvalueCutoff = 0.01,
            qvalueCutoff = 0.05,
            readable=T,
            minGSSize = 10,
            maxGSSize = 500
            )
  if(dim(go.all)[1] == 0){
    # p <- barplot(go.all, showCategory = dim(go.all)[1]) + ggtitle(paste("scRNA GO ALL", compare[1], "vs.", compare[2], "(p.adjust<0.01 & qvalue <0.05)")) # show all results
    # ggsave(paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/allDGEs.", change, ".", compare[1], ".vs.", compare[2], "-GO-enrichment.pdf"), p, width = 10, height = 15)
    go.all <-enrichGO(gene = eg$ENTREZID, 
          keyType = 'ENTREZID', 
          OrgDb = org.Mmu.eg.db, 
          ont = 'ALL', 
          pAdjustMethod = 'BH', 
          pvalueCutoff = 0.05, 
          qvalueCutoff = 0.05,
          readable=T,
          minGSSize = 10,
          maxGSSize = 500
          )
  if(dim(go.all)[1] == 0){
    # p <- barplot(go.all, showCategory = dim(go.all)[1]) + ggtitle(paste("scRNA GO ALL", compare[1], "vs.", compare[2], "(p.adjust<0.05, qvalue<0.05)")) # show all results
    # ggsave(paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/allDGEs.", change, ".", compare[1], ".vs.", compare[2], "-GO-enrichment.pdf"), p, width = 10, height = 15)
    go.all <-enrichGO(gene = eg$ENTREZID, 
        keyType = 'ENTREZID', 
        OrgDb = org.Mmu.eg.db, 
        ont = 'ALL', 
        pAdjustMethod = 'BH', 
        pvalueCutoff = 0.05, 
        qvalueCutoff = 0.05,
        readable=T,
        minGSSize = 1,
        maxGSSize = 5000
        )
    # p <- barplot(go.all, showCategory = dim(go.all)[1]) + ggtitle(paste("scRNA GO ALL", compare[1], "vs.", compare[2], "(p.adjust<0.05, qvalue<0.05 & GSSsize in 1:5000)")) # show all results
  }
  }

  return(go.all)
}
## ##

macaque.combined <- readRDS("./result/2.samples_integration/scRNA_integration/scRNA_integrated.rds")

clusters <- levels(macaque.combined$celltype)
Idents(macaque.combined) <- macaque.combined$state

compare1 <- c("oldC", "youngC")
compare2 <- c("oldP", "youngP")

## set compare condition
compare <- compare1
# compare <- compare2

###############################################################################################################
## sample levels DGEs
all.cells.DGEs <- FindMarkers(macaque.combined, ident.1 = compare[1], ident.2 = compare[2], verbose = FALSE, avg_log2FC.threshold = 0, min.pct = 0, min.cells.feature = 1, min.cells.group = 1)

# vocalno plot for global gene signature
all.cells.DGEs$Change = as.factor(ifelse(all.cells.DGEs$p_val_adj < 0.01 & abs(all.cells.DGEs$avg_log2FC) > 0.5,
                              ifelse(all.cells.DGEs$avg_log2FC > 0.5 ,'UP','DOWN'),'STABLE'))
all.cells.DGEs$gene <- rownames(all.cells.DGEs)
all.cells.DGEs <- all.cells.DGEs %>% filter(grepl("ENSMMUG", gene) == "FALSE") # 移除没有gene symbol的feature
table(all.cells.DGEs$Change)
# oldC vs youngC: 198 down, 672 up, 11353 stable
# oldP vs youngP: 356 down, 1161 up, 11556 stable 
# specify th gene labels that you want to show
all.cells.DGEs$label <- ifelse(all.cells.DGEs$p_val_adj< 0.00001& abs(all.cells.DGEs$avg_log2FC) >= 2.5, all.cells.DGEs$gene,"")
table(all.cells.DGEs$label)
# another choice
# show.labels <- c("MASP1", "SLC35F4")
# all.cells.DGEs$label <- ifelse(all.cells.DGEs$gene %in% show.labels, all.cells.DGEs$gene, "")

p1 <- ggplot(all.cells.DGEs, aes(x=avg_log2FC, y=-log10(p_val_adj),color=Change)) + 
      geom_point(alpha=0.4, size=2) + 
      theme_bw(base_size = 12) + 
      xlab("Log2(Fold change)") +
      ylab("-Log10(P.adj)") +
      theme(plot.title = element_text(size=15,hjust = 0.5)) + 
      scale_colour_manual(values = c('#3685c7','#302e2e','brown')) +
      geom_hline(yintercept = -log10(0.01), lty = 4) +
      geom_vline(xintercept = c(-0.5, 0.5), lty = 4)+
      labs(title = paste0("scRNA DGEs ", compare[1]," vs. ", compare[2]))+
      geom_label_repel(data = all.cells.DGEs, aes(label = label),
                      size = 3,box.padding = unit(0.5, "lines"),
                      point.padding = unit(0.8, "lines"),
                      segment.color = "black",
                      show.legend = FALSE, max.overlaps = 10000)
ggsave(paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/vocalno_plot_of_DGEs_", compare[1], "_vs_", compare[2], ".pdf"), p1)

write.csv(all.cells.DGEs, paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/allDGEs.", compare[1], ".vs.", compare[2], ".csv"), row.names = TRUE, quote = F)
saveRDS(all.cells.DGEs, paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/allDGEs.", compare[1], ".vs.", compare[2], ".rds"))

# enrichment analysis
go.all.up <- GOall_enrichment(all.cells.DGEs, change = "UP")
p1 <- barplot(go.all.up, showCategory = dim(go.all.up)[1]) + ggtitle(paste("scRNA GO ALL UP", compare[1], "vs.", compare[2]))

go.all.do <- GOall_enrichment(all.cells.DGEs, change = "DOWN")
p2 <- barplot(go.all.do, showCategory = dim(go.all.do)[1]) + ggtitle(paste("scRNA GO ALL DOWN", compare[1], "vs.", compare[2]))

ggsave(paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/samples_levels_", compare[1], ".vs", compare[2],"_DGEs_GO_all_enrichment.pdf"), p1+p2)

###############################################################################################################

## celltypes level
all.celltypes.DGEs <- data.frame()
for(i in 1:12){
  macaque.celltype <- subset(macaque.combined, celltype == clusters[i])
  # to get global gene signature of retina cells related with aging
  DGEs <- FindMarkers(macaque.celltype, ident.1 = compare[1], ident.2 = compare[2], 
                      verbose = FALSE, avg_log2FC.threshold = 0, min.pct = 0.1)
  dim(DGEs)
  DGEs$gene <- rownames(DGEs)
  DGEs <- DGEs %>% filter(grepl("ENSMMU",gene)=="FALSE")
  DGEs$cluster <- clusters[i]

  all.celltypes.DGEs <- rbind(all.celltypes.DGEs, DGEs)

  DGEs_sig <- DGEs %>% filter(p_val_adj<0.01 & abs(avg_log2FC)>0.5) %>% arrange(desc(avg_log2FC))
  # write.xlsx(DGEs_sig %>% select(-cluster), paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/celltype-DGEs-",compare[1],".vs.",compare[2],"-findmarkers.xlsx"), sheetName = clusters[i], row.names = T, append = TRUE)

}
all.celltypes.DGEs$sig=""
all.celltypes.DGEs$sig[abs(all.celltypes.DGEs$avg_log2FC) > 0.5 & all.celltypes.DGEs$p_val_adj < 0.01] = "sig"
all.celltypes.DGEs$sig2=paste(all.celltypes.DGEs$cluster,all.celltypes.DGEs$sig,sep = "_")
all.celltypes.DGEs$sig2[str_detect(all.celltypes.DGEs$sig2,"_$")]="not_sig"
all.celltypes.DGEs$sig2=str_replace(all.celltypes.DGEs$sig2,"_sig","")

#控制顺序
all.celltypes.DGEs$sig2=factor(all.celltypes.DGEs$sig2,levels = c("not",sort(unique(all.celltypes.DGEs$cluster))))
all.celltypes.DGEs$cluster=factor(all.celltypes.DGEs$cluster,levels = sort(unique(all.celltypes.DGEs$cluster)))
all.celltypes.DGEs=all.celltypes.DGEs%>%arrange(cluster,sig2)
#控制范围
all.celltypes.DGEs$avg_log2FC[all.celltypes.DGEs$avg_log2FC > 3]=3
all.celltypes.DGEs$avg_log2FC[all.celltypes.DGEs$avg_log2FC < c(-3)]= -3

## show all clusters DEGs in vocalno plot
color_ct=c(brewer.pal(12, "Set3")[-c(2,9)],
           brewer.pal(5, "Set1")[2],
           brewer.pal(3, "Dark2")[1])
names(color_ct)=sort(unique(as.character(all.celltypes.DGEs$cluster)))

# all.celltypes.DGEs %>% ggplot(aes(x=cluster,y=avg_log2FC,color=sig2))+geom_jitter(width = 0.25,size=0.5)+
#   scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
#   scale_y_continuous("Arep VS Brep, average log2FC",expand = c(0.02,0))+
#   theme_bw()+
#   theme(
#     panel.grid = element_blank(),
#     legend.position = "none",
#     axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
#     axis.text.y.left = element_text(size = 14,color = "black"),
#     axis.title.x.bottom = element_blank(),
#     axis.title.y.left = element_text(size = 16)
#   )
all.celltypes.DGEs$padj_log10_neg= -log10(all.celltypes.DGEs$p_val_adj)
all.celltypes.DGEs$padj_log10_neg=ifelse(all.celltypes.DGEs$avg_log2FC > 0,
                                        all.celltypes.DGEs$padj_log10_neg,
                                        -all.celltypes.DGEs$padj_log10_neg)
plot.list=list()
for (ci in sort(unique(as.character(all.celltypes.DGEs$cluster)))) {
  tmpdf=all.celltypes.DGEs %>% filter(cluster == ci)
  minabs=abs(min(tmpdf$padj_log10_neg))
  maxabs=max(tmpdf$padj_log10_neg)
  thre=0
  if(minabs < maxabs) {
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > minabs] = minabs
    thre=minabs
  }
  if(minabs > maxabs) {
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-maxabs)] = -maxabs
    thre=maxabs
  }
  if(minabs == maxabs & maxabs == Inf) {
    thre = min(
      abs(
        range(
          tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < Inf & tmpdf$padj_log10_neg > -Inf]
        )
      )
    )
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg < (-thre)] = -thre
    tmpdf$padj_log10_neg[tmpdf$padj_log10_neg > thre] = thre
  }
  
  plotdata = tmpdf
  tmpdf=tmpdf%>%filter(sig2 != "not") #这里我取了logFC最极端的几个gene来标注文本，实际处理中不一定这样做
  tmpdf=tmpdf%>%arrange(desc(avg_log2FC))
  tmpdf.a=head(tmpdf%>%filter(avg_log2FC > 0),5)
  tmpdf.a$d=thre*2*0.05+(-thre)-tmpdf.a$padj_log10_neg
  tmpdf.b=tail(tmpdf%>%filter(avg_log2FC < 0),5)
  tmpdf.b$d=thre*2*0.95-thre  - tmpdf.b$padj_log10_neg
  textdata.down = tmpdf.b
  textdata.up   = tmpdf.a
  
  ###画图
  tmpplot=plotdata%>%ggplot(aes(x=padj_log10_neg,y=avg_log2FC))+
    geom_point(aes(color=sig2),size=1)+
    geom_hline(yintercept = c(-0.5,0.5),linetype="dashed")+
    geom_text_repel(data = textdata.down,
                    mapping = aes(label=gene),
                    nudge_x=textdata.down$d,
                    direction = "y", hjust = 1,segment.size = 0.2)+
    geom_text_repel(data = textdata.up,
                    mapping = aes(label=gene),
                    nudge_x=textdata.up$d,
                    direction = "y", hjust = 0,segment.size = 0.2)+
    labs(title = ci)+
    scale_color_manual(values = c(color_ct,"not"="#dee1e6"))+
    scale_y_continuous("Arep VS Brep, average log2FC",expand = c(0.02,0),limits = c(-3,3))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      
      axis.ticks.x.bottom = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = 14,color = "black"),
      axis.title.y.left = element_text(size = 16),
      
      plot.title = element_text(size = 16,hjust = 0.5)
    )
  
  index=which(ci == sort(unique(as.character(all.celltypes.DGEs$cluster))))
  if (index!=1) {
    tmpplot=tmpplot+theme(
      axis.title.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y.left = element_blank()
    )
  }
  if (index == length(sort(unique(as.character(all.celltypes.DGEs$cluster))))) {
    segment.df=data.frame(x=c(0 - thre / 5,0 + thre / 5),
                          xend=c(-thre,thre),
                          y=c(-3,-3),
                          yend=c(-3,-3))
    tmpplot=tmpplot+geom_segment(data = segment.df,
                                 mapping = aes(x=x,xend=xend,y=y,yend=yend),
                                 arrow = arrow(length=unit(0.3, "cm")))
    
  }
  plot.list[[get("index")]]=tmpplot
}
wrap_plots(plot.list,ncol = 12)&theme(plot.margin = unit(c(0,0,0,0),"cm"))
ggsave(paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/all.celltypes.DGEs.", compare[1], ".vs.", compare[2], "-volcano-plot.pdf"), width = 30, height = 18)

## enrichment analysis for each celltypes DGEs
all.celltypes.DGEs.up <- all.celltypes.DGEs %>%
                      filter(p_val_adj < 0.01) %>% filter(avg_log2FC > 0.5) 
eg <- bitr(all.celltypes.DGEs.up$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mmu.eg.db)
eg <- merge(all.celltypes.DGEs.up, eg, by.x = "gene", by.y = "SYMBOL")
gcSample=split(eg$ENTREZID, eg$cluster)
go.all.celltypes.up <- compareCluster(
  gcSample, 
  fun="enrichGO", 
  OrgDb="org.Mmu.eg.db",
  # organism = "mcc",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)
table(go.all.celltypes.up@compareClusterResult$Cluster)
p.up <- dotplot(go.all.celltypes.up, showCategory=dim(go.all.celltypes)[1], font.size = 7) + ggtitle("genes(p.adjust<0.01 & log2FC>0.5) - p.adjust<0.01 & qvalue<0.05 - GSSsize(10, 500)")

all.celltypes.DGEs.do <- all.celltypes.DGEs %>%
                      filter(p_val_adj < 0.01) %>% filter(avg_log2FC < -0.5) 
eg <- bitr(all.celltypes.DGEs.do$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mmu.eg.db)
eg <- merge(all.celltypes.DGEs.do, eg, by.x = "gene", by.y = "SYMBOL")
gcSample=split(eg$ENTREZID, eg$cluster)
go.all.celltypes.do <- compareCluster(
  gcSample, 
  fun="enrichGO", 
  OrgDb="org.Mmu.eg.db",
  # organism = "mcc",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)
table(go.all.celltypes.do@compareClusterResult$Cluster)
p.do <- dotplot(go.all.celltypes.do, showCategory=dim(go.all.celltypes)[1], font.size = 7) + ggtitle("genes(p.adjust<0.01 & log2FC<-0.5) - p.adjust<0.01 & qvalue<0.05 - GSSsize(10, 500)")

ggsave(paste0("./result/3.senescence_related_DEGs/scRNA_DEGs/all.celltypes.DGEs.", compare[1], ".vs.", compare[2], "-GO-dotplot-plot.pdf"), p.up+p.do, width = 10, height = 20)

# if above enrichment codes dont have great performance, than try following codes and adjust parameters
# go.all.celltypes.1 <- compareCluster(
#   gcSample, 
#   fun="enrichGO", 
#   OrgDb="org.Mmu.eg.db",
#   # organism = "mcc",
#   ont = "ALL",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05,
#   minGSSize = 1,
#   maxGSSize = 5000
# )
# table(go.all.celltypes.1@compareClusterResult$Cluster)
# table(go.all.celltypes.1@compareClusterResult$pvalue < 0.05)
# table(go.all.celltypes.1@compareClusterResult$p.adjust < 0.05)

# go.all.celltypes_sim <- clusterProfiler::simplify(go.all.celltypes.1, cutoff=0.7, by="p.adjust", select_fun=min) # remove reaboundant terms
# table(go.all.celltypes_sim@compareClusterResult$Cluster)
# table(go.all.celltypes_sim@compareClusterResult$pvalue < 0.05)
# table(go.all.celltypes_sim@compareClusterResult$p.adjust < 0.05)


