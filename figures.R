###############################################################################################
## R code to reproduce scRNAseq figures (primary and supplementary) from
## "Identification of expanded HER2 specific T cells following vaccination using scRNAseq"
##
## Author: Chaitanya Acharya
## Updated: Jun 16,2020

## Preliminary analysis was run in dataCROSBY.R on Macbook Pro macOS High Sierra v10.13.6
##

## Load the required R libraries
lib.list = c("data.table","plyr","dplyr","ggplot2","GSVA","pheatmap","ggthemes","ggpubr","clustree","devtools","tidyr","gtools","gridExtra","Matrix","rhdf5","car","Seurat","umap","reticulate","RColorBrewer")

for(i in 1:length(lib.list)){
    require(lib.list[i],character.only=T,quietly = T,warn.conflicts=F)
}

## Load utility functions
source("/code/util_funcs.R")

## Define color palette
my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
col_breaks = c( seq(-1,0,length=50),
                seq(0.01,0.8,length=50),
                seq(0.81,1,length=50))

### Path to data input/results
input_path = "../data/"
output_path = "../results/"

### Load previous workspace

load(paste0(input_path,"main_run.rds"))

###############################################
### FIGURES
###############################################

## Original groupings

pdf(paste0(output_path,"original_clusters.pdf"),width=10,height=7)
DimPlot(dat_int,combine = FALSE,label=T)
DimPlot(dat_int,combine = FALSE,label=T,group.by="Rx")
DimPlot(dat_int,combine = FALSE,label=T,group.by="CLONES")

ggbarplot(macro_ave,x="clusters",y="scores",xlab="",ylab="Macrophage Z-score",fill="black")+geom_hline(yintercept=quantile(macro_ave$scores,0.5),col="red",linetype="dashed")+geom_hline(yintercept=0,col="black")

## T cells
FeaturePlot(dat_int, c("Cd3d","Cd3g","Cd3e","Cd2"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## CD8+ T cells
FeaturePlot(dat_int, c("Cd8a","Gzma"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"),label=T)
## CD4+ T cells and Treg
FeaturePlot(dat_int, c("Cd4","Foxp3"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"),label=T)
## NK cells
FeaturePlot(dat_int, c("Klrc1","Klrc3"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Plasma cells
FeaturePlot(dat_int, c("Slamf7","Igkc"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Macrophages iNOS, CD80, MHCII, TLR2,4
FeaturePlot(dat_int, c("Apoe","H2.Aa","C1qa","C1qb","C1qc"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Naive cells
FeaturePlot(dat_int, c("Il7r","Il2rg","Ccr7","Sell"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
dev.off()


## Plotting clonal expansion

pdf(paste0(output_path,"clonal_expansion.pdf"),width=10,height=7)
DimPlot(dat_int,combine = FALSE,label=T)
DimPlot(dat_int,combine = FALSE,label=F,group.by="CLONES",cols=c("lightgrey","steelblue"))
DimPlot(dat_int,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","steelblue"))
DimPlot(dat_int,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","magenta"),split.by="Rx")
DimPlot(dat_int,combine = FALSE,label=F,group.by="new_clusters",cols=cluster_cols)
dev.off()

## New groupings

pdf(paste0(output_path,"new_groupings.pdf"),width=10,height=7)
grid.table(table(dat_int1@meta.data$orig.ident,dat_int1@meta.data$Rx))
clustree(dat_int1,prefix="integrated_snn_res.",layout="sugiyama")
clustree(dat_int1, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
PlotClusterTree(dat_int1)
DimPlot(dat_int1,combine = FALSE,label=T,pt.size=2)
DimPlot(dat_int1,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","steelblue"),pt.size=1)
DimPlot(dat_int1,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","magenta"),split.by="Rx",pt.size=1)
## T cells
FeaturePlot(dat_int1, c("Cd3d","Cd3g","Cd3e","Cd2"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## CD8+ T cells
FeaturePlot(dat_int1, c("Cd8a","Gzma"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"),label=T)
## CD4+ T cells and Treg
FeaturePlot(dat_int1, c("Cd4","Foxp3"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"),label=T)
## NK cells
FeaturePlot(dat_int1, c("Klrc1","Klrc3"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Plasma cells
FeaturePlot(dat_int1, c("Slamf7","Igkc"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Macrophages iNOS, CD80, MHCII, TLR2,4
FeaturePlot(dat_int1, c("Apoe","H2-Aa","C1qa","C1qb","C1qb"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Naive cells
FeaturePlot(dat_int1, c("Il7r","Il2rg","Ccr7","Sell"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
## Checkpoint genes
FeaturePlot(dat_int1, c("Pdcd1","Ctla4","Havcr2","Cd274"), min.cutoff = "q05", max.cutoff = "q95", sort.cell=T,cols = c("grey", "blue"))
dev.off()

## Geneset barplots

pdf(paste0(output_path,"genesets_barplots.pdf"),width=10,height=7)
for(i in 1:nrow(genesets_zscore)){
  print(i)
  print(ggbarplot(zscores,x="Rx",y=names(zscores)[i],ylab=names(zscores)[i],add="mean_se",fill="black")+stat_compare_means(method = "anova",label.y = 0.05)+geom_hline(yintercept=0,col="black"))
  print(ggbarplot(zscores,x="cluster",y=names(zscores)[i],add="mean_se",fill="Rx",ylab=names(zscores)[i],position=position_dodge())+geom_hline(yintercept=0,col="black")+geom_vline(xintercept=seq(1.5,14.5,1)))
}
dev.off()

## New groupings with no Tregs
## Activation and exhausttion signatures

pdf(paste0(output_path,"activation_exh_signatures.pdf"))
FeaturePlot(dat_int2, c("Exhaustion_Sig","Activation_Sig"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"))
FeaturePlot(dat_int2, c("Exhaustion_Sig","Activation_Sig"),min.cutoff = "q05", max.cutoff = "q95",sort.cell=T,cols=c("grey", "blue"),combine=T,split.by="Rx")
p1=ggbarplot(dat_int2@meta.data,x="Rx",y="Exhaustion_Sig",xlab="Treatment",ylab="Exhaustion_Sig",fill="black",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.02) + geom_hline(yintercept=0,col="black")
p2=ggbarplot(dat_int2@meta.data,x="Rx",y="Activation_Sig",xlab="Treatment",ylab="Activation_Sig",fill="black",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.1) + geom_hline(yintercept=0,col="black")
CombinePlots(list(p1,p2))
p3=ggbarplot(dat_int2@meta.data,x="seurat_clusters",y="Exhaustion_Sig",xlab="Treatment",ylab="Exhaustion_Sig",fill="black",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.02) + geom_hline(yintercept=0,col="black")
p4=ggbarplot(dat_int2@meta.data,x="seurat_clusters",y="Activation_Sig",xlab="Treatment",ylab="Activation_Sig",fill="black",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.05) + geom_hline(yintercept=0,col="black")
CombinePlots(list(p3,p4))
dev.off()



