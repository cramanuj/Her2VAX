require(Seurat); require(data.table); require(plyr); require(dplyr); require(ggplot2)
require(GSVA); require(pheatmap); require(ggthemes); require(ggpubr); require(clustree)
require(psych); require(RColorBrewer); require(limma); require(biomaRt); library(tidyr)
require(matrixStats)

source("~/Research/scripts/r_scripts/useful_functions.R")
source("~/Research/scripts/r_scripts/plotfns.R")
source("~/Research/Zack/read_hi5file.R")

file1 = read_h5("~/Research/Zack/ANALYSIS/GEX/GEX_401AB_result/filtered_feature_bc_matrix.h5")
file2 = read_h5("~/Research/Zack/ANALYSIS/GEX/GEX_403A_result/filtered_feature_bc_matrix.h5")
file3 = read_h5("~/Research/Zack/ANALYSIS/GEX/GEX_403B_result/filtered_feature_bc_matrix.h5")
file4 = read_h5("~/Research/Zack/ANALYSIS/GEX/GEX_464_result/filtered_feature_bc_matrix.h5")
file5 = read_h5("~/Research/Zack/ANALYSIS/GEX/GEX_471_result/filtered_feature_bc_matrix.h5")
file6 = read_h5("~/Research/Zack/NEW_DATA/ANALYSIS_RESULTS/GEX-663/filtered_feature_bc_matrix.h5")
file7 = read_h5("~/Research/Zack/NEW_DATA/ANALYSIS_RESULTS/GEX-668/filtered_feature_bc_matrix.h5")
file8 = read_h5("~/Research/Zack/NEW_DATA/ANALYSIS_RESULTS/GEX-669/filtered_feature_bc_matrix.h5")
file9 = read_h5("~/Research/Zack/NEW_DATA/ANALYSIS_RESULTS/GEX-672/filtered_feature_bc_matrix.h5")

cc_genes = scan("~/Research/Josh/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt",what="")
s_genes = convertHumanGeneList(cc_genes[1:43])
g2m_genes = convertHumanGeneList(cc_genes[44:97])

data_list = list(file1$mat,file2$mat, file3$mat, file4$mat,file5$mat,file6$mat,file7$mat,file8$mat,file9$mat)
files = c("GEX_401AB","GEX_403A","GEX_403B","GEX_464","GEX_471","GEX-663","GEX-668","GEX-669","GEX-672")
seu_list = gene_list = vector(mode="list",length=length(files))
for(i in 1:length(seu_list)){
  print(i)
  dat = data.frame(data_list[[i]],check.names=F)
  dat$Genes = as.character(file1$feat$name)
  dat = setDT(dat)[, lapply(.SD, sum), by = Genes]
  dat = as.data.frame(dat)
  rownames(dat) = dat[,1]
  dat = dat[,-1]
  names(dat)=paste(files[i],colnames(dat),sep="-")
  dat = dat[-c(grep("Rpl",rownames(dat)),grep("Rps",rownames(dat)),grep("^mt-",rownames(dat))),]
  seu = CreateSeuratObject(raw.data=dat,min.cells=5)
  seu@meta.data$group = files[i]
  # seu = FilterCells(seu, subset.names = c("nGene","nUMI"), low.thresholds = c(round(quantile(seu@meta.data$nGene,0.1)),-Inf), high.thresholds = c(round(quantile(seu@meta.data$nGene,0.95)),round(quantile(seu@meta.data$nGene,0.99))))
  seu = FilterCells(seu, subset.names = "nGene", low.thresholds = round(quantile(seu@meta.data$nGene,0.1)), high.thresholds = Inf)
  seu = NormalizeData(seu,display.progress = F)
  seu = CellCycleScoring(seu, s.genes = s_genes, g2m.genes = g2m_genes, set.ident = F)
  seu@meta.data$CC.Difference = seu@meta.data$S.Score - seu@meta.data$G2M.Score
  seu = ScaleData(seu, display.progress = F,vars.to.regress = c("nUMI"))
  seu = FindVariableGenes(seu,display.progress = F, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR)
  genes_use = head(rownames(seu@hvg.info), 2000)
  seu_list[[i]] = seu
  gene_list[[i]] = genes_use
}
rm(file1,file2,file3,file4,file5,file6,file7,file8,file9,file10)
llply(1:length(seu_list),function(i) dim(seu_list[[i]]@data))

#######################################################################################
## Gene selection for CCA
## We find genes that are highly variable in at least two datasets
#######################################################################################

genes_use <- names(which(table(unlist(gene_list)) > 1))
for (i in 1:length(seu_list)) {
  genes_use <- genes_use[genes_use %in% rownames(seu_list[[i]]@scale.data)]
}

#######################################################################################
### Run multi-set Canonical Correlation Analysis (CCA) on the common genes across all the genesets
### Calculate the ratio of total variance explained by PPCA vs total variance explained by CCA,
###       and filter cells based on these values
### CCA Align the matrices and compute the alignment metric score
### Map a t-SNE plot and find clusters
#######################################################################################

pdf("CCAnew_plots.pdf",width=10,height=7)
cca_out = RunMultiCCA(seu_list,genes.use=genes_use,num.ccs = 10)
DimPlot(object = cca_out, reduction.use = "cca", group.by = "group", pt.size = 1, do.return = F)
VlnPlot(object = cca_out, features.plot = "CC1", group.by = "group", do.return = F)
MetageneBicorPlot(cca_out, grouping.var = "group", dims.eval = 1:20, display.progress = FALSE)
DimHeatmap(object = cca_out, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
cca_out = CalcVarExpRatio(cca_out,reduction.type = "pca", grouping.var = "group", dims.use = 1:20)
cca_out = SubsetData(cca_out, subset.name = "var.ratio.pca",accept.low = 0.40)
cca_out = AlignSubspace(cca_out, reduction.type = "cca", grouping.var = "group", dims.align = 1:10,num.genes=50)
metric = rep(0,10)
for(i in 1:10){
  metric[i]=CalcAlignmentMetric(cca_out,reduction.use = "cca.aligned",dims.use = 1:i, grouping.var =  "group")
  cat("Number of PCs: ",i," Alignment metric: ",metric[i],"\n")
}
# cca_out@meta.data$Rx = ifelse(cca_out@meta.data$group %in% c("GEX_401AB","GEX_464"),"Control","Vaccine+aPD1")
# cca_out@meta.data$Rx = ifelse(cca_out@meta.data$group %in% c("GEX_401AB","GEX_464"),"Control",ifelse(cca_out@meta.data$group %in% c("GEX_403A","GEX_403B"),"aPD1","Vaccine_aPD1"))
max_pc = 8
cca_out = FindClusters(cca_out, reduction.type = "cca.aligned", resolution = c(0.4,0.6,0.8,1), dims.use = 1:max_pc,print.output = F)
clustree(cca_out,prefix="res.",layout="sugiyama")
clustree(cca_out, prefix = "res.", node_colour = "sc3_stability")
cca_out = SetAllIdent(cca_out,id="res.0.6")
cca_out = BuildClusterTree(cca_out)
PlotClusterTree(cca_out)

cca_out = RunTSNE(cca_out, reduction.use = "cca.aligned", dims.use = 1:max_pc, do.fast = T,check_duplicates=F)
tsne_plot = TSNEPlot(object = cca_out, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5)
TSNEPlot(object = cca_out, do.label = F,group.by="group")
TSNEPlot(object = cca_out, do.label = F,group.by="Rx")

# cca_out = RunUMAP(cca_out,reduction.use = "cca.aligned",dims.use = 1:max_pc)
cca_out = RunUMAP(cca_out,reduction.use = "cca.aligned",dims.use = 1:5,n_neighbors = 50L,min_dist=0.5,metric="cosine",umap.method = "umap-learn")
umap_plot = DimPlot(object = cca_out, reduction.use = 'umap',do.label=T,label.size=5)
DimPlot(object = cca_out, reduction.use = 'umap',group.by="group")
DimPlot(object = cca_out, reduction.use = 'umap',group.by="Rx")
ggarrange(tsne_plot,umap_plot,nrow=1,ncol=2)
dev.off()

##############################
# set.seed(123)
# com_graph = MUDAN::getComMembership(cca_out@dr$cca.aligned@cell.embeddings, k=30, method=igraph::cluster_louvain,verbose=T)
# plotEmbedding(emb=cca_out@dr$umap@cell.embeddings, groups=com_graph, main='Graph-based Community Detection', xlab="Component 1", ylab="Component 2", mark.clusters=TRUE, alpha=0.8, mark.cluster.cex=2,verbose=T,show.legend=F)
# plotEmbedding(emb=cca_out@dr$tsne@cell.embeddings, groups=com_graph, main='Graph-based Community Detection', xlab="Component 1", ylab="Component 2", mark.clusters=TRUE, alpha=0.8, mark.cluster.cex=2,verbose=T,show.legend=F)
##############################


# FeaturePlot(cca_out, features.plot = c("Ptprc","Cd3d","Cd4","Cd8a"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
# FeaturePlot(cca_out, features.plot = c("Gzmb","Ifng","Prf1","Gzma","Pdcd1","Ctla4"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

pdf("feature_plots.pdf",width=10,height=7)
DimPlot(object = cca_out, reduction.use = 'umap',do.label=T,label.size=5,no.legend=F)
## T cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd3d","Cd3g","Cd3e","Cd2"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## CD8+ T cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd8a","Gzma"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## CD4+ T cells and Treg
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd4","Foxp3"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## NK cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Klrc1","Klrc3"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## B cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd79a"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## Plasma cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Slamf7","Igkc"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## Dendritic cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Flt3"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
### Macrophages iNOS, CD80, MHCII, TLR2,4
FeaturePlot(cca_out, no.legend=F,features.plot = c("Apoe","H2-Aa","C1qa","C1qb","C1qb"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

FeaturePlot(cca_out, no.legend=F,features.plot = c("Pdcd1","Ctla4","Havcr2","Cd274"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
dev.off()

##########################################################
## Cluster-specific markers
##########################################################
markers = FindAllMarkers(object = cca_out,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = cca_out, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F,cex.row = 5,group.cex = 7,col.low="blue",col.mid="white",col.high="red")

##########################################################
## Macrophage activity
##########################################################
macro_geneset = c("Apoe","H2-Aa","C1qa","C1qb","C1qc")
macro = as.matrix(cca_out@data[macro_geneset,])
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
macro_score = as.numeric(scale(apply(macro,2,gm_mean)))
cca_out@meta.data$macrophages = macro_score

my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
col_breaks = c( seq(-1,0,length=50),
                seq(0.01,0.8,length=50),
                seq(0.81,1,length=50))
# ggplot(ES4_mean, aes(cluster, Genesets)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradientn(colours = my_palette,breaks=col_breaks)
# #
# ES5_mean = acharya_melt %>% dplyr::group_by(cluster,Genesets,group) %>% summarize_all(list(winsor.mean))
# ggplot(ES5_mean, aes(cluster, Genesets)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradientn(colours = my_palette,breaks=col_breaks) +
#       theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + facet_grid(~group)
#
# ES6_mean = aggregate(value ~ Genesets + group,ES5_mean,FUN="mean")
# ES6_mean = ES6_mean[order(ES6_mean$value,decreasing=T),]
# ES6_mean$Genesets = factor(ES6_mean$Genesets,levels=unique(split(ES6_mean,ES6_mean$group)$Control$Genesets))
# ggplot(ES6_mean,aes(x=Genesets,y=value,fill=group)) + geom_bar(stat="identity",width=0.6)+coord_flip()+theme_classic()+scale_fill_manual(values=genespring.colors(3))
#
# col_breaks = c( seq(-1,0,length=50),               # for red
#                 seq(0.01,0.8,length=50),           # for yellow
#                 seq(0.81,1,length=50))             # for green
# my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
# tmp = cca_out@meta.data[,c("group","Phase","res.0.6")]
# colnames(tmp)[3]=c("cluster")
# tmp$group=as.factor(tmp$group); tmp$Phase=as.factor(tmp$Phase); tmp$cluster=as.factor(tmp$cluster)
# # tmp$cluster = factor(tmp$cluster,levels=c(0:nlevels(tmp$cluster)))
# tmp = tmp[order(tmp$cluster),]
# tmp_dat = Matrix::t(scale(Matrix::t(acharya_gsva)))
# tmp_dat = tmp_dat[,match(rownames(tmp),colnames(tmp_dat))]
# mat_col = list(group=genespring.colors(5),Phase=c("red","blue","yellow"),cluster=matlab.colors(nlevels(tmp$cluster)))
# names(mat_col[[1]])=levels(tmp$group)
# names(mat_col[[2]])=levels(tmp$Phase)
# names(mat_col[[3]])=levels(tmp$cluster)
# pheatmap(tmp_dat,color=my_palette,breaks=col_breaks,cluster_row=T,cluster_col=F,annotation_col=tmp,annotation_colors =mat_col,show_colnames = F,cutree_rows=4,treeheight_row = 0,gaps_col=cumsum(table(tmp$cluster)))
#
# cca_out@meta.data = cbind(cca_out@meta.data,data.frame(t(Matrix::t(scale(Matrix::t(azizi_gsva)))),check.names=F),data.frame(t(Matrix::t(scale(Matrix::t(acharya_gsva)))),check.names=F))
# FeaturePlot(cca_out, no.legend=F,features.plot = c("TCELL_CYTOLYTIC_ACT","Treg","TCELL_ACTIVATION","CP_INHIBITORY"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
# FeaturePlot(cca_out, no.legend=F,features.plot = c("Type I Interferon response","Type II Interferon Response","Pro-inflammatory","Anti-inflammatory"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
# FeaturePlot(cca_out, no.legend=F,features.plot = c("M1 Macrophage Polarization","M2 Macrophage Polarization"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

#### Subsetting on macrophage activity
cca_out2 = SubsetData(cca_out,ident.remove = c(3,7,8,9,10,11,12),subset.raw=T)
cca_out2@meta.data = cca_out2@meta.data[,-c(11:14)]

cca_out2 = FindClusters(cca_out2, reduction.type = "cca.aligned", resolution = c(0.4,0.6,0.8,1), dims.use = 1:max_pc,print.output = F,force.recalc = TRUE)
clustree(cca_out2,prefix="res.",layout="sugiyama")
clustree(cca_out2, prefix = "res.", node_colour = "sc3_stability")
cca_out2 = SetAllIdent(cca_out2,id="res.0.8")
cca_out2 = BuildClusterTree(cca_out2)
PlotClusterTree(cca_out2)

cca_out2 = RunTSNE(cca_out2, reduction.use = "cca.aligned", dims.use = 1:max_pc, do.fast = T,check_duplicates=F)
tsne_plot = TSNEPlot(object = cca_out2, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5)
TSNEPlot(object = cca_out2, do.label = F,group.by="group")
TSNEPlot(object = cca_out2, do.label = F,group.by="Rx")

# cca_out = RunUMAP(cca_out,reduction.use = "cca.aligned",dims.use = 1:max_pc)
cca_out2 = RunUMAP(cca_out2,reduction.use = "cca.aligned",dims.use = 1:5,n_neighbors = 50L,min_dist=0.5,metric="cosine",umap.method = "umap-learn")
umap_plot = DimPlot(object = cca_out2, reduction.use = 'umap',do.label=T,label.size=5)
DimPlot(object = cca_out, reduction.use = 'umap',group.by="group")
DimPlot(object = cca_out, reduction.use = 'umap',group.by="Rx")
ggarrange(tsne_plot,umap_plot,nrow=1,ncol=2)
dev.off()

###########################################################
## ssGSEA using Acharya genesets
###########################################################
acharya_gs = read.gmt.file("~/Research/pathways/mouse_pathways/acharya_genesets_mouse.gmt")
acharya_gs$geneset.descriptions = acharya_gs$geneset.descriptions[c(1:16,37:38)]
acharya_gs$geneset.names = acharya_gs$geneset.names[c(1:16,37:38)]
acharya_gs$genesets = acharya_gs$genesets[c(1:16,37:38)]
acharya_gsva = gsva(as.matrix(cca_out2@data),acharya_gs$genesets,method="ssgsea",ssgsea.norm = T)
# es = apply(es, 2, function(x, es) x / (range(es)[2] - range(es)[1]), es)
# acharya_out = acharya_gsva[,match(rownames(cca_out2@meta.data),colnames(acharya_gsva))]
# mmin = rowMins(acharya_out)
# mmax = rowMaxs(acharya_out)
# scores = acharya_out/(mmax - mmin)
# rownames(scores) = rownames(a)
# acharya_out = data.frame(t(scores),check.names=F)
acharya_out = data.frame(t(Matrix::t(scale(Matrix::t(acharya_out)))),check.names=F)
acharya_out$cluster = cca_out2@meta.data$res.0.6
acharya_out$group = as.factor(cca_out2@meta.data$Rx)
acharya_melt = melt(acharya_out)

##########################################################
## CD8+ T cells
##########################################################
cca_out2@meta.data$CD8 = cca_out2@scale.data[grep("Cd8a",rownames(cca_out2@scale.data)),]
cd8 = cca_out2@meta.data[,c("res.0.6","CD8")]
colnames(cd8)[1]=c("cluster")
cd8$cluster = factor(cd8$cluster,levels=0:(length(table(cd8$cluster))-1))
ggbarplot(cd8,x="cluster",y="CD8",add="mean_se",fill="grey") + stat_compare_means(method = "anova",label.y = 0.5)

# cca_out@meta.data$macrophages = acharya_out$MACROPHAGE_ACTIVITY
# cca_out@meta.data$dendritic_cells = acharya_out$DENDRITIC_CELLS
# FeaturePlot(cca_out, no.legend=F,features.plot = c("macrophages","dendritic_cells"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

################################################
## TCR Overlay
################################################
tmp_barcode = rownames(cca_out2@meta.data)
tmp_barcode = laply(1:length(tmp_barcode),function(i) paste(strsplit(tmp_barcode[i],split="-")[[1]][-1],collapse="-"))
cca_out2@meta.data$new_barcode = tmp_barcode

TCRfiles = gsub("GEX","TCR",files)
TCRfiles[4] = paste0(TCRfiles[4],"A"); TCRfiles[5] = paste0(TCRfiles[5],"A");
tcr_dat=list()
for(i in 1:length(TCRfiles)){
  print(i)
  tcr_dat[[i]] = read.csv(paste0("/Users/ca31/Research/Zack/ANALYSIS/TCR/",TCRfiles[i],"_result/filtered_contig_annotations.csv"),header=T)
  tcr_dat[[i]]$group = TCRfiles[i]
}
tcr_dat = do.call("rbind",tcr_dat)
cca_out@meta.data$CLONES = (cca_out@meta.data$new_barcode %in% tcr_dat$barcode) + 0
FeaturePlot(cca_out, no.legend=F,features.plot = c("CLONES"), reduction.use="umap",nCol=1,cols.use = c("lightgrey", "steelblue"), pt.size = 1)

# cca_out1 = SubsetData(cca_out,ident.remove = macro_mean[which(macro_mean$value>quantile(acharya_out$MACROPHAGE_ACTIVITY,0.5)),"cluster"])
cca_out1 = SubsetData(cca_out,cells.use = rownames(cca_out@meta.data[cca_out@meta.data$macrophages<quantile(cca_out@meta.data$macrophages,0.5),]))
tcr_dat_new = tcr_dat[!duplicated(tcr_dat$barcode),]
tcr_dat_new$new_clonotype = paste(tcr_dat_new$raw_clonotype_id,tcr_dat_new$group,sep="_")
tcr_merged = base::merge(cca_out1@meta.data,tcr_dat_new,by.x=names(cca_out1@meta.data)[17],by.y=names(tcr_dat_new)[1],all=T)
tcr_merged = tcr_merged[match(cca_out1@meta.data$new_barcode,tcr_merged$new_barcode),]
dim(tcr_merged)
cca_out1@meta.data = cbind(cca_out1@meta.data,tcr_merged[,18:ncol(tcr_merged)])
#####
# clone_table = data.frame(table(tcr_dat_new$new_clonotype,tcr_dat_new$group))
# clone_tab_melt = melt(clone_table)
# clone_tab_melt = clone_tab_melt[clone_tab_melt$value!=0,]
# colnames(clone_tab_melt)[1:2]=c("clonotypes","group")
# clone_tab_melt1 = clone_tab_melt[clone_tab_melt$group==names(table(clone_tab_melt$group))[1],]
# clone_tab_melt1 = clone_tab_melt1[order(clone_tab_melt1$value,decreasing=T),]
# clone_tab_melt1$clonotypes = factor(clone_tab_melt1$clonotypes,levels=clone_tab_melt1$clonotypes)
# ggbarplot(clone_tab_melt1,x="clonotypes",y="value",fill="steelblue",xlim=c(1:50),ylim=c(1:max(clone_tab_melt1$value)),ylab="Frequency of clonotypes")+theme(axis.title.x=element_blank())
# clone_tab_melt2 = clone_tab_melt[clone_tab_melt$group==names(table(clone_tab_melt$group))[2],]
# clone_tab_melt2 = clone_tab_melt2[order(clone_tab_melt2$value,decreasing=T),]
# clone_tab_melt2$clonotypes = factor(clone_tab_melt2$clonotypes,levels=clone_tab_melt2$clonotypes)
# ggbarplot(clone_tab_melt2,x="clonotypes",y="value",fill="steelblue",xlim=c(1:50),ylim=c(1:max(clone_tab_melt2$value)),ylab="Frequency of clonotypes")+theme(axis.title.x=element_blank())
# clone_tab_melt3 = clone_tab_melt[clone_tab_melt$group==names(table(clone_tab_melt$group))[3],]
# clone_tab_melt3 = clone_tab_melt3[order(clone_tab_melt3$value,decreasing=T),]
# clone_tab_melt3$clonotypes = factor(clone_tab_melt3$clonotypes,levels=clone_tab_melt3$clonotypes)
# ggbarplot(clone_tab_melt3,x="clonotypes",y="value",fill="steelblue",xlim=c(1:50),ylim=c(1:max(clone_tab_melt3$value)),ylab="Frequency of clonotypes")+theme(axis.title.x=element_blank())
# clone_tab_melt4 = clone_tab_melt[clone_tab_melt$group==names(table(clone_tab_melt$group))[4],]
# clone_tab_melt4 = clone_tab_melt4[order(clone_tab_melt4$value,decreasing=T),]
# clone_tab_melt4$clonotypes = factor(clone_tab_melt4$clonotypes,levels=clone_tab_melt4$clonotypes)
# ggbarplot(clone_tab_melt4,x="clonotypes",y="value",fill="steelblue",xlim=c(1:50),ylim=c(1:max(clone_tab_melt4$value)),ylab="Frequency of clonotypes")+theme(axis.title.x=element_blank())
# clone_tab_melt5 = clone_tab_melt[clone_tab_melt$group==names(table(clone_tab_melt$group))[5],]
# clone_tab_melt5 = clone_tab_melt5[order(clone_tab_melt5$value,decreasing=T),]
# clone_tab_melt5$clonotypes = factor(clone_tab_melt5$clonotypes,levels=clone_tab_melt5$clonotypes)
# ggbarplot(clone_tab_melt5,x="clonotypes",y="value",fill="steelblue",xlim=c(1:50),ylim=c(1:max(clone_tab_melt5$value)),ylab="Frequency of clonotypes")+theme(axis.title.x=element_blank())
# ###
clone_table = table(tcr_dat_new$new_clonotype,tcr_dat_new$group)
colSums(clone_table>5)
clones = rownames(which(clone_table>5,arr.ind=T))
clones = clones[-grep("None",clones)]
clones = data.frame("clonotype"=laply(1:length(clones),function(i) paste(strsplit(clones[i],split="_")[[1]][1],sep="_")),"group"=laply(1:length(clones),function(i) paste(strsplit(clones[i],split="_")[[1]][3],sep="_")))
clones$clonotype = factor(clones$clonotype,levels=paste0("clonotype",1:18))

cca_out3 = cca_out
cca_out3@meta.data = cca_out@meta.data[cca_out@meta.data$new_clonotype %in% clones,]
cca_out3 = SubsetData(cca_out3,cells.use = rownames(cca_out3@meta.data),subset.raw=T)
cca_out3@meta.data$new_clonotype = as.factor(as.character( cca_out3@meta.data$new_clonotype))
cca_out3@meta.data$raw_clonotype_id = as.factor(as.character( cca_out3@meta.data$raw_clonotype_id))
cca_out3@meta.data$raw_clonotype_id = factor(cca_out3@meta.data$raw_clonotype_id,levels=paste0("clonotype",1:18))
DimPlot(object = cca_out3, reduction.use = 'umap',group.by="raw_clonotype_id",do.label=F,pt.size=2) + facet_grid(~cca_out3@meta.data$Rx)

#############################################################################################
cca_temp = cca_out1
cca_temp@meta.data = cca_temp@meta.data[cca_temp@meta.data$CLONES==1,]
cca_temp = SubsetData(cca_temp,cells.use=rownames(cca_temp@meta.data),subset.raw=T)
cca_out2 = cca_temp
# cca_out2 = SubsetData(cca_out,ident.use = c(3),subset.raw=T)
cca_out2@meta.data = cca_out2@meta.data[,-c(10:12)]
cca_out2 = FindClusters(cca_out2, reduction.type = "cca.aligned", resolution = c(0.4,0.6,0.8,1), dims.use = 1:10,print.output = F,force.recalc=T)
clustree(cca_out2,prefix="res.",layout="sugiyama")
clustree(cca_out2, prefix = "res.", node_colour = "sc3_stability")
cca_out2 = SetAllIdent(cca_out2,id="res.0.4")
cca_out2 = BuildClusterTree(cca_out2)
PlotClusterTree(cca_out2)

cca_out2 = RunUMAP(cca_out2,reduction.use = "cca.aligned",dims.use = 1:10)
cca_out2@meta.data$Rx = factor(cca_out2@meta.data$Rx,levels=c("Control","aPD1","Vaccine_aPD1"))
DimPlot(object = cca_out2, reduction.use = 'umap',do.label=T,facet.by="group",pt.size=3,label.size=5)
DimPlot(object = cca_out2, reduction.use = 'umap',do.label=F,pt.size=3,label.size=5) + facet_grid(~cca_out2@meta.data$Rx)

FeaturePlot(cca_out2, features.plot = c("Ptprc","Cd3d","Cd4","Cd8a"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1) + facet_grid(~cca_out2@meta.data$Rx)
FeaturePlot(cca_out2, features.plot = c("Gzmb","Ifng","Prf1","Gzma","Pdcd1"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
# FeaturePlot(cca_out2, features.plot = c("TCELL_CYTOLYTIC_ACT","Treg","TCELL_ACTIVATION","CP_INHIBITORY"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

markers2 = FindAllMarkers(object = cca_out2,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
write.table(markers2,"CD8_clusters_withTCR.txt",sep="\t",col.names=NA,quote=F)
top10new = markers2 %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = cca_out2, genes.use = top10new$gene, slim.col.label = TRUE, remove.key = F,cex.row = 5,group.cex = 7,col.low="blue",col.mid="white",col.high="red")

#################################################
## ssGSEA using genesets from Azizi et al
#################################################

azizi_gs = read.gmt.file("~/Research/pathways/mouse_pathways/azizi_genesets_mouse.gmt")
azizi_gsva = gsva(as.matrix(cca_out2@data),azizi_gs$genesets,method="ssgsea")
azizi_out = data.frame(t(Matrix::t(scale(Matrix::t(azizi_gsva)))),check.names=F)
azizi_out$cluster = cca_out2@meta.data$res.0.4
azizi_out$Rx = as.factor(cca_out2@meta.data$Rx)
azizi_melt = melt(azizi_out)
names(azizi_melt)[3] = c("Genesets")
ES1_mean = azizi_melt[,-2] %>% dplyr::group_by(Genesets,cluster) %>% summarize_all(list(winsor.mean))
my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
col_breaks = c( seq(-1,0,length=50),               # for red
                seq(0.01,0.8,length=50),           # for yellow
                seq(0.81,1,length=50))
tmp = cca_out2@meta.data[,c("res.0.6","Rx")]
colnames(tmp)[1:2]=c("cluster","Rx")
tmp = tmp[order(tmp$Rx),]
tmp$cluster = as.factor(tmp$cluster)
tmp_dat = Matrix::t(scale(Matrix::t(azizi_gsva)))
tmp_dat = tmp_dat[,match(rownames(tmp),colnames(tmp_dat))]
mat_col = list(cluster=matlab.colors(nlevels(tmp$cluster)),Rx=genespring.colors(3))
names(mat_col[[2]])=levels(tmp$Rx)
names(mat_col[[1]])=levels(tmp$cluster)
pheatmap(tmp_dat,color=my_palette,breaks=col_breaks,cluster_row=T,cluster_col=F,annotation_col=tmp,annotation_colors =mat_col,show_colnames = F,cutree_rows=3,treeheight_row = 0,gaps_col=cumsum(table(tmp$Rx)))

ggplot(ES1_mean, aes(cluster, Genesets)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradientn(colours = my_palette,breaks=col_breaks)+theme(legend.position = "none")+ylab("")

ES2_mean = azizi_melt %>% dplyr::group_by(cluster,Genesets,Rx) %>% summarize_all(list(winsor.mean))
ggplot(ES2_mean, aes(cluster, Genesets)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradientn(colours = my_palette,breaks=col_breaks) +
      theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + facet_grid(~Rx)+theme(legend.position = "none")+ylab("")

# ES3_mean = ES2_mean %>% arrange(group,value,.by_group=T)
# ES3_mean$Genesets = factor(ES3_mean$Genesets,levels=unique(ES3_mean$Genesets))
# ggbarplot(ES3_mean,x="Genesets",y="value",add="mean",fill="group",position=position_dodge(),palette=genespring.colors(2))
# ggplot(ES3_mean,aes(x=Genesets,y=value))+geom_bar(stat='identity', aes(fill=Genesets),width=.5)+coord_flip() + theme_classic()

ES3_mean = aggregate(value ~ Genesets + Rx,ES2_mean,FUN="mean")
ES3_mean = ES3_mean[order(ES3_mean$value,decreasing=T),]
ES3_mean$Genesets = factor(ES3_mean$Genesets,levels=unique(split(ES3_mean,ES3_mean$Rx)$Control$Genesets))
ggplot(ES3_mean,aes(x=Genesets,y=value,fill=Rx)) + geom_bar(stat="identity",width=0.6)+coord_flip()+theme_classic()+scale_fill_manual(values=genespring.colors(3))

ES3_mean = aggregate(value ~ Genesets + Rx,ES2_mean,FUN="mean")
ES3_mean = ES3_mean[order(ES3_mean$value,decreasing=T),]
ES4_mean = ES3_mean %>% arrange(Genesets) %>% group_split(Rx) %>% data.frame()
ES4_mean = ES4_mean[,-c(2,4,5,7,8)]
colnames(ES4_mean)[2:4]=c("Control","aPD1","Vaccine_aPD1")
ggbarplot(melt(ES4_mean),x="Genesets",y="value",add="mean_se",fill="variable",position=position_dodge()) + stat_compare_means(method = "anova",label.y = 0.2)
tmp = cca_out2@meta.data[,c("group","Phase","res.0.4")]
colnames(tmp)[3]=c("cluster")
tmp$group=as.factor(tmp$group); tmp$Phase=as.factor(tmp$Phase); tmp$cluster=as.factor(tmp$cluster)
tmp = tmp[order(tmp$cluster),]
tmp_dat = Matrix::t(scale(Matrix::t(azizi_gsva)))
tmp_dat = tmp_dat[,match(rownames(tmp),colnames(tmp_dat))]
mat_col = list(group=genespring.colors(5),Phase=c("red","blue","yellow"),cluster=matlab.colors(nlevels(tmp$cluster)))
names(mat_col[[1]])=levels(tmp$group)
names(mat_col[[2]])=levels(tmp$Phase)
names(mat_col[[3]])=levels(tmp$cluster)
pheatmap(tmp_dat,color=my_palette,breaks=col_breaks,cluster_row=T,cluster_col=F,annotation_col=tmp,annotation_colors =mat_col,show_colnames = F,cutree_rows=4,treeheight_row = 0,gaps_col=cumsum(table(tmp$Rx)))

###
##
IFNG_sig = convertHumanGeneList(c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG"))
IMM_sig = convertHumanGeneList(c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB"))
Tcell_inflamedGEP_sig = convertHumanGeneList(c("TIGIT","CD27","CD8A","PDCD1LG2","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E"))
# cca_out2 = AddModuleScore(cca_out2,genes.list = list(IFNG_sig,IMM_sig,Tcell_inflamedGEP_sig),enrich.name=c("IFNG_sig","IMM_sig","TCELL_INF_GEP"))
# colnames(cca_out2@meta.data)[c(16:18)] = c("IFNG_sig","IMM_sig","TCELL_INF_GEP")
imm_gsva = gsva(as.matrix(cca_out2@data),list(IFNG_sig,IMM_sig,Tcell_inflamedGEP_sig),method="ssgsea")
imm_gsva = data.frame(t(Matrix::t(scale(Matrix::t(imm_gsva)))),check.names=F)
colnames(imm_gsva)[1:3] = c("IFNG_sig","IMM_sig","TCELL_INF_GEP")
cca_out2@meta.data = cbind(cca_out2@meta.data,imm_gsva)
FeaturePlot(cca_out2, features.plot = c("IFNG_sig","IMM_sig","TCELL_INF_GEP"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
cca_out2@meta.data$Samples = rownames(cca_out2@meta.data)
norm_int_sig_data = cca_out2@meta.data[,c("Rx","res.0.6","IFNG_sig","IMM_sig","TCELL_INF_GEP")]
ctrl = melt(norm_int_sig_data[norm_int_sig_data$Rx=="Control",])
pd1 = melt(norm_int_sig_data[norm_int_sig_data$Rx=="aPD1",])
vac = melt(norm_int_sig_data[norm_int_sig_data$Rx=="Vaccine_aPD1",])
ggbarplot(melt(cca_out2@meta.data[,c("Rx","res.0.6","IFNG_sig","IMM_sig","TCELL_INF_GEP")]),x="res.0.6",y="value",fill="variable",add="mean_se",position=position_dodge(0.8),facet.by="Rx")
col_breaks = c( seq(-1,0,length=50),               # for red
                seq(0.01,0.8,length=50),           # for yellow
                seq(0.81,1,length=50))             # for green
my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
tmp = norm_int_sig_data[,c(1,2)]
colnames(tmp)[1:2]=c("Rx","cluster")
tmp = tmp[order(tmp$Rx),]
tmp$cluster = as.factor(tmp$cluster)
tmp_dat= t(norm_int_sig_data[,-c(1,2)])
tmp_dat = tmp_dat[,match(rownames(tmp),colnames(tmp_dat))]
mat_col = list(Rx=genespring.colors(3),cluster=matlab.colors(nlevels(tmp$cluster)))
names(mat_col[[1]])=levels(tmp$Rx)
names(mat_col[[2]])=levels(tmp$cluster)
pheatmap(tmp_dat,color=my_palette,breaks=col_breaks,cluster_row=T,cluster_col=F,annotation_col=tmp,annotation_colors =mat_col,show_colnames = F,cutree_rows=3,treeheight_row = 0,gaps_col=cumsum(table(tmp$Rx)))

# cca_out2 = AddModuleScore(cca_out2,genes.list = azizi_gs$genesets[c(3,5,16)],enrich.name=names(azizi_gs$genesets[c(3,5,16)]))
cca_out2@meta.data = cbind(cca_out2@meta.data,azizi_out[,c(3,5,16)])

norm_int_sig_data = cca_out2@meta.data[,c("Rx","res.0.6","IFNG_sig","IMM_sig","TCELL_INF_GEP","Anti-inflammatory","Pro-inflammatory","Type II Interferon Response")]
colnames(norm_int_sig_data)[6:8] = c("Anti-inflammatory","Pro-inflammatory","Type II Interferon Response")
tmp = norm_int_sig_data[,c(1,2)]
colnames(tmp)[1:2]=c("Rx","cluster")
tmp = tmp[order(tmp$Rx),]
tmp$cluster = as.factor(tmp$cluster)
tmp_dat= t(norm_int_sig_data[,-c(1,2)])
tmp_dat = tmp_dat[,match(rownames(tmp),colnames(tmp_dat))]
mat_col = list(Rx=genespring.colors(3),cluster=matlab.colors(nlevels(tmp$cluster)))
names(mat_col[[1]])=levels(tmp$Rx)
names(mat_col[[2]])=levels(tmp$cluster)
pheatmap(tmp_dat,color=my_palette,breaks=col_breaks,cluster_row=T,cluster_col=F,annotation_col=tmp,annotation_colors =mat_col,show_colnames = F,cutree_rows=4,treeheight_row = 0,gaps_col=cumsum(table(tmp$Rx)))
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Type II Interferon Response")]),x="Rx",y="value",xlab="Treatment",ylab="Type II Interferon Response",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.2)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Pro-inflammatory")]),x="Rx",y="value",xlab="Treatment",ylab="Pro-inflammatory2",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.2)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Anti-inflammatory")]),x="Rx",y="value",xlab="Treatment",ylab="Anti-inflammatory1",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.2)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","IFNG_sig")]),x="Rx",y="value",xlab="Treatment",ylab="IFNG_sig",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.2)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","IMM_sig")]),x="Rx",y="value",xlab="Treatment",ylab="IMM_sig",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.2)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","TCELL_INF_GEP")]),x="Rx",y="value",xlab="Treatment",ylab="TCELL_INF_GEP",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.2)

################################################
## Compute activation and exhaustion signatures
################################################

## ACTIVATION vs EXHAUSTION
lag3 = t(data.frame("Lag3"=cca_out2@scale.data[grep("Lag3$",rownames(cca_out2@scale.data)),]))
ifng = t(data.frame("Ifng"=cca_out2@scale.data[grep("Ifng$",rownames(cca_out2@scale.data)),]))
temp = cca_out2@scale.data[match(cca_out2@var.genes,rownames(cca_out2@scale.data)),]
exh_out = act_out = matrix(0,nrow(temp),2)
for(i in 1:nrow(exh_out)){
  a_cor = cor.test(temp[i,],lag[1,])
  exh_out[i,1]=a_cor$estimate; exh_out[i,2]=a_cor$p.val
  b_cor = cor.test(temp[i,],ifng[1,])
  act_out[i,1]=b_cor$estimate; act_out[i,2]=b_cor$p.val
}
rownames(exh_out) = rownames(act_out) = rownames(temp)
colnames(exh_out) = colnames(act_out) = c("COR","PVAL")

exh_out = data.frame(exh_out)
exh_out$ADJ_P = p.adjust(exh_out$PVAL)
exh_out = exh_out[order(exh_out$ADJ_P,decreasing=F),]
exh_out = rownames(exh_out[exh_out$COR>0,])[1:50]

act_out = data.frame(act_out)
act_out$ADJ_P = p.adjust(act_out$PVAL)
act_out = act_out[order(act_out$ADJ_P,decreasing=F),]
act_out = rownames(act_out[act_out$COR>0,])[1:50]
ex_ac_gsva = gsva(as.matrix(cca_out2@data),list(exh_out,act_out),method="ssgsea")
ex_ac_gsva = data.frame(t(Matrix::t(scale(Matrix::t(ex_ac_gsva)))),check.names=F)
colnames(ex_ac_gsva) = c("Exhaustion_Sig","Activation_Sig")
cca_out2@meta.data = cbind(cca_out2@meta.data,ex_ac_gsva)
FeaturePlot(cca_out2, features.plot = c("Exhaustion_Sig","Activation_Sig"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Exhaustion_Sig")]),x="Rx",y="value",xlab="Treatment",ylab="Exhaustion_Sig",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.5)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Activation_Sig")]),x="Rx",y="value",xlab="Treatment",ylab="Activation_Sig",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.5)
ggscatter(cca_out2@meta.data,x="Exhaustion_Sig",y="Activation_Sig",size=1,col="light blue",color="Rx",pt.size=2)

### DYSFUNCTION vs CYTOTOXIC
dys_genes = convertHumanGeneList(c("LAG3","HAVCR2","PDCD1","PTMS","FAM3C","IFNG","AKAP5","CD7","PHLDA1","ENTPD1","SNAP47","TNS3","CXCL13",
 "RDH10","DGKH","KIR2DL4","LYST","MIR155HG","RAB27A","CSF1","CTLA4","TNFRSF9","CD27","CCL3","ITGAE","PAG1","TNFRSF1B","GALNT1","GBP2",
 "MYO7A"))
cyt_genes = convertHumanGeneList(c("FGFBP2","CX3CR1","FCGR3A","S1PR5","PLAC8","FGR","C1orf21","SPON2","CD300A","TGFBR3","PLEK","S1PR1","EFHD2","KLRF1","FAM65B",
"C1orf162","STK38","SORL1","FCRL6","TRDC","EMP3","CCND3","KLRG1","BIN2","SELL","KLRB1","SAMD3","ARL4C","IL7R","GNLY"))

dys_cyt_ssgsea = gsva(as.matrix(cca_out2@data),list(dys_genes,cyt_genes),method="ssgsea")
dys_cyt_ssgsea = data.frame(t(Matrix::t(scale(Matrix::t(dys_cyt_ssgsea)))),check.names=F)
colnames(dys_cyt_ssgsea) = c("Dysfunction_Sig","Cytotoxic_Sig")
cca_out2@meta.data = cbind(cca_out2@meta.data,dys_cyt_ssgsea)
FeaturePlot(cca_out2, features.plot = c("Dysfunction_Sig","Cytotoxic_Sig"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Dysfunction_Sig")]),x="Rx",y="value",xlab="Treatment",ylab="Dysfunction_Sig",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.5)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Cytotoxic_Sig")]),x="Rx",y="value",xlab="Treatment",ylab="Cytotoxic_Sig",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.5)

### Rx expansion vs loss
Rx_expanded = markers2[markers2$cluster=="0" & markers2$avg_logFC>0,"gene"][1:50]
Rx_loss = markers2[markers2$cluster=="5" & markers2$avg_logFC>0,"gene"][1:50]
Rx_exp_loss_ssgsea = gsva(as.matrix(cca_out2@data),list(Rx_expanded,Rx_loss),method="ssgsea")
Rx_exp_loss_ssgsea = data.frame(t(Matrix::t(scale(Matrix::t(Rx_exp_loss_ssgsea)))),check.names=F)
colnames(Rx_exp_loss_ssgsea) = c("Rx_Expanded","Rx_Loss")
cca_out2@meta.data = cbind(cca_out2@meta.data,t(Rx_exp_loss_ssgsea))
names(cca_out2@meta.data)[49:50]=c("Rx_Expanded","Rx_Loss")
FeaturePlot(cca_out2, features.plot = c("Rx_Expanded","Rx_Loss"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Rx_Expanded")]),x="Rx",y="value",xlab="Treatment",ylab="Rx_Expansion",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.5)
ggbarplot(melt(cca_out2@meta.data[,c("Rx","Rx_Loss")]),x="Rx",y="value",xlab="Treatment",ylab="Rx_Loss",fill="light blue",add="mean_se") + stat_compare_means(method = "anova",label.y = 0.5)
#######
# ggscatter(cca_out2@meta.data,x="Exhaustion_Sig",y="Activation_Sig",size=1,col="light red",color="IFNG_sig")
# cca_out2@meta.data$HAVCR2 = as.numeric(scale(as.numeric(cca_out2@scale.data[grep("Havcr2$",rownames(cca_out2@scale.data)),])))
# cca_out2@meta.data$IFNG = as.numeric(scale(as.numeric(cca_out2@scale.data[grep("Ifng$",rownames(cca_out2@scale.data)),])))
# ggscatter(cca_out2@meta.data,x="Exhaustion_Sig",y="Activation_Sig",size=1,color="HAVCR2")
# ggscatter(cca_out2@meta.data,x="Exhaustion_Sig",y="Activation_Sig",size=1,color="IFNG")+facet_grid(~Rx)

DimPlot(object = cca_out2, reduction.use = 'umap',do.label=T,label.size=5,pt.size=2)
cca_out2@meta.data$new_groups = ifelse(cca_out2@meta.data$"res.0.6" %in% c(5,7),"Exhausted_Component","Activation_Component")

mod = model.matrix(~Rx+group+new_groups,cca_out2@meta.data)
mod_dat = data.frame(as.matrix(cca_out2@data[match(cca_out2@var.genes,rownames(cca_out2@data)),]),check.names=F)
fit = eBayes(lmFit(mod_dat,mod))
top_gene_tab = topTable(fit,number=nrow(mod_dat),coef="new_groupsExhausted_Component")
top_sig = top_gene_tab[top_gene_tab$adj.P.Val<=0.01,]
top_sig$Genes = rownames(top_sig)
top_sig$DIR = ifelse(top_sig$logFC<0,"Down Reg","Up Reg")
top_sig = top_sig[order(top_sig$logFC,decreasing=T),]
top_sig$Genes = factor(top_sig$Genes,levels=top_sig$Genes)
ggplot(rbind(top_sig[1:25,],top_sig[(nrow(top_sig)-25):nrow(top_sig),])
      ,aes(x=Genes,y=logFC,label=logFC))+geom_bar(stat='identity', aes(fill=DIR),width=.5)+coord_flip() + theme_classic()

#######################################
## Run diffusion component analysis
##################################################

# DE_genes = as.character(top_sig$Genes)
# cca_out3 = cca_out2
# cca_out3@raw.data = cca_out2@raw.data[match(DE_genes,rownames(cca_out2@raw.data)),]
# cca_out3@data = cca_out2@data[match(DE_genes,rownames(cca_out2@data)),]
# cca_out3@scale.data = cca_out2@scale.data[match(DE_genes,rownames(cca_out2@scale.data)),]
# cca_out3 = RunDiffusion(cca_out3,reduction.use = "cca.aligned",dims.use = 1:5)
# DimPlot(object = cca_out3, reduction.use = 'dm')
# DimPlot(object = cca_out3, reduction.use = 'dm',group.by="Rx")
#

# clones$Rx = ifelse(clones$group %in% c("401AB","464A"),"Control",ifelse(clones$group %in% c("403A","403B"),"aPD1","Vaccine_aPD1"))
# clones$Rx = factor(clones$Rx,levels=c("Control","aPD1","Vaccine_aPD1"))
# clones = clones[order(clones$clonotype),]
# table(clones$clonotype,clones$group)
# clones_melt1 = melt(table(clones$clonotype,clones$group))
# colnames(clones_melt1)[1:2]=c("clonotypes","group")
# plot1 = ggplot(clones_melt1, aes(clonotypes, group)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low ="black", high ="magenta", mid ="black",space = "Lab") +
#       theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + coord_flip()
# clones_melt2 = melt(table(clones$clonotype,clones$Rx))
# colnames(clones_melt2)[1:2]=c("clonotypes","Rx")
# plot2 = ggplot(clones_melt2, aes(clonotypes, Rx)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low ="black", high ="magenta", mid ="black",space = "Lab") +
#       theme_minimal() + theme(axis.text.x=element_text(size=8,angle=90)) + coord_flip()
# ggarrange(plot1,plot2,ncol=2,nrow=1)



save.image("datERIKA.rds")
