###########################################################################################
### R code to analyze single cell RNAseq data
### Authored by Chaitanya Acharya
### Updated on: June 16, 2020
###########################################################################################
lib.list = c( "Seurat","data.table","plyr","dplyr","ggplot2","GSVA","pheatmap","ggthemes","ggpubr","clustree",
              "psych","RColorBrewer","limma","biomaRt","tidyr","gtools","gridExtra","Matrix","rhdf5","car")
for(i in 1:length(lib.list)){
  if(any(installed.packages()[,1]==lib.list[i])){
    require(lib.list[i],character.only=T,quietly = T,warn.conflicts=F)
  }
  else{
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(lib.list[i]);
    require(lib.list[i],character.only=T,quietly = T,warn.conflicts=F)
  }
}

## Load utility functions
source("util_funcs.R")

# source("~/Research/scripts/r_scripts/useful_functions.R")
# source("~/Research/scripts/r_scripts/plotfns.R")
# source("~/Research/Zack/read_hi5file.R")

## Define color palette
my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
col_breaks = c( seq(-1,0,length=50),
                seq(0.01,0.8,length=50),
                seq(0.81,1,length=50))

### Path to data input/results
input_path = "data/"
output_path = "results/"

#########################################################################################################
### T Cell Receptor processing
#########################################################################################################
TCRfiles = list.files(pattern="filtered_contig_annotations.csv",recursive=T,include.dirs=T,full.names=T)
tcr_dat=list()
for(i in 1:length(TCRfiles)){
  print(i)
  a = read.csv(TCRfiles[i],header=T)
  if(length(grep("None",a$raw_clonotype_id))>0){
    a = a[-grep("None",a$raw_clonotype_id),]
  }
  a$group = gsub("_result","",strsplit(TCRfiles[i],"/")[[1]][4])
  tcr_dat[[i]] = a
}
tcr_dat = do.call("rbind",tcr_dat)
tcr_dat$new_clonotype = paste(tcr_dat$raw_clonotype_id,tcr_dat$group,sep="_")

tcr_dat_new = tcr_dat[!duplicated(tcr_dat$barcode),]
tcr_dat_new$raw_clonotype_id = factor(tcr_dat_new$raw_clonotype_id,levels= mixedsort(names(table(tcr_dat_new$raw_clonotype_id))))
clone_table  = table(tcr_dat_new$raw_clonotype_id,tcr_dat_new$group)
write.table(clone_table,paste0(output_path,"clonotype_count_table.txt"),sep="\t",col.names=NA,quote=F)
write.table(round(prop.table(clone_table,2)*100,2),paste0(output_path,"clonotype_prop_table.txt"),sep="\t",col.names=NA,quote=F)

######################################################################################################
### scRNAseq processing
###
######################################################################################################

file1 = read_h5(paste0(input_path,"GEX/GEX_401AB/filtered_feature_bc_matrix.h5"))
file2 = read_h5(paste0(input_path,"GEX/GEX_403A/filtered_feature_bc_matrix.h5"))
file3 = read_h5(paste0(input_path,"GEX/GEX_403B/filtered_feature_bc_matrix.h5"))
file4 = read_h5(paste0(input_path,"GEX/GEX_464/filtered_feature_bc_matrix.h5"))
file5 = read_h5(paste0(input_path,"GEX/GEX_471/filtered_feature_bc_matrix.h5"))
file6 = read_h5(paste0(input_path,"GEX/GEX_663/filtered_feature_bc_matrix.h5"))
file7 = read_h5(paste0(input_path,"GEX/GEX_668/filtered_feature_bc_matrix.h5"))
file8 = read_h5(paste0(input_path,"GEX/GEX_669/filtered_feature_bc_matrix.h5"))
file9 = read_h5(paste0(input_path,"GEX/GEX_672/filtered_feature_bc_matrix.h5"))
data_list = list(file1$mat,cbind(file2$mat, file3$mat), file4$mat,file5$mat,file6$mat,file7$mat,file8$mat,file9$mat)
files = c("GEX_401AB","GEX_403AB","GEX_464","GEX_471","GEX_663","GEX_668","GEX_669","GEX_672")
g2m.features = scan(paste0(input_path,"g2m.features.txt"),what="",quiet=T)
s.features = scan(paste0(input_path,"s.features.txt"),what="",quiet=T)

### Create a SEURAT object from each gene expression matrix within the list
### Each SEURAT object is preprocessed separately and highly variable genes are computed
seu_list = vector(mode="list",length=length(data_list))
for(i in 1:length(seu_list)){
  print(i)
  tmp = data.frame(data_list[[i]],check.names=F)
  tmp$Genes = as.character(file1$feat$name)
  tmp = setDT(tmp)[, lapply(.SD, sum), by = Genes] %>%data.frame(check.names=F)
  rownames(tmp) = tmp[,1]
  tmp = tmp[,-1]
  seu = CreateSeuratObject(counts = tmp, min.cells = 3, min.features = ncol(tmp)*0.02, project = files[i])
  seu[["percent.mt"]]=PercentageFeatureSet(seu,pattern="^mt")
  seu = subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < quantile(seu@meta.data$nFeature_RNA,0.99) & nCount_RNA < quantile(seu@meta.data$nCount_RNA,0.95) & percent.mt < 10)
  seu = NormalizeData(seu,verbose=F)
  seu = FindVariableFeatures(seu,selection.method="vst",nfeatures=5000,verbose=F)
  seu = CellCycleScoring(seu,g2m.features=g2m.features,s.features=s.features)
  seu = ScaleData(seu,vars.to.regress=c("nCount_RNA","S.Score","G2M.Score","percent.mt"),verbose=F)
  seu_list[[i]] = seu
}
cell_names = llply(seu_list,colnames)
rm(file1,file2,file3,file4,file5,file6,file7,file8,file9)
gc()

sum(laply(1:length(seu_list),function(i) ncol(seu_list[[i]]@assays$RNA@data)))
pheno = c("Isotype","aPD1","Isotype","Her2VAC_aPD1","aPD1","Her2VAC","Her2VAC_aPD1","Isotype")

#######################################################################################
### 1) Select data integration featues
### 2) Find intergration anchors
### 3) Integrate data using reference datasets
### 4) Scale intergrated data
#######################################################################################

# pdf("CCA_plots.pdf",width=10,height=7)
seu_features = SelectIntegrationFeatures(object.list = seu_list, nfeatures = 5000)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset = which(names(seu_list) %in% c("GEX_401AB","GEX_464","GEX_672"))
anchors = FindIntegrationAnchors(object.list = seu_list, normalization.method = "LogNormalize", anchor.features = seu_features, reference = reference_dataset,dims = 1:30)
dat_int = IntegrateData(anchorset = anchors, dims = 1:30)
dat_int = ScaleData(dat_int, verbose = FALSE)
dat_int@meta.data$Rx = factor(rep(pheno,as.numeric(table(dat_int@meta.data$orig.ident))))
dat_int@meta.data$Rx = factor(dat_int@meta.data$Rx,levels=c("Isotype","aPD1","Her2VAC","Her2VAC_aPD1"))

### We remove the following genes in order to eliminate any biased clustering
REMOVEgenes = c(VariableFeatures(dat_int)[grep("Tra",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Trb",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Trg",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Trd",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Tcrg",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Igh",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Igk",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Igl",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Rps",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Rpl",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("^mt",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Mrps",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Mrpl",VariableFeatures(dat_int))],
                VariableFeatures(dat_int)[grep("Lactb",VariableFeatures(dat_int))])
VariableFeatures(dat_int) = VariableFeatures(dat_int)[!VariableFeatures(dat_int) %in% REMOVEgenes]
dat_int@assays$integrated@data = dat_int@assays$integrated@data[match(VariableFeatures(dat_int),rownames(dat_int@assays$integrated@data)),]
tmp_barcode = rownames(dat_int@meta.data)
tmp_barcode = laply(1:length(tmp_barcode),function(i) paste(strsplit(tmp_barcode[i],split="_")[[1]][1],collapse="-"))
dat_int@meta.data$new_barcode = tmp_barcode
dat_int = subset(dat_int,cells=rownames(dat_int@meta.data)[which(!duplicated(tmp_barcode))])

dat_int@meta.data$CLONES = (dat_int@meta.data$new_barcode %in% tcr_dat$barcode) + 0
dat_int@meta.data$CLONES = ifelse(dat_int@meta.data$CLONES>0,"TCR","No TCR")

##########################################################
### 1) Dimensionality reduction usinng PCA+UMAP, and
### 2) Graoh-based clustering
##########################################################
pdf(paste0(output_path,"original_clusters.pdf"),width=10,height=7)
dat_int = RunPCA(object = dat_int, verbose = FALSE)
dat_int = RunUMAP(object = dat_int, dims = 1:20, umap.method="umap-learn",n.neighbors = 50L,min.dist=1,metric="cosine")
dat_int = FindNeighbors(dat_int, reduction = "pca", dims = 1:10,features = VariableFeatures(object = dat_int))
dat_int = FindClusters(dat_int,resolution=c(0.4,0.6),random.seed=123)
DimPlot(dat_int,combine = FALSE,label=T)
DimPlot(dat_int,combine = FALSE,label=T,group.by="Rx")
DimPlot(dat_int,combine = FALSE,label=T,group.by="CLONES")

########################
### Macrophages
### Macrophages (unexpectedly) dominated the cell types.
### Our objective here is to find clusters dominated by macrophages.
#########################
macrophages = c("Apoe","H2-Aa","C1qa","C1qb","C1qc")
macro_zscore = gsva(as.matrix(dat_int@assays$integrated@data),list(macrophages),method="zscore")

macro_out = data.frame("scores"=as.numeric(scale(as.numeric(macro_zscore))),"clusters"=Idents(dat_int))
macro_ave = aggregate(scores ~ clusters,macro_out,"mean")
ggbarplot(macro_ave,x="clusters",y="scores",xlab="",ylab="Macrophage Z-score",fill="black")+geom_hline(yintercept=quantile(macro_ave$scores,0.5),col="red",linetype="dashed")+geom_hline(yintercept=0,col="black")
keep_out = as.numeric(as.character(macro_ave[macro_ave$scores>quantile(macro_ave$scores,0.5),"clusters"]))

###########################################
## Feature plots to label cell types
###########################################
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

########################################
### Clonal expansion
########################################
tcr_merged = base::merge(dat_int@meta.data,tcr_dat_new,by.x=names(dat_int@meta.data)[9],by.y=names(tcr_dat_new)[1],all=T)
dim(tcr_merged)
tcr_merged_new = tcr_merged[match(dat_int@meta.data$new_barcode,tcr_merged$new_barcode),]
rownames(tcr_merged_new)=rownames(dat_int@meta.data)
dat_int@meta.data = tcr_merged_new
clone_table = table(dat_int@meta.data$new_clonotype,dat_int@meta.data$group)
clone_melt = reshape2::melt(clone_table)
clone_melt = clone_melt[clone_melt$value>=10,]
dat_int@meta.data$clonal_expansion = (dat_int@meta.data$new_clonotype %in% clone_melt$Var1) + 0
dat_int@meta.data$clonal_expansion = ifelse(dat_int@meta.data$clonal_expansion>0,">10","<10")
pdf(paste0(output_path,"clonal_expansion.pdf"),width=10,height=7)
DimPlot(dat_int,combine = FALSE,label=T)
DimPlot(dat_int,combine = FALSE,label=F,group.by="CLONES",cols=c("lightgrey","steelblue"))
DimPlot(dat_int,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","steelblue"))
DimPlot(dat_int,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","magenta"),split.by="Rx")


## Cluster relabeling
new_clusters = Idents(dat_int)
new_clusters = recode(new_clusters,"c('0','2','13')='CD4 T cells';c('3','1','8')='CD8 T cells';c('5','7','10','11')='Macrophages';'16'='Tumor cells'; '4'='Tregs';'9'='Dendritic cells'; '6'='Activated CD8 T cells';c('12','14','15')=' '")
dat_int@meta.data$new_clusters = new_clusters
cluster_cols=c("grey","red","blue3","chartreuse4","darkorchid3","burlywood","darkgoldenrod1","chocolate4")
DimPlot(dat_int,combine = FALSE,label=F,group.by="new_clusters",cols=cluster_cols)
dev.off()

## Proportion test
prop.test(table(dat_int@meta.data$new_clusters,dat_int@meta.data$Rx)[2,],colSums(table(dat_int@meta.data$new_clusters,dat_int@meta.data$Rx)))
cluster_tab = table(dat_int@meta.data$new_clusters,dat_int@meta.data$Rx)
cd8_cluster_tab = cluster_tab[c(2,4),]
prop.test(x=cd8_cluster_tab[1,c(1,4)],n=colSums(cd8_cluster_tab)[c(1,4)])

########################################
### Filtering data without macrophages
### PCA + UMAP + Clustering
########################################
pdf(paste0(output_path,"new_groupings.pdf"),width=10,height=7)

dat_int1 = subset(dat_int,idents=setdiff(dat_int@meta.data$seurat_clusters,keep_out),slot='count')

dat_int1 = subset(dat_int1,subset=CLONES%in%c("TCR"))
dat_int1@meta.data = dat_int1@meta.data[,-c(11:13)]
dat_int1@meta.data$CLONES = as.factor(dat_int1@meta.data$CLONES)
dat_int1 = RunPCA(object = dat_int1, verbose = FALSE)
dat_int1 = RunUMAP(object = dat_int1, dims = 1:20, umap.method="umap-learn",n.neighbors = 50L,min.dist=1,metric="cosine")
dat_int1 = FindNeighbors(dat_int1, reduction = "pca", dims = 1:10,features = VariableFeatures(object = dat_int1))
dat_int1 = FindClusters(dat_int1,resolution=c(0.4,0.6,0.8),random.seed=123)
plot.new()
grid.table(table(dat_int1@meta.data$orig.ident,dat_int1@meta.data$Rx))
clustree(dat_int1,prefix="integrated_snn_res.",layout="sugiyama")
clustree(dat_int1, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
Idents(object=dat_int1)<-"seurat_clusters"
dat_int1 = BuildClusterTree(dat_int1)
PlotClusterTree(dat_int1)

DimPlot(dat_int1,combine = FALSE,label=T,pt.size=2)
DimPlot(dat_int1,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","steelblue"),pt.size=1)
DimPlot(dat_int1,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","magenta"),split.by="Rx",pt.size=1)

##########################################################
## Feature plots
##########################################################

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

########################################################################
## Compute single sample enrichment scores for immune-related genesets
########################################################################

azizi_gs = read.gmt.file(paste0(input_path,"azizi_genesets_mouse.gmt"))$genesets[c(1:5,14:16,18)]
acharya_gs = read.gmt.file(paste0(input_path,"acharya_genesets_mouse.gmt"))$genesets[c(1:12,15:16,37:38)]
IFNG_sig = c("Ido1","Cxcl10","Cxcl9","H2-Ea-ps","Stat1","Ifng")
genesets = c(azizi_gs,acharya_gs,"IFNG_sig"=list(IFNG_sig))
names(genesets)=gsub("-","_",names(genesets));
names(genesets)=gsub(" ","_",names(genesets));
genesets_zscore = gsva(as.matrix(dat_int1@assays$integrated@data),genesets,method="zscore")
genesets_zscore =  data.frame(genesets_zscore,check.names=F)
zscores = data.frame(t(genesets_zscore))
zscores$Rx = dat_int1@meta.data$Rx
zscores$cluster = dat_int1@meta.data$seurat_clusters
pdf(paste0(output_path,"genesets_barplots.pdf"),width=10,height=7)
for(i in 1:nrow(genesets_zscore)){
  print(i)
  print(ggbarplot(zscores,x="Rx",y=names(zscores)[i],ylab=names(zscores)[i],add="mean_se",fill="black")+stat_compare_means(method = "anova",label.y = 0.05)+geom_hline(yintercept=0,col="black"))
  print(ggbarplot(zscores,x="cluster",y=names(zscores)[i],add="mean_se",fill="Rx",ylab=names(zscores)[i],position=position_dodge())+geom_hline(yintercept=0,col="black")+geom_vline(xintercept=seq(1.5,14.5,1)))
}
dev.off()

############################################################
## ACTIVATION vs EXHAUSTION in all cells except for Tregs
############################################################
dat_int2 = subset(dat_int1,idents=setdiff(dat_int1@meta.data$integrated_snn_res.0.8,c(3)))
dat_int2@meta.data = dat_int2@meta.data[,c(1:30)]
dat_int2 = RunPCA(object = dat_int2, verbose = FALSE)
dat_int2 = RunUMAP(object = dat_int2, dims = 1:20, umap.method="umap-learn",n.neighbors = 50L,min.dist=1,metric="cosine")
dat_int2 = FindNeighbors(dat_int2, reduction = "pca", dims = 1:10,features = VariableFeatures(object = dat_int2))
dat_int2 = FindClusters(dat_int2,resolution=c(0.4,0.6),random.seed=123)
DimPlot(dat_int2,combine = FALSE,label=T,pt.size=2)
DimPlot(dat_int2,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","steelblue"),pt.size=1)
DimPlot(dat_int2,combine = FALSE,label=F,group.by="clonal_expansion",cols=c("lightgrey","magenta"),split.by="Rx",pt.size=1)

tim3 = t(data.frame("Tim3"=dat_int2@assays$integrated@data[grep("Havcr2$",rownames(dat_int2@assays$integrated@data)),]))
Tnf = t(data.frame("Tnf"=dat_int2@assays$integrated@data[grep("Tnf$",rownames(dat_int2@assays$integrated@data)),]))
temp = as.matrix(dat_int2@assays$integrated@data)

exh_out = act_out = matrix(0,nrow(temp),2)
for(i in 1:nrow(exh_out)){
  a_cor = cor.test(as.numeric(temp[i,]),as.numeric(tim3[1,]))
  exh_out[i,1]=a_cor$estimate; exh_out[i,2]=a_cor$p.val
  b_cor = cor.test(as.numeric(temp[i,]),as.numeric(Tnf[1,]))
  act_out[i,1]=b_cor$estimate; act_out[i,2]=b_cor$p.val
}
rownames(exh_out)  = rownames(act_out) = rownames(temp)
colnames(exh_out)  = colnames(act_out) = c("COR","PVAL")

exh_out = data.frame(exh_out)
exh_out$ADJ_P = p.adjust(exh_out$PVAL)
exh_out = exh_out[order(exh_out$ADJ_P,decreasing=F),]
exh_out = rownames(exh_out[exh_out$COR>0,])[1:50]

act_out = data.frame(act_out)
act_out$ADJ_P = p.adjust(act_out$PVAL)
act_out = act_out[order(act_out$ADJ_P,decreasing=F),]
act_out = rownames(act_out[act_out$COR>0,])[1:50]
dat_int2@meta.data=dat_int2@meta.data[,-c(34:35)]
dat_int2 = AddModuleScore(dat_int2,features=list(exh_out,act_out),ctrl=10,name=c("Exhaustion_score","Activation_score"))
colnames(dat_int2@meta.data)[34:35]=c("Exhaustion_Sig","Activation_Sig")

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


## session info
sessionInfo()

save.image("paste0(datCRz.rds")
