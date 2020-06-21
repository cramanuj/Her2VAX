################
## Read HD5 files
##################

read_h5 = function(myh5_path){
    mydat = h5read(myh5_path,name="/matrix")
    raw_data = mydat$data
    indices = mydat$indices
    indptr = mydat$indptr
    shp = mydat$shape
    barcodes = mydat$barcodes
    gene_names = do.call("cbind",mydat$features)
    gene_names = data.frame(gene_names[,c(4:5)])
    sparse_mat = sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(raw_data[]), dims = shp[], giveCsparse = FALSE)
    rownames(sparse_mat) <- gene_names$id
    colnames(sparse_mat) <- barcodes
    mat = as.matrix(sparse_mat)
    dim(mat)
    out = list(mat=mat,feat=gene_names)
    return(out)
    h5closeAll()
}

##########################
## Read gene-set files
##########################

read.gmt.file = function(filename){
	a = scan(filename, what = list("", ""), sep = "\t", quiet=T,quote = NULL, fill = T, flush = T, multi.line = F)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(filename, what = "", sep = "\t", quote = NULL,quiet=T)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
      while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] !=
          geneset.descriptions[i])) {
          ii = ii + 1
    	}
      ox[i] = ii
      ii = ii + 1
  }
  genesets = vector("list", nn)
  for (i in 1:(nn - 1)) {
      i1 = ox[i] + 2
      i2 = ox[i + 1] - 1
      geneset.descriptions[i] = dd[ox[i] + 1]
      genesets[[i]] = dd[i1:i2]
  }
  geneset.descriptions[nn] = dd[ox[nn] + 1]
  genesets[[nn]] = dd[(ox[nn] + 2):n]
  out = list(genesets = genesets, geneset.names = geneset.names, geneset.descriptions = geneset.descriptions)
	names(out$genesets) = out$geneset.names
	for(i in 1:length(out$genesets)){
		aa = out$genesets[[i]]
		aa = aa[aa!=""]
		out$genesets[[i]] = aa
	}
	cat("Number of genesets read: ",length(out$genesets),"\n")
	return(out)
}
