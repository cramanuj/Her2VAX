read_h5 = function(myh5_path){
    require(rhdf5); require(Matrix)
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
