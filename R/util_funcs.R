get_top_n_genes <- function(exprData, n= 3){
  top_ns_exprs <- apply(exprData, MARGIN = 2,function(x){
    x <- sort(x, decreasing = TRUE)
    n_genes <- names(x[1:n])
    names(x) <- NULL
    row.data <- c(n_genes, round(x[1:n],3))
    row.data <- data.frame(c=row.data)
    rownames(row.data) <- c(paste0(rep('gene',n), 1:n), paste0(rep('expr',n), 1:n))
    row.data
  })

  top_ns_exprs <- bind_cols(top_ns_exprs)
  top_ns_exprs <- t(top_ns_exprs)

  rownames(top_ns_exprs) <- colnames(exprData)
  top_ns_exprs <- as.data.frame(top_ns_exprs)
  top_ns_exprs$barcode <- colnames(exprData)
  top_ns_exprs
}

prepare_expr_grid <- function(genes_list, n.col = 15){
  genes_arr <- c()
  cluster_arr <- c()
  row_arr <- c()
  col_arr <- c()

  current_row = 0
  for(i in 1:length(genes_list)){
    clust.i_genes <- genes_list[[i]]
    clust_genes_count <- length(clust.i_genes)

    clust.i <- rep(i, clust_genes_count)

    clust.i_row <- current_row + seq(1,ceiling(clust_genes_count/n.col))
    clust.i_row <- sort(rep(clust.i_row, n.col))[1:clust_genes_count]
    current_row <- clust.i_row[length(clust.i_row)]

    clust.i_col <- rep(1:n.col, ceiling(clust_genes_count/n.col))[1:clust_genes_count]

    genes_arr <- c(genes_arr, clust.i_genes)
    cluster_arr <- c(cluster_arr, clust.i)
    row_arr <- c(row_arr, clust.i_row)
    col_arr <- c(col_arr, clust.i_col)
  }

  return(data.frame(cluster = cluster_arr, gene = genes_arr,
                    row = row_arr, col = col_arr,
                    expr = 0))

}
