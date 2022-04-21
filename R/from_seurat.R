#' Title
#'
#' @param seuratObj A Seurat object
#' @param data.path path for 10x output directory
#' @param export.path path for output export
#' @param clusterNames names of the clusters
#' @param clusterGenes genes related to the clusters. This should be a list of arrays.
#' @param top_5_genes whether to be shown top fir genes at the spot on mouseover
#' @param multiSamplePattern the patterns to separate the multiple samples
#' @param verbose whether to print the process on colsole
#'
#' @return
#' @export
#' @examples
prepare_from_seurat <- function(seuratObj, data.path, export.path,
                                    clusterNames = NULL,
                                    clusterGenes = NULL,
                                    top_5_genes = TRUE,
                                    multiSamplePattern = "_",
                                    verbose = FALSE){

  stopifnot('seuratObj must be a Seurat object' = "SeuratObject" %in% class(seuratObj))

  seurat.idents <- unique(seuratObj$orig.ident)
  sel.ident <- NA
  if(length(seurat.idents) > 1){
    cat("Multiple indents found in the Seurat object. \n")
    cat(paste(1:length(seurat.idents), seurat.idents), sep = "\n")
    sel.ident.id <- readline(prompt="Select the matching indent: ")
    sel.ident.id <- as.integer(sel.ident.id)
    if(sel.ident.id > 0 & sel.ident.id <= length(seurat.idents)){
      sel.ident <- seurat.idents[sel.ident.id]
    }else{
      stop("Improper indent selection...")
    }
  }else{
    sel.ident = seurat.idents[0]
  }

  seuratObj.temp <- subset(seuratObj,subset = orig.ident == sel.ident)

  #step 0
  #create the export directory
  data_id = basename(data.path)
  export.path = file.path(export.path, data_id)
  dir.create(export.path, showWarnings = FALSE)
  #Step 1
  #checking scalefactors_json.json, this file is a must have one else error
  if(! file.exists(file.path(data.path,"spatial","scalefactors_json.json"))){
    if(file.exists(file.path(data.path,"spatial","scalefactors_json.json.gz"))){
      require(R.utils,quietly = TRUE)
      gunzip(file.path(data.path,"spatial","scalefactors_json.json.gz"), remove=FALSE)
    }
    else{
      stop(paste0("scalefactors_json.json file not found at ",
                  file.path(data.path,"spatial","scalefactors_json.json")))
    }}
  file.copy(from=file.path(data.path,"spatial","scalefactors_json.json"),
            to=export.path,
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)

  #checking tissue_positions_list.sv, this file is a must have one else error
  if(! file.exists(file.path(data.path,"spatial","tissue_positions_list.csv"))){
    if(file.exists(file.path(data.path,"spatial","tissue_positions_list.csv.gz"))){
      require(R.utils,quietly = TRUE)
      gunzip(file.path(data.path,"spatial","tissue_positions_list.csv.gz"), remove=FALSE)
    }
    else{
      stop(paste0("tissue_positions_list.csv file not found at ",
                  file.path(data.path,"spatial","scalefactors_json.json")))
    }}
  file.copy(from=file.path(data.path,"spatial","tissue_positions_list.csv"),
            to=data.path,
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)
  #
  # #Step 2
  # #Checking the presence of filtered_feature_bc_matrix.h5, if this file is not present then...can we create one
  # if(! file.exists(file.path(data.path,"filtered_feature_bc_matrix.h5"))){
  #   #filtered_feature_bc_matrix.h5 doesn't exist
  #   if(verbose)print("filtered_feature_bc_matrix.h5 doesn't exist. Trying to create it from the barcode matrix...")
  #
  #   #Creating the filtered_feature_bc_matrix.h5
  #   if((file.exists(file.path(data.path,"filtered_feature_bc_matrix","barcodes.tsv.gz")) |
  #       file.exists(file.path(data.path,"filtered_feature_bc_matrix","barcodes.tsv"))) &
  #      (file.exists(file.path(data.path,"filtered_feature_bc_matrix","features.tsv")) |
  #       file.exists(file.path(data.path,"filtered_feature_bc_matrix","features.tsv.gz"))) &
  #      (file.exists(file.path(data.path,"filtered_feature_bc_matrix","matrix.mtx"))) |
  #      (file.exists(file.path(data.path,"filtered_feature_bc_matrix","matrix.mtx.gz")))){
  #
  #     write10X_to_h5(data.path)
  #   }else{
  #     stop(paste0("Either filtered_feature_bc_matrix.h5 file exptected at ",data_id," folder or barcodes.tsv.gz,
  #                 barcodes.tsv.gz, matrix.mtx.gz files are expected at ",file.path(data.path,"filtered_feature_bc_matrix")))
  #
  #   }
  # }else{
  #   if(verbose)print("filtered_feature_bc_matrix.h5 file found")
  # }
  #
  #Step 3
  #Check either tissue_lowres_image.png or tissue_hires_image.png image present
  if(!(file.exists(file.path(data.path,"spatial","tissue_lowres_image.png")) |
       file.exists(file.path(data.path,"spatial","tissue_hires_image.png")))){
    stop(paste0("Expected either " ,file.path(data.path,"spatial","tissue_lowres_image.png")," or ",
                file.path(data.path,"spatial","tissue_hires_image.png"),
                ", but not found"))
  }

  #Step 3.a
  #Checking if tissue_lowres_image.png present or not, if not then trying to use tissue_hires_image.png
  if((!file.exists(file.path(data.path,"spatial","tissue_lowres_image.png"))) &
     file.exists(file.path(data.path,"spatial","tissue_hires_image.png"))){

    #copy tissue_hires_image.png as tissue_lowres_image.png
    file.copy(file.path(data.path,"spatial","tissue_hires_image.png"),
              file.path(data.path,"spatial","tissue_lowres_image.png"), overwrite = FALSE )


    #updating the scalefactors_json.json file
    require(rjson)
    scalefactors <- fromJSON(file = file.path(data.path,"spatial","scalefactors_json.json"))

    scalefactors$tissue_lowres_scalef = scalefactors$tissue_hires_scalef
    write(toJSON(scalefactors), file.path(export.path,"spatial","scalefactors_json.json"))
  }

  #Step 3.b
  #If tissue_hires_image.png is not present, then trying to use tissue_lowres_image.png
  if((!file.exists(file.path(data.path,"spatial","tissue_hires_image.png"))) &
     file.exists(file.path(data.path,"spatial","tissue_lowres_image.png"))){

    #copy tissue_lowres_image.png as tissue_hires_image.png
    file.copy(file.path(data.path,"spatial","tissue_lowres_image.png"),
              file.path(data.path,"spatial","tissue_hires_image.png"), overwrite = FALSE )


    #updating the scalefactors_json.json file
    scalefactors$tissue_hires_scalef = scalefactors$tissue_lowres_scalef
    write(toJSON(scalefactors), file.path(export.path,"spatial","scalefactors_json.json"))
  }

  file.copy(from=file.path(data.path,"spatial","tissue_hires_image.png"),
            to=file.path(export.path,"tissue_hires_image.png"),
            overwrite = TRUE,
            copy.mode = TRUE)

  #creating tissue info
  tissue_info <- data.frame(sample = c(NA),age = c(NA), gender = c(NA), age_cat = c(NA))
  write.csv(tissue_info, file.path(export.path,"sample_info.csv"), row.names = FALSE,
            quote = TRUE)


  #creating metadata
  metadata <- seuratObj.temp@meta.data
  if("seurat_clusters" %in% names(metadata)){
    metadata$cluster <- metadata$seurat_clusters
  }else{
    metadata$cluster <- 1
  }

  metadata$barcode <-  str_split(rownames(metadata),pattern =  multiSamplePattern, simplify = TRUE)[,1]

  spot_info <- read.csv(file.path(data.path,"/spatial/tissue_positions_list.csv"),header = FALSE)
  colnames(spot_info) <- c("barcode","in_tissue","array_row",
                           "array_col","pxl_row_in_fullres","pxl_col_in_fullres")
  metadata <- merge(x = metadata, y = spot_info, by = "barcode")

  write.csv(metadata,
            file.path(export.path,"metadata.csv"),
            quote = TRUE, row.names = FALSE)


  #creating top-5 genes
  normalized_counts <- GetAssayData(seuratObj.temp, slot = "data")
  colnames(normalized_counts) <-  str_split(colnames(normalized_counts),
                                            pattern = multiSamplePattern,
                                            simplify = TRUE)[,1]
  top_5_genes <- get_top_n_genes(normalized_counts, 5)

  write.csv(top_5_genes,
            file.path(export.path,"top_5_genes.csv"),
            quote = TRUE, row.names = FALSE)


  #creating normalized_counts data
  normalized_counts <- as.data.frame(t(as.matrix(normalized_counts)))
  normalized_counts <- round(normalized_counts, 3)

  final_de_genes <- unique(unlist(clusterGenes))
  normalized_counts.store <- matrix(data=rep(0,nrow(normalized_counts) * length(final_de_genes)),
                                    nrow = nrow(normalized_counts))
  normalized_counts.store <- as.data.frame(normalized_counts.store)
  rownames(normalized_counts.store) <- rownames(normalized_counts)
  colnames(normalized_counts.store) <- final_de_genes

  final_de_genes <- final_de_genes[final_de_genes %in% colnames(normalized_counts)]
  normalized_counts.store[,final_de_genes] <- normalized_counts[,final_de_genes]
  normalized_counts.store$barcode = rownames(normalized_counts.store)
  write.csv(normalized_counts.store,
            file.path(export.path,"normalized_marker_expr.csv"),
            row.names = FALSE, quote = TRUE)

  #creating cluster  info
  unique_clusters <- sort(unique(metadata$cluster))
  require(ggplot2)
  a <- scale_color_hue(h.start=180)
  clust_cols <- a$palette(length(unique_clusters))#number of colors

  if(is.null(clusterNames)) clusterNames <- rep("undefined",length(unique_clusters))
  if(is.null(clusterGenes)) clusterGenes <- rep("undefined",length(unique_clusters))

  clusterInfo.df <- data.frame(cluster = unique_clusters,
                               color = clust_cols,
                               name = clusterNames,
                               genes = unlist(lapply(clusterGenes, toString)))
  write.csv(clusterInfo.df,
            file.path(export.path,"cluster_info.csv"),
            row.names = FALSE, quote = TRUE)

  expr_grid_orientation <- prepare_expr_grid(clusterGenes, n.col = 15)
  write.csv(expr_grid_orientation,
            file.path(export.path,"marker_grid_position.csv"),
            row.names = FALSE, quote = TRUE)

  remove(seuratObj.temp)
  if(verbose)print("All required files are found/prepared.")
}
