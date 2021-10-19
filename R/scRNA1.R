#' load counts
#' @description Load counts from datapath and return a seurat object
#' @param datapath character, indicating data path of counts, if it is 10X data,
#' please input the path of folder which containing matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz
#' @param datatype 10X:datatype = 0, counts:datatype = 1, sparse matrix:datatype = 2
#' @importFrom utils read.table
#' @importFrom utils read.delim
#' @return return a seurat object
#' @export
#'
#' @examples \dontrun{load_counts('D:/scRNA/10X',datatype=1)}
load_counts <- function(datapath, datatype = 0) {
  if (datatype[1] == 1) {
    RAW1 <- read.table(datapath, sep = "\t", header = TRUE, row.names = 1)
  }
  if (datatype[1] == 0) {
    RAW1 <- Seurat::Read10X(data.dir = datapath)
  }
  if (datatype[1] == 2) {
    matrix_dir <- datapath
    DirFile1 <- list.files(datapath)
    if ('barcodes.tsv.gz' %in% DirFile1) {
      barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
      features.path <- paste0(matrix_dir, "/features.tsv.gz")
      matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
    }else if('barcodes.tsv' %in% DirFile1){
      barcode.path <- paste0(matrix_dir, "/barcodes.tsv")
      features.path <- paste0(matrix_dir, "/features.tsv")
      matrix.path <- paste0(matrix_dir, "/matrix.mtx")
    }
    RAW1 <- readMM(file = matrix.path)
    feature.names <- read.delim(features.path,
                                header = FALSE,
                                stringsAsFactors = FALSE
    )
    barcode.names <- read.delim(barcode.path,
                                header = FALSE,
                                stringsAsFactors = FALSE
    )
    colnames(RAW1) <- barcode.names$V1
    rownames(RAW1) <- make.unique(feature.names$V2, sep = "_")
    rm(feature.names)
    rm(barcode.names)
  }
  sampproj <- Seurat::CreateSeuratObject(counts = RAW1, min.cells = 5)
  return(sampproj)
}

#' Calculate pseudotime of cells
#' @description Use monocle to calculate the pseudotime and return a monocle object
#' @param seurat_object seurat object
#' @param reverse TRUE or FALSE, whether to reverse the pseudotime, default is TURE
#' @param gene.use vector, indicating the variable genes to calculate pseudotime
#' @import monocle
#' @import pbapply
#' @import ROCR
#' @importFrom monocle newCellDataSet
#' @importFrom BiocGenerics estimateSizeFactors
#' @importFrom BiocGenerics estimateDispersions
#' @importFrom monocle reduceDimension
#' @importFrom monocle orderCells
#' @importFrom VGAM negbinomial.size
#' @importFrom DDRTree DDRTree
#' @return return monocle object which contain pseudotime
#' @export
#'
#' @examples load(system.file("extdata", "test_seurat.rda", package = "IReNA"))
#' get_pseudotime(test_seurat)
get_pseudotime <- function(seurat_object, reverse = FALSE,gene.use = NULL) {
  if (is.null(gene.use)) {
    seurat_object <- Seurat::FindVariableFeatures(seurat_object)
    seurat_object2 <- seurat_object[seurat_object@assays$RNA@var.features,]
  } else{seurat_object2 <- seurat_object[gene.use,]}
  data <- as(as.matrix(seurat_object@assays$RNA@counts), "sparseMatrix")
  pd <- new("AnnotatedDataFrame", data = seurat_object@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)
  monocle_cds <- monocle::newCellDataSet(data,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.5,
                                         expressionFamily = VGAM::negbinomial.size()
  )
  cds <- BiocGenerics::estimateSizeFactors(monocle_cds)
  data2 <- as(as.matrix(seurat_object2@assays$RNA@counts), "sparseMatrix")
  pd2 <- new("AnnotatedDataFrame", data = seurat_object2@meta.data)
  fData2 <- data.frame(gene_short_name = row.names(data2), row.names = row.names(data2))
  fd2 <- new("AnnotatedDataFrame", data = fData2)
  monocle_cds2 <- monocle::newCellDataSet(data2,
                                          phenoData = pd2,
                                          featureData = fd2,
                                          lowerDetectionLimit = 0.5,
                                          expressionFamily = VGAM::negbinomial.size()
  )
  cds2 <- BiocGenerics::estimateSizeFactors(monocle_cds2)
  cds2 <- monocle::reduceDimension(cds2,
                                   max_components = 2,
                                   method = "DDRTree"
  )
  cds2 <- monocle::orderCells(cds2, reverse = reverse)
  cds@dim_reduce_type <- cds2@dim_reduce_type
  cds@cellPairwiseDistances <- cds2@cellPairwiseDistances
  cds@minSpanningTree <- cds2@minSpanningTree
  cds@reducedDimS <- cds2@reducedDimS
  cds@reducedDimW <- cds2@reducedDimW
  cds@reducedDimK <- cds2@reducedDimK
  cds@phenoData@data <- cds2@phenoData@data
  return(cds)
}



#' add pseudotime in monocle object to the metadata of the seurat object
#'
#' @param seurat_object seurat object, which has same cells and genes as monocle_object
#' @param monocle_object monocle object, which include pesudotime and has same
#' cells and genes as seurat_object
#' @return seurat_object with pseudotime
#' @export
#'
#' @examples load(system.file("extdata", "test_seurat.rda", package = "IReNA"))
#' monocle_object = get_pseudotime(test_seurat)
#' add_pseudotime_DEG_filter(seurat_object = test_seurat,monocle_object = monocle_object)
add_pseudotime <- function(seurat_object, monocle_object){
  se <- seurat_object
  mo <- monocle_object
  #### add_pseduotime
  pse <- mo@phenoData@data
  meta1 <- se@meta.data
  filter1 <- rownames(meta1)
  pse <- pse[filter1, ]
  se[["Pseudotime"]] <- pse$Pseudotime
  se[["State"]] <- pse$State
  return(se)
}


diffgenetest_pseudotime <- function(monocle_object){
  mo <- detectGenes(monocle_object, min_expr = 3)
  mo <- estimateDispersions(mo)
  diff1 <- monocle::differentialGeneTest(mo,
                                         fullModelFormulaStr = "~Pseudotime",
                                         relative_expr = TRUE
  )
  ed <- c()
  for (i in diff1$num_cells_expressed) {
    a <- log(i) * 0.95 - log(i) * 0.05
    ed <- c(ed, a)
  }
  diff1$expression_difference <- ed
  return(diff1)
}
