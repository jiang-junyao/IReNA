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
get_pseudotime <- function(seurat_object, reverse = FALSE) {
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
  sc_cds <- BiocGenerics::estimateSizeFactors(monocle_cds)
  cds <- monocle::reduceDimension(sc_cds,
                                  max_components = 2,
                                  method = "DDRTree"
  )
  cds <- monocle::orderCells(cds, reverse = reverse)
  return(cds)
}

#' add pseudtoime and identify DEGs
#' @description Add pseudotime to seurat object, and filter differential
#' expressed genes according to pseudotime
#' @param seurat_object seurat object, which has same cells and genes as monocle_object
#' @param monocle_object monocle object, which include pesudotime and has same
#' cells and genes as seurat_object
#' @param DEG logic, indicating whether filter differential expressed genes for
#' seurat object
#' @param qvalue numeric, indicating q-value to indentify differentially
#' expressed genes
#' @param nce numeric, indicating num cells expressed to indentify differentially
#' expressed genes
#' @param ed numeric, indicating expression difference to indentify differentially
#' expressed genes
#' @param normlize1 logic, indicating whether normalize the data in seurat object
#'
#' @return return a seurat object
#' @export
#'
#' @examples load(system.file("extdata", "test_seurat.rda", package = "IReNA"))
#' monocle_object = get_pseudotime(test_seurat)
#' add_pseudotime_DEG_filter(seurat_object = test_seurat,monocle_object = monocle_object, DEG = FALSE, normlize1 = FALSE)
#' #add_pseudotime_DEG_filter(seurat_object = test_seurat,monocle_object = monocle_object, DEG = TRUE, qvalue = 0.001, nce = 0.1, ed = 0.1)
add_pseudotime_DEG_filter <- function(seurat_object, monocle_object, DEG = TRUE,
                                      qvalue = 0.05, nce = 0.1, ed = 0.1,
                                      normlize1 = TRUE) {
  se <- seurat_object
  mo <- monocle_object
  #### add_pseduotime
  pse <- mo@phenoData@data
  meta1 <- se@meta.data
  filter1 <- rownames(meta1)
  pse <- pse[filter1, ]
  se[["Pseudotime"]] <- pse$Pseudotime
  if (normlize1 == TRUE) {
    se <- Seurat::NormalizeData(object = se)
  }
  if (DEG == TRUE) {
    print("return Seurat object has Pseudotime DEGs")
    mo <- estimateDispersions(mo)
    mo <- detectGenes(mo, min_expr = 3)
    diff1 <- monocle::differentialGeneTest(mo,
                                           fullModelFormulaStr = "~Pseudotime",
                                           relative_expr = TRUE
    )
    sig_genes <- subset(diff1, qval < qvalue)
    sig_genes <- subset(sig_genes, num_cells_expressed > nce)
    b <- c()
    for (i in sig_genes$num_cells_expressed) {
      a <- log(i) * 0.95 - log(i) * 0.05
      b <- c(b, a)
    }
    sig_genes$expression_difference <- b
    sig_genes <- subset(sig_genes, expression_difference > ed)
    pbmc <- subset(se, features = rownames(sig_genes))
    return(pbmc)
  } else {
    return(se)
  }
}

