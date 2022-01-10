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
#' add_pseudotime(seurat_object = test_seurat,monocle_object = monocle_object)
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

#' extract expressed transcription factors
#' @description extract expressed transcription factors according to the proportion of cells that
#' transcription factors express
#' @param seurat_object seurat object
#' @param TFs character, indicating transcription factors names
#' @param cells_quantile numeric, indicating proportion of cells that
#' transcription factors express (threshold to get expressed transcription factors)
#'
#' @return
#' @export
#'
#' @examples load(system.file("extdata", "test_seurat.rda", package = "IReNA"))
#' TFs <- c('CLK1','TCEB3','SOX4')
#' extract_expressed_TFs(test_seurat,TFs)
extract_expressed_TFs <- function(seurat_object,TFs,cells_quantile = 0.05){
  TFs <- TFs[TFs%in%rownames(seurat_object)]
  matrix_tf <- seurat_object@assays$RNA@counts[TFs,]
  if (cells_quantile==0) {
    TfExp <- matrix_tf[rowSums(as.matrix(matrix_tf))>0,]
  }else{
    quantile_exp <- ncol(matrix_tf)/(1/cells_quantile)
    TfExp <- matrix_tf[rowSums(as.matrix(matrix_tf))>quantile_exp,]}
  return(TfExp)
}

#' filter original regulation
#' @description This function is used to further filter regulation generated by
#' Genin3, PIDC, or other methods.
#' @param potential_regulation data.frame consists of three columns, first column
#' contains source genes, second column contains target genes, third column contains
#' weight.
#' @param motif motif file, you can choose our bulit-in motif database of
#' 'mus musculus', 'homo sapiens', 'zebrafish' and 'chicken' by 'motif = Tranfac201803_Mm_MotifTFsF',
#' 'motif = Tranfac201803_Hs_MotifTFsF', 'motif = Tranfac201803_Zf_MotifTFsF',
#' 'motif = Tranfac201803_Ch_MotifTFsF' respectively, or you can upload your own motif data base, but the formata use be the same as our built-in motif database.
#' @param tf_threshold the numbers of targets with highest weight of each
#' transcription factor. Indicating the threshold to filter regulation.
#' @param target_threshold numbers of top transcription factors for each target.
#' Indicating the threshold to filter regulation.
#'
#' @return
#' @export
#'
#' @examples
filter_original_regulation <-function(potential_regulation,motif,
                                      tf_threshold = 500,target_threshold = 100){
  motifgene <- c()
  for (i in 1:nrow(motif)) {
    gene1 <- strsplit(motif[i,5],';')[[1]]
    motifgene <- c(motifgene,gene1)
  }
  potential_regulation <- potential_regulation[potential_regulation[,1]%in%motifgene,]
  TF <- as.character(potential_regulation[,1][!duplicated(potential_regulation[,1])])
  for (i in TF) {
    link_tf <- potential_regulation[potential_regulation[,1]==i,]
    if (nrow(link_tf)>tf_threshold) {
      TOP <- link_tf[order(link_tf[,3],decreasing = T),][1:tf_threshold,]
    }else{TOP <- link_tf[order(link_tf[,3],decreasing = T),]}
    if (i==TF[1]) {
      link_top <- TOP
    }else{link_top <- rbind(link_top,TOP)}
  }
  target <- as.character(link_top[,2][!duplicated(link_top[,2])])
  for (i in target) {
    link_target <- link_top[link_top[,2]==i,]
    if (nrow(link_target)>target_threshold) {
      TOP <- link_target[order(link_target[,3],decreasing = T),][1:target_threshold,]
    }else{TOP <- link_target}
    if (i==target[1]) {
      link_final <- TOP
    }else{link_final <- rbind(link_final,TOP)}
  }
  return(link_final)
}


#' Add regulation type for each gene pair
#' @description This function calculate the Spearman's correlation for each gene
#' pair as the regulation types. Correlation > 0 is positive regulation, Correlation
#' < 0 is negative regulation. Then, this function add the module and gene Symbol
#' for Kmeans_result to each gene in the gene pairs.
#'
#' @param Kmeans_result Kmeans result data.frame, row names should be ENSEMBEL ID,
#' and the first column should be gene Symbol ID, the second column should be KmeansGroup
#' @param potential_regulation First column is source gene, second column is target
#' gene, third column is weight. Several methods are recommended to generate
#' potential regulation: Genin3, PIDC, GRNBOOST2
#' @param start_column numeric, indicating the start column of expression value,
#' defalut is 3
#' @importFrom reshape2 melt
#' @return
#' @export
#'
#' @examples
add_regulation_type <- function(Kmeans_result,potential_regulation,start_column=4){
  source1 <- match(potential_regulation[,1],rownames(Kmeans_result))
  target1 <- match(potential_regulation[,2],rownames(Kmeans_result))
  source2 <- Kmeans_result[source1,c(1,2)]
  colnames(source2) <- c('TFSymbol','TFGroup')
  target2 <- Kmeans_result[target1,c(1,2)]
  colnames(target2) <- c('TargetSymbol','TargetGroup')
  colnames(potential_regulation) <- c('TF','Target','Weight')
  regulatory_relationships_Gen <- cbind(potential_regulation,source2,target2)
  regulatory_relationships_Gen <- regulatory_relationships_Gen[,c(1,4,5,2,6,7,3)]
  cor1 <- sparse.cor(t(Kmeans_result[,start_column:ncol(Kmeans_result)]))
  cor2 <- reshape2::melt(cor1)
  correlationIndex <-  match(paste(regulatory_relationships_Gen[,1],regulatory_relationships_Gen[,4]),
                             paste(cor2[,1],cor2[,2]))
  regulatory_relationships_Gen$Correlation <- cor2[correlationIndex,3]
  return(regulatory_relationships_Gen)
}
