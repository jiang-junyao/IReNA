#' ID change
#' @description Function to transfer Symbol ID to ENSEMBEL or ENSEMBEL to Symbol
#' @param Gene1 vector, indicating original gene ids
#' @param GeneInf1 data.frame, correspondence file of gene ID, row names should be ENSEMBLE ID, first column should be Symbol ID, if GeneInf1 = NULL this function will be built-in corresponding gene ID file.
#' @param Spec1 If you don’t have a gene ID corresponding file, you can also use our built-in corresponding gene ID file, 'Mm' for mus musculus
#'
#' @return return a data.frame contain original gene ids and conversed gene ids
#' @export
#'
#' @examples load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
#'
#' Converse_GeneIDSymbol(rownames(test_clustering), Spec1 = 'Hs')
Converse_GeneIDSymbol <- function(Gene1, GeneInf1 = NULL, Spec1 = "") {
  if (is.null(GeneInf1)) {
    if (Spec1[1] == "") {
      stop("Error: please input species for variable Spec1")
    } else if (Spec1 == "Mm") {
      GeneInf1 <- MmscRNA_genes
    } else if (Spec1 == 'Zf') {
      GeneInf1 <- ZfscRNA_genes
    } else if (Spec1 == 'Ch') {
      GeneInf1 <- ChscRNA_genes
    } else if (Spec1 == 'Hs') {
      GeneInf1 <- HsscRNA_genes
    }
  }

  if (grepl("ENS", Gene1[1])) {
    GeneID01 <- Gene1
    Ind1 <- match(GeneID01, rownames(GeneInf1))
    if (length(Ind1[!is.na(Ind1)]) > 0) {
      GeneID1 <- GeneID01[!is.na(Ind1)]
      Symb1 <- as.character(GeneInf1[Ind1[!is.na(Ind1)], 3])
      if (length(GeneID01[is.na(Ind1)]) > 0) {
        message(paste("Warning: no", paste(GeneID01[is.na(Ind1)], collapse = " ")))
      }
    } else {
      stop("Error: no such gene")
    }
  } else {
    Symb01 <- Gene1
    Ind1 <- apply(as.matrix(Symb01), 1, function(x1) {
      if (is.element(x1, GeneInf1[, 1])) {
        x2 <- match(x1, GeneInf1[, 1])
      } else if (is.element(x1, GeneInf1[, 3])) {
        x2 <- match(x1, GeneInf1[, 3])
      } else {
        x2 <- NA
      }
      return(x2)
    })

    if (length(Ind1[!is.na(Ind1)]) > 0) {
      GeneID1 <- rownames(GeneInf1[Ind1[!is.na(Ind1)], ])
      Symb1 <- Symb01[!is.na(Ind1)]
      if (length(Symb01[is.na(Ind1)]) > 0) {
        message(paste("Warning: no", paste(Symb01[is.na(Ind1)], collapse = " ")))
      }
    } else {
      stop("Error: no such gene")
    }
  }
  GeneID2 <- cbind(GeneID1, Symb1)
  colnames(GeneID2) <- c("EnsemblID", "Symbol")

  return(GeneID2)
}
