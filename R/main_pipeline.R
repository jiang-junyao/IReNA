
#' irena_step1
#'
#' @param obj seurat object, where meta.data should contain Pseudotime column
#' @param metacell_num numeric, indicating excepted metacell number
#' @param module_number numeric,indicating module number
#'
#' @return return a data.frame, and the first column is K-means group
#' @export
#'
#' @examples
irena_step1 <- function(obj,metacell_num = 50,module_number=4){
  expression_profile <- get_SmoothByBin_PseudotimeExp(obj, Bin = metacell_num,
                                                      FC = F)
  expression_profile = filter_expression_profile(expression_profile,
                                                 filterfc = F)
  clustering <- clustering_Kmeans(expression_profile, K1=module_number)
  return(clustering)
}

#' Title
#'
#' @param exp
#' @param grn_method
#' @param start_column
#' @param nCores
#' @param candidate_tfs
#'
#' @return
#' @export
#'
#' @examples
irena_step2 <- function(exp,grn_method = 'pearson correlation',start_column=2,
                        nCores=5,candidate_tfs){
  if (grn_method == 'pearson correlation') {
    print('Using pearson correlation to infer basic grn')
    cor1 <- sparse.cor(t(exp[,start_column:ncol(exp)]))
    cor2 <- reshape2::melt(cor1)
    cor2 <- cor2[cor2[,3]!=1,]
    cor2 <- cor2[cor2$Var1 %in% candidate_tfs,]
    colnames(cor2) <- c('TF','Target','Correlation')
    SourceIdx <- match(cor2[,1],rownames(exp))
    TargetIdx <- match(cor2[,2],rownames(exp))
    cor2$TFGroup <- exp[SourceIdx,"KmeansGroup"]
    cor2$TargetGroup <- exp[TargetIdx,"KmeansGroup"]
    return(cor2)
  }else if (grn_method == 'GENIE3'){
    weightMat <- GENIE3::GENIE3(as.matrix(exp),nCores = nCores)
    weightMat <- getLinkList(weightMat)
    ### add regulation type for each gene pair
    regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS,regulation)
    ### check whether source genes are transcription factors
    motifTF <- c()
    for (i in 1:nrow(motif1)) {
      TF <- strsplit(motif1[i,5],';')[[1]]
      motifTF <- c(motifTF,TF)
    }
    regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]
  }
  return(regulatory_relationships)
}
