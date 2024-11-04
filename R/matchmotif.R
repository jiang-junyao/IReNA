
Str_to_GR <- function(x){
  sp = strsplit(x,split='-')
  chr = sapply(sp,function(x) x[[1]])
  s = as.numeric(sapply(sp,function(x) x[[2]]))
  e = as.numeric(sapply(sp,function(x) x[[3]]))
  GR_out = GRanges(chr,IRanges(s,e))
  return(GR_out)
}

Seqnames <- function(x){
  out= as.character(seqnames(x))
  return(out)
}

Start <- function(x){
  out= as.numeric(start(x))
  return(out)
}

End <- function(x){
  out= as.numeric(end(x))
  return(out)
}

Must_to_GR <- function(x){
  library('pbapply')
  chr_all = pblapply(x,Seqnames)
  start_all = pblapply(x,Start)
  end_all = pblapply(x,End)
  len_all = pblapply(x,function(x) length(x))
  #####
  chr_all = as.character(unlist(chr_all))
  start_all = as.numeric(unlist(start_all))
  end_all = as.numeric(unlist(end_all))
  ##### #####
  names_all = rep(names(x),len_all)
  GR_out = GRanges(chr_all,IRanges(start_all,end_all),motifs=names_all)
  return(GR_out)
}

motifs_select2 <- function(motif,gene){
  index <- c()
  if (stringr::str_sub(gene[1],1,3)=='ENS') {
    col_idx = 5
  }else{col_idx = 4}
  for (i in 1:nrow(motif)) {
    judge <- c()
    gene1 <- strsplit(motif[i,col_idx],';')[[1]]
    for (j in gene1) {
      if (j %in% gene) {
        judge <- c(judge,'YSE')
      }
    }
    if ('YSE' %in% judge) {
      index <- c(index,i)
    }
  }
  motif1 <- motif[index,]
  return(motif1)
}

#' Title
#'
#' @param GR
#' @param gene.use
#' @param motifdb
#' @param pvalue.cutoff
#' @param BSdb
#'
#' @return
#' @export
#'
#' @examples
identify_region_tfs <- function(GR,gene.use,motifdb,
                                pvalue.cutoff = 5e-05,BSdb){
  motif_use = motifs_select2(motifdb,gene.use)
  PWM = Transfac_PWMatrixList
  PWM = PWM[motif_use$Accession]
  matched_motif <- motifmatchr::matchMotifs(PWM,
                               GR,genome = BSdb,
                               out='positions',p.cutoff = pvalue.cutoff)
  matched_motif <- Must_to_GR(matched_motif)
  overlapped_region <- findOverlaps(GR,matched_motif)
  enriched_tf <- c()
  if (stringr::str_sub(gene.use[1],1,3)=='ENS') {
    col_idx = 5
  }else{col_idx = 4}
  for (i in unique(overlapped_region@from)) {
    all_motif <- matched_motif$motifs[overlapped_region[overlapped_region@from==i]@to]
    all_motif <- motifdb[motifdb$Accession%in%all_motif,]
    all_tf <- paste(unique(unlist(strsplit(all_motif[,col_idx],';')))
                    ,collapse = ';')
    names(all_tf) <- gene.use[i]
    enriched_tf <- c(enriched_tf,all_tf)
  }
  regulation <- data.frame('TF'=enriched_tf,'Target'=names(enriched_tf))
  return(regulation)
}

overlap_peak_motif <- function(peak,motif,motifdb){
  overlaped = findOverlaps(peak,motif)
  peak_motif = cbind(as.data.frame(peak[overlaped@from]),as.data.frame(motif[overlaped@to]))
  peak_motif$TF = motifdb[match(peak_motif$motifs,motifdb$Accession),4]
  return(peak_motif)
}

make_tf_target <- function(atac_out){
  tf = atac_out$TF
  tf = paste0(tf,'#',atac_out$symbol)
  tf_target = unlist(map(tf,~paste_gene(.x)))
  return(tf_target)
}

paste_gene <- function(gene){
  tf = strsplit(gene,'#')[[1]][1]
  target = strsplit(gene,'#')[[1]][2]
  tf = unlist(strsplit(tf,';'))
  return(paste0(tf,'-',target))
}
