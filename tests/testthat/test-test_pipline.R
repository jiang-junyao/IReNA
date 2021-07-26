test_that("merge footprints & get fasta", {
  merged_footprint<-merge_footprints(fdr005)
  fastadir<-'D:\\GIBH\\IReNA2 R package\\IReNA2\\Public\\GRCm38Chr\\Genome\\GRCm38Chr.fasta'
  fasta<-getfasta(merged_footprint,fastadir)
  expect_equal(nrow(fasta),17624)
  expect_equal(nrow(merged_footprint),8812)
})

test_that("overlap & annotate", {
  peak_bed <- get_bed(peak)
  overlapped <- overlap_footprints_peaks(combined, peak_bed)
  expect_equal(nrow(overlapped),9)
  expect_equal(nrow(overlapped),500)
  expect_equal(ncol(overlapped),3)
  motif1 <- Tranfac201803_Mm_MotifTFsF
  library(TxDb.Mmusculus.UCSC.mm10.knownGene )
  txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene
  list1<-get_related_genes(overlapped,motif=motif1,txdb = txdb,Species = 'Mm')
  expect_equal(list1[[2]][["V1"]],c(rep('Intron',4),'Distal Intergenic','Intron','Distal Intergenic','Intron','Promoter'))
  peak_genes<-get_peaks_genes(list1,MmscRNA_PHx_Exp_NewF)
})
