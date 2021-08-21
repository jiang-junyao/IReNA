test_that("merge footprints & get fasta", {
  merged_footprint <- merge_footprints(fdr005)
  fastadir <- "D:\\GIBH\\IReNA2 R package\\IReNA2\\Public\\GRCm38Chr\\Genome\\GRCm38Chr.fasta"
  fasta <- getfasta(merged_footprint, fastadir)
  expect_equal(nrow(fasta), 1770)
  expect_equal(nrow(merged_footprint), 885)
})

test_that("overlap & annotate", {
  peak_bed <- get_bed(peak)
  overlapped <- overlap_footprints_peaks(combined, test_peak)
  expect_equal(nrow(overlapped), 7)
  expect_equal(ncol(overlapped), 10)
})

test_that("FOS", {
  Kmeans_clustering_ENS <- add_ENSID(test_clustering, Spec1='Hs')
  regulatory_relationships <- Footprints_FOS(wig_list, Candid,0.1)
  expect_equal(nrow(regulatory_relationships), 538)
})
