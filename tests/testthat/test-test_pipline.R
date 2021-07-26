test_that("merge footprints & get fasta", {
  merged_footprint<-merge_footprints(fdr005)
  fastadir<-'D:\\GIBH\\IReNA2 R package\\IReNA2\\Public\\GRCm38Chr\\Genome\\GRCm38Chr.fasta'
  fasta<-getfasta(merged_footprint,fastadir)
  expect_equal(nrow(fasta),17624)
  expect_equal(nrow(merged_footprint),8812)
})

test_that("merge footprints & get fasta", {
  merged_footprint<-merge_footprints(fdr005)
  fastadir<-'D:\\GIBH\\IReNA2 R package\\IReNA2\\Public\\GRCm38Chr\\Genome\\GRCm38Chr.fasta'
  fasta<-getfasta(merged_footprint,fastadir)
  expect_equal(nrow(fasta),17624)
  expect_equal(nrow(merged_footprint),8812)
})
