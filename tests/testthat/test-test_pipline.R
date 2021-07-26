test_that("merge footprints & get fasta", {
  merged_footprint <- merge_footprints(fdr005)
  fastadir <- "D:\\GIBH\\IReNA2 R package\\IReNA2\\Public\\GRCm38Chr\\Genome\\GRCm38Chr.fasta"
  fasta <- getfasta(merged_footprint, fastadir)
  expect_equal(nrow(fasta), 17624)
  expect_equal(nrow(merged_footprint), 8812)
})

test_that("overlap & annotate", {
  peak_bed <- get_bed(peak)
  overlapped <- overlap_footprints_peaks(combined, peak_bed)
  expect_equal(nrow(overlapped), 9)
  expect_equal(ncol(overlapped), 10)
})

test_that("FOS", {
  regulatory_relationships <- Footprints_FOS(wig_list, Candid, MmscRNA_PHx_Exp_NewF)
  expect_equal(nrow(regulatory_relationships), 191)
})
