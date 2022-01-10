
test_that("overlap & annotate", {
  load(system.file("extdata", "test_peak.rda", package = "IReNA"))
  load(system.file("extdata", "combined.rda", package = "IReNA"))
  peak_bed <- get_bed(test_peak)
  overlapped <- overlap_footprints_peaks(combined, peak_bed)
  expect_equal(nrow(overlapped), 7)
  expect_equal(ncol(overlapped), 10)
})

test_that("FOS", {
  load(system.file("extdata", "test_clustering.rda", package = "IReNA"))
  load(system.file("extdata", "wig_list.rda", package = "IReNA"))
  load(system.file("extdata", "Candid.rda", package = "IReNA"))
  Kmeans_clustering_ENS <- add_ENSID(test_clustering, Spec1='Hs')
  regulatory_relationships <- Footprints_FOS(wig_list, Candid,0.1)
  expect_equal(nrow(regulatory_relationships), 747)
})
