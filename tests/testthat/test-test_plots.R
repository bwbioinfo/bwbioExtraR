library(testthat)
library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(R6)

# Assuming RegionPlot is in your package namespace
# Ensure this script is located in tests/testthat/ directory

# Test Initialization of RegionPlot
test_that("RegionPlot initializes correctly", {
  plotter <-
    RegionPlot$new(
      system.file(
        package = "bwbioExtraR",
        "extdata",
        "Test.hg38-chr17.bam"
        ),
      "chr17",
      7673300,
      7675000
      )

  expect_true(is.R6(plotter))
  expect_true(inherits(plotter, "RegionPlot"))
  expect_equal(plotter$chr, "chr17")
  expect_equal(plotter$from, 7673300)
  expect_equal(plotter$to, 7675000)
})

# Test CreateTracks Method
test_that("createTracks method creates tracks correctly", {
  plotter <- RegionPlot$new(
    system.file(
      package = "bwbioExtraR",
      "extdata",
      "Test.hg38-chr17.bam"
    ),
    "chr17",
    7673300, 7675000)
  plotter$createTracks()

  expect_s4_class(plotter$alignmentsTrack, "AlignmentsTrack")
  expect_s4_class(plotter$seqTrack, "SequenceTrack")
  expect_s4_class(plotter$ideoTrack, "IdeogramTrack")
  expect_s4_class(plotter$axisTrack, "GenomeAxisTrack")
  expect_s4_class(plotter$geneTrack, "GeneRegionTrack")
  expect_s4_class(plotter$txTrack, "GeneRegionTrack")
  expect_s4_class(plotter$txTrackRef, "GeneRegionTrack")
})

# Test Plot Method
test_that("plot method generates a PDF", {
  plotter <- RegionPlot$new(
    system.file(
      package = "bwbioExtraR",
      "extdata",
      "Test.hg38-chr17.bam"
    ),
    "chr17",
    7673300,
    7675000
    )
  output_file <- tempfile()

  plotter$plot(output_file)

  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})
