library(R6)
library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

#' RegionPlot R6 Class
#'
#' This class creates a genomic region plot using various tracks from the Gviz package.
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{new(bgFile, chr, from, to)}}{Initialize the RegionPlot object.}
#'   \item{\code{createTracks()}}{Create the necessary tracks for visualization.}
#'   \item{\code{plot(output_file = "~/Downloads/regionPlot.pdf")}}{Plot the tracks to a PDF file.}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{bgFile}}{The path to the background file.}
#'   \item{\code{chr}}{Chromosome name.}
#'   \item{\code{from}}{Start position of the region.}
#'   \item{\code{to}}{End position of the region.}
#' }
#'
#' @examples
#' \dontrun{
#' plotter <- RegionPlot$new("path_to_my_file", "chr17", 7673300, 7675000)
#' plotter$plot()
#' }

RegionPlot <- R6Class("RegionPlot",
                      public = list(
                        bgFile = NULL,
                        chr = NULL,
                        from = NULL,
                        to = NULL,
                        alignmentsTrack = NULL,
                        seqTrack = NULL,
                        ideoTrack = NULL,
                        axisTrack = NULL,
                        geneTrack = NULL,
                        txTrack = NULL,
                        txTrackRef = NULL,

                        #' Initialize the RegionPlot object
                        #'
                        #' @param bgFile The path to the background file.
                        #' @param chr Chromosome name.
                        #' @param from Start position of the region.
                        #' @param to End position of the region.
                        initialize = function(bgFile, chr, from, to) {
                          self$bgFile <- bgFile
                          self$chr <- chr
                          self$from <- from
                          self$to <- to
                          self$createTracks()
                        },

                        #' Create the necessary tracks for visualization
                        createTracks = function() {
                          # Create AlignmentsTrack
                          self$alignmentsTrack <- AlignmentsTrack(
                            start = self$from,
                            end = self$to,
                            range = self$bgFile,
                            chromosome = self$chr,
                            showIndels = TRUE
                          )

                          # Create SequenceTrack
                          self$seqTrack <- SequenceTrack(
                            Hsapiens,
                            chromosome = self$chr,
                            cex = 0.5,
                            min.height = 8
                          )

                          # Create IdeogramTrack
                          self$ideoTrack <- IdeogramTrack(
                            genome = "hg38",
                            chromosome = self$chr,
                            showId = FALSE
                          )

                          # Create GenomeAxisTrack
                          self$axisTrack <- GenomeAxisTrack()

                          # Extract gene data from the TxDb object
                          txdb <- TxDb.Hsapiens.UCSC.hg38.refGene
                          genes <- genes(txdb)

                          # Map gene IDs to symbols using the annotation package
                          geneSymbols <- mapIds(
                            org.Hs.eg.db,
                            keys = genes$gene_id,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first"
                          )
                          # Add gene symbols to the GRanges object
                          mcols(genes)$symbol <- geneSymbols

                          # Create GeneRegionTrack with gene symbols
                          self$geneTrack <- GeneRegionTrack(
                            genes,
                            chromosome = self$chr,
                            start = self$from,
                            end = self$to,
                            gene = geneSymbols,
                            collapseTranscripts = "gene",
                            symbol = geneSymbols,
                            fill = "darkgreen",
                            col = "darkgreen",
                            shape = "arrow"
                          )

                          # Create GeneRegionTrack for knownGene
                          self$txTrack <- GeneRegionTrack(
                            TxDb.Hsapiens.UCSC.hg38.knownGene,
                            chromosome = self$chr,
                            start = self$from,
                            end = self$to,
                            fill = "purple4",
                            col = "purple4"
                          )

                          # Create GeneRegionTrack for refGene
                          self$txTrackRef <- GeneRegionTrack(
                            TxDb.Hsapiens.UCSC.hg38.refGene,
                            chromosome = self$chr,
                            start = self$from,
                            end = self$to,
                            col = "darkorange",
                            fill = "darkorange"
                          )

                          # Set display parameters for alignments track
                          displayPars(self$alignmentsTrack) <- list(fontsize = 15, showTitle = FALSE)
                        },

                        #' Plot the tracks to a PDF file
                        #'
                        #' @param output_file The path to the output PDF file. Default is "~/Downloads/regionPlot.pdf".
                        plot = function(output_file = "~/Downloads/regionPlot.pdf") {
                          pdf(output_file, width = 20, height = 10)
                          plotTracks(
                            list(
                              self$ideoTrack,
                              self$axisTrack,
                              self$alignmentsTrack,
                              self$seqTrack,
                              self$txTrack,
                              self$txTrackRef,
                              self$geneTrack
                            ),
                            from = self$from,
                            to = self$to,
                            transcriptAnnotation = "transcript",
                            exponent = 6,
                            cex = 0.5
                          )
                          dev.off()
                        }
                      )
)

# Example usage:
# plotter <- RegionPlot$new("path_to_my_file", "chr17", 7673300, 7675000)
# plotter$plot()
