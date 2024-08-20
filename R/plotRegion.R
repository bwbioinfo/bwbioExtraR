#' RegionPlot R6 Class
#'
#' The `RegionPlot` class provides functionality for creating and visualizing a genomic region plot using various tracks from the Gviz package. This class is particularly useful for visualizing specific regions of the human genome (hg38).
#'
#' @importFrom R6 R6Class
#' @importFrom Gviz AlignmentsTrack
#' @importFrom Gviz SequenceTrack
#' @importFrom Gviz IdeogramTrack
#' @importFrom Gviz GenomeAxisTrack
#' @importFrom Gviz GeneRegionTrack
#' @importFrom Gviz plotTracks
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg38.refGene TxDb.Hsapiens.UCSC.hg38.refGene
#' @importFrom GenomicRanges mcols
#'
#' @examples
#' \dontrun{
#' plotter <- RegionPlot$new("path_to_my_file.bam", "chr17", 7673300, 7675000)
#' plotter$plot(output_file = "regionPlot.pdf")
#' }
#'
#' @name RegionPlot

RegionPlot <- R6Class("RegionPlot",
                      public = list(
                        #' @field bgFile The path to the background file.
                        bgFile = NULL,
                        #' @field chr Chromosome name.
                        chr = NULL,
                        #' @field from Start position of the region.
                        from = NULL,
                        #' @field to End position of the region.
                        to = NULL,
                        #' @field alignmentsTrack AlignmentsTrack object.
                        alignmentsTrack = NULL,
                        #' @field seqTrack SequenceTrack object.
                        seqTrack = NULL,
                        #' @field ideoTrack IdeogramTrack object.
                        ideoTrack = NULL,
                        #' @field axisTrack GenomeAxisTrack object.
                        axisTrack = NULL,
                        #' @field geneTrack GeneRegionTrack object.
                        geneTrack = NULL,
                        #' @field txTrack GeneRegionTrack object for knownGene.
                        txTrack = NULL,
                        #' @field txTrackRef GeneRegionTrack object for refGene.
                        txTrackRef = NULL,
                        #' Initialize the RegionPlot object
                        #' @description
                        #' This method initializes the RegionPlot object with the
                        #' specified background file, chromosome, and genomic region.
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
                        #' @description
                        #' This method creates the AlignmentsTrack, SequenceTrack,
                        #' IdeogramTrack, GenomeAxisTrack, GeneRegionTrack, and
                        #' GeneRegionTrack objects for visualization.
                        #' @return NULL
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
                        #' @description
                        #' This method plots the tracks to a PDF file for visualization.
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
