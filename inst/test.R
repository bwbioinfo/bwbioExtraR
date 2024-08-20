library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

bgFile <- "path_to_my_file"

chr = "chr17"
from <- 7673300
to <- 7675000

alignmentsTrack <-
AlignmentsTrack(
  start = from,
  end = to,
  range = bgFile,
  chromosome = "chr17",
  showIndels = TRUE
)

seqTrack <-
SequenceTrack(
  Hsapiens,
  chromosome = "chr17",
  cex = 0.5,
  min.height = 8
)

ideoTrack <-
IdeogramTrack(
  genome = "hg38",
  chromosome = "chr17",
  showId = FALSE
)

axisTrack <- GenomeAxisTrack()

txdb <- TxDb.Hsapiens.UCSC.hg38.refGene
# Extract gene data from the TxDb object
genes <- genes(txdb)
# Map gene IDs to symbols using the annotation package
geneSymbols <-
  mapIds(
    org.Hs.eg.db,
    keys = genes$gene_id,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)
# Add gene symbols to the GRanges object
mcols(genes)$symbol <- geneSymbols
# Create the GeneRegionTrack with gene symbols
geneTrack <-
  GeneRegionTrack(
    genes,
    chromosome = chr,
    start = from,
    end = to,
    gene = geneSymbols,
    collapseTranscripts = "gene",
    symbol = geneSymbols,
    fill = "darkgreen",
    col = "darkgreen",
    shape = "arrow"
  )

txTrack <- GeneRegionTrack(
    TxDb.Hsapiens.UCSC.hg38.knownGene,
    chromosome = chr,
    start = from,
    end = to,
    fill = "purple4",
    col = "purple4"
  )

txTrackRef <- GeneRegionTrack(
    TxDb.Hsapiens.UCSC.hg38.refGene,
    chromosome = chr,
    start = from,
    end = to,
    col = "darkorange",
    fill = "darkorange"
  )

displayPars(alignmentsTrack) <- list(fontsize = 15, showTitle = FALSE )

pdf(
    "~/Downloads/regionPlot.pdf",
    width = 20,
    height = 10
  )

plotTracks(
    list(
      ideoTrack,
      axisTrack,
      alignmentsTrack,
      seqTrack,
      # knownTranscripts,
      txTrack,
      txTrackRef,
      geneTrack
    ),
    from = from, to = to,
    transcriptAnnotation = "transcript",
    exponent = 6, cex = 0.5
  )

dev.off()
