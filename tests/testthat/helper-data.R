# Toy data with known detection rates (number of cells with expression > 0):
#
#               gene1   gene2
#   A (10)          4      10
#   B (20)         15       0
#   C  (5)          1       5
#
# Categories appear in the order B, A, C. The maximum is 4 for gene1 and 6 for
# gene2.
toy_data <- function() {
  data.frame(
    cell_type = rep(c("B", "A", "C"), times = c(20, 10, 5)),
    gene1 = c(seq_len(15) / 5, rep(0, 5), 1:4, rep(0, 6), 2.5, rep(0, 4)),
    gene2 = c(rep(0, 20), seq(0.5, 5, by = 0.5), c(1, 2, 3, 4, 6))
  )
}

# The data of the first layer drawn with the given geom.
geom_data <- function(p, geom) {
  for (layer in p$layers) {
    if (inherits(layer$geom, geom)) {
      return(layer$data)
    }
  }
  stop("No layer with geom ", geom)
}

detection_data <- function(p) {
  d <- geom_data(p, "GeomRect")
  d[order(d$gene, d$category), ]
}

# A small Seurat object with log-normalised data and a "cell_type" column.
toy_seurat <- function() {
  counts <- rbind(t(as.matrix(toy_data()[c("gene1", "gene2")])) * 5, other = 1)
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  object <- SeuratObject::SetAssayData(object, layer = "data",
                                       new.data = log1p(counts))
  object$cell_type <- toy_data()$cell_type
  SeuratObject::Idents(object) <- factor(object$cell_type,
                                         levels = c("C", "B", "A"))
  object
}

toy_sce <- function(label = TRUE) {
  counts <- t(as.matrix(toy_data()[c("gene1", "gene2")])) * 5
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = log1p(counts))
  )
  object$cell_type <- toy_data()$cell_type

  if (label) {
    SingleCellExperiment::colLabels(object) <- object$cell_type
  }

  object
}
