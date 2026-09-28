expected_labels <- c("40%", "75%", "20%", "100%", "0%", "100%")

test_that("Seurat objects use the identities and the data layer by default", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  object <- toy_seurat()

  p <- dotlin_plot(object, c("gene1", "gene2"))
  d <- detection_data(p)

  expect_equal(levels(d$category), c("C", "B", "A"))
  expect_equal(p$labels$x, "Identity")
  expect_equal(d$label[order(d$gene, match(d$category, c("A", "B", "C")))],
               expected_labels)
  expect_equal(max(geom_data(p, "GeomPoint")$expression), log1p(30))
})

test_that("Seurat objects accept a metadata column, layer and assay", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  object <- toy_seurat()

  p <- dotlin_plot(object, c("gene1", "gene2"), "cell_type")
  expect_equal(detection_data(p)$label, expected_labels)
  expect_equal(p$labels$x, "cell_type")

  p <- dotlin_plot(object, "gene2", "cell_type", layer = "counts",
                   assay = "RNA")
  expect_equal(max(geom_data(p, "GeomPoint")$expression), 30)
})

test_that("Seurat layers split by sample are combined", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  object <- toy_seurat()
  object[["RNA"]] <- split(object[["RNA"]], f = object$cell_type)

  d <- detection_data(dotlin_plot(object, c("gene1", "gene2"), "cell_type"))

  expect_equal(d$label, expected_labels)
})

test_that("Seurat objects with a v3 assay are supported", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  object <- toy_seurat()
  counts <- SeuratObject::LayerData(object, layer = "counts")
  assay <- SeuratObject::CreateAssayObject(counts = counts)
  assay <- SeuratObject::SetAssayData(assay, layer = "data",
                                      new.data = log1p(counts))
  SeuratObject::Key(assay) <- "old_"
  object[["OLD"]] <- assay

  d <- detection_data(dotlin_plot(object, c("gene1", "gene2"), "cell_type",
                                  assay = "OLD"))

  expect_equal(d$label, expected_labels)
})

test_that("Seurat genes are not confused with metadata columns", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  object <- toy_seurat()
  object$gene1 <- 100

  d <- detection_data(dotlin_plot(object, "gene1", "cell_type"))

  expect_equal(d$label, c("40%", "75%", "20%"))
})

test_that("Seurat objects give informative errors", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  object <- toy_seurat()

  expect_error(dotlin_plot(object, "gene1", assay = "ADT"), "Assay 'ADT'")
  expect_error(dotlin_plot(object, "gene1", layer = "dat"), "Layer 'dat'")
  expect_error(dotlin_plot(object, "gene1", layer = "scale.data"),
               "Available layers: 'counts', 'data'")
  expect_error(dotlin_plot(object, c("gene1", "nope")), "'nope'")
  expect_error(dotlin_plot(object, "gene1", "nope"), "cell metadata")
})

test_that("SingleCellExperiment objects default to colLabels and logcounts", {
  skip_if_not_installed("SingleCellExperiment")
  object <- toy_sce()

  p <- dotlin_plot(object, c("gene1", "gene2"))
  expect_equal(detection_data(p)$label, expected_labels)
  expect_equal(p$labels$x, "label")
  expect_equal(max(geom_data(p, "GeomPoint")$expression), log1p(30))

  p <- dotlin_plot(object, "gene2", "cell_type", layer = "counts")
  expect_equal(max(geom_data(p, "GeomPoint")$expression), 30)
  expect_equal(p$labels$x, "cell_type")
})

test_that("SingleCellExperiment objects give informative errors", {
  skip_if_not_installed("SingleCellExperiment")

  expect_error(dotlin_plot(toy_sce(label = FALSE), "gene1"), "must be supplied")
  expect_error(dotlin_plot(toy_sce(), "gene1", layer = "normcounts"),
               "Available assays: 'counts', 'logcounts'")
  expect_error(dotlin_plot(toy_sce(), "gene1", assay = "counts"),
               "choose the assay with `layer`")
  expect_error(dotlin_plot(toy_sce(), "nope"), "'nope'")
  expect_error(dotlin_plot(toy_sce(), "gene1", "nope"), "colData")
})
