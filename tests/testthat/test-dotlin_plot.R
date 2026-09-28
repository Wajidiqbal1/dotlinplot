test_that("dotlin_plot returns a ggplot that renders without warnings", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type")

  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("non-zero proportions and labels are correct", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type")
  d <- detection_data(p)

  expect_equal(as.character(d$gene), rep(c("gene1", "gene2"), each = 3))
  expect_equal(as.character(d$category), rep(c("A", "B", "C"), 2))
  expect_equal(d$nonzero_count, c(4, 15, 1, 10, 0, 5))
  expect_equal(d$nonzero_prop, c(0.4, 0.75, 0.2, 1, 0, 1))
  expect_equal(d$label, c("40%", "75%", "20%", "100%", "0%", "100%"))
})

test_that("bars are scaled to the maximum expression of each gene", {
  d <- detection_data(dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type"))

  max_expression <- rep(c(4, 6), each = 3)
  expect_equal(d$bar_ymin, -0.15 * max_expression)
  expect_equal(d$bar_ymax, -0.05 * max_expression)
  expect_equal(d$label_y, -0.225 * max_expression)

  expect_equal(d$x, rep(1:3, 2))
  expect_equal(d$bar_xmin, d$x - 0.325)
  expect_equal(d$bar_xmax, d$x + 0.325)
  expect_equal(d$detection_xmax, d$bar_xmin + 0.65 * d$nonzero_prop)
})

test_that("bars of genes without expression are still drawn", {
  cells <- toy_data()
  cells$gene3 <- 0

  d <- detection_data(dotlin_plot(cells, "gene3", "cell_type"))

  expect_equal(d$label, rep("0%", 3))
  expect_true(all(d$bar_ymin < d$bar_ymax))
})

test_that("points show all non-zero values, violins need min_nonzero", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type")

  points <- geom_data(p, "GeomPoint")
  expect_equal(nrow(points), 35)
  expect_true(all(points$expression > 0))

  violins <- geom_data(p, "GeomViolin")
  expect_setequal(unique(paste(violins$gene, violins$category)),
                  c("gene1 B", "gene2 A"))

  violins <- geom_data(dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type",
                                 min_nonzero = 5), "GeomViolin")
  expect_setequal(unique(paste(violins$gene, violins$category)),
                  c("gene1 B", "gene2 A", "gene2 C"))
})

test_that("a single non-zero value does not produce a violin", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type", min_nonzero = 1)

  violins <- geom_data(p, "GeomViolin")
  expect_setequal(unique(as.character(violins$category)), c("A", "B"))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("points are jittered horizontally only", {
  set.seed(1)
  p <- dotlin_plot(toy_data(), "gene1", "cell_type")
  points <- ggplot2::layer_data(p, 2)

  expect_equal(sort(points$y), sort(geom_data(p, "GeomPoint")$expression))
  expect_true(all(abs(points$x - round(points$x)) <= 0.1))
})

test_that("genes keep the order given and duplicates are dropped", {
  p <- dotlin_plot(toy_data(), c("gene2", "gene1", "gene2"), "cell_type")

  expect_equal(levels(detection_data(p)$gene), c("gene2", "gene1"))
  expect_equal(nrow(detection_data(p)), 6)
})

test_that("categories follow factor levels, or sorted values otherwise", {
  cells <- toy_data()
  p <- dotlin_plot(cells, "gene1", "cell_type")
  expect_equal(levels(detection_data(p)$category), c("A", "B", "C"))

  cells$cell_type <- factor(cells$cell_type, c("C", "A", "B", "unused"))
  p <- dotlin_plot(cells, "gene1", "cell_type")
  expect_equal(levels(detection_data(p)$category), c("C", "A", "B"))

  cells$cell_type <- rep(c(10, 2, 1), times = c(20, 10, 5))
  p <- dotlin_plot(cells, "gene1", "cell_type")
  expect_equal(levels(detection_data(p)$category), c("1", "2", "10"))
})

test_that("x-axis breaks and labels match the categories", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type", category_order = c("C", "A"))
  x_scale <- ggplot2::layer_scales(p)$x

  expect_equal(x_scale$get_breaks(), c(1, 2))
  expect_equal(x_scale$get_labels(), c("C", "A"))
})

test_that("category_order sets the order and filters categories", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type", category_order = c("C", "A"))
  d <- detection_data(p)

  expect_equal(levels(d$category), c("C", "A"))
  expect_equal(d$label, c("20%", "40%"))
  expect_setequal(unique(as.character(geom_data(p, "GeomPoint")$category)),
                  c("A", "C"))
})

test_that("a category named \"\" is handled", {
  cells <- toy_data()
  cells$cell_type[cells$cell_type == "C"] <- ""

  p <- dotlin_plot(cells, "gene1", "cell_type",
                 palette = c("red", "green", "blue"))
  d <- detection_data(p)
  expect_equal(d$label, c("20%", "40%", "75%"))
  expect_equal(ggplot2::layer_data(p, 4)$fill, c("red", "green", "blue"))

  p <- dotlin_plot(cells, "gene1", "cell_type",
                 palette = c(A = "green", B = "blue", "red"))
  expect_equal(ggplot2::layer_data(p, 4)$fill, c("red", "green", "blue"))
})

test_that("cells without a category are left out", {
  cells <- toy_data()
  cells$cell_type[cells$cell_type == "C"] <- NA

  d <- detection_data(dotlin_plot(cells, "gene1", "cell_type"))

  expect_equal(levels(d$category), c("A", "B"))
})

test_that("palettes can be named, unnamed or generated", {
  default <- dotlin_plot(toy_data(), "gene1", "cell_type")
  expect_equal(ggplot2::layer_data(default, 4)$fill,
               scales::hue_pal()(3))

  named <- c(C = "red", B = "green", A = "blue", D = "black")
  p <- dotlin_plot(toy_data(), "gene1", "cell_type", palette = named)
  expect_equal(ggplot2::layer_data(p, 4)$fill, c("blue", "green", "red"))

  p <- dotlin_plot(toy_data(), "gene1", "cell_type",
                 palette = list(A = "blue", B = "green", C = "red"))
  expect_equal(ggplot2::layer_data(p, 4)$fill, c("blue", "green", "red"))

  p <- dotlin_plot(toy_data(), "gene1", "cell_type",
                 palette = c("blue", "green", "red", "black"))
  expect_equal(ggplot2::layer_data(p, 4)$fill, c("blue", "green", "red"))
})

test_that("labels and title are set", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type", title = "Dotlin plot")

  expect_equal(p$labels$title, "Dotlin plot")
  expect_equal(p$labels$x, "cell_type")
  expect_equal(p$labels$fill, "cell_type")
  expect_equal(p$labels$y, "Expression")
})

test_that("y-axis breaks are non-negative", {
  expect_equal(nonnegative_breaks(c(-0.9, 4.2)), c(0, 1, 2, 3, 4))
})

test_that("negative and missing values give warnings", {
  cells <- toy_data()
  cells$gene1[1] <- -1
  expect_warning(dotlin_plot(cells, "gene1", "cell_type"), "negative values")

  cells <- toy_data()
  cells$gene1[1] <- NA
  expect_warning(p <- dotlin_plot(cells, "gene1", "cell_type"), "treated as zero")
  expect_equal(detection_data(p)$label, c("40%", "70%", "20%"))
})

test_that("invalid arguments give informative errors", {
  cells <- toy_data()

  expect_error(dotlin_plot(cells, "gene1", "cell_type", min_nonzero = 0),
               "positive integer")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", min_nonzero = 2.5),
               "positive integer")
  expect_error(dotlin_plot(cells, character(), "cell_type"), "at least one gene")
  expect_error(dotlin_plot(cells, c("gene1", "nope"), "cell_type"), "'nope'")
  expect_error(dotlin_plot(cells, "gene1"), "must be supplied")
  expect_error(dotlin_plot(cells, "gene1", "nope"), "not a column")
  expect_error(dotlin_plot(cells, "gene1", c("a", "b")), "single column name")
  expect_error(dotlin_plot(cells, "cell_type", "cell_type"), "not numeric")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", layer = "counts"),
               "do not apply")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", layer = c("a", "b")),
               "single layer name")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", assay = c("a", "b")),
               "single assay name")
  expect_error(dotlin_plot(cells, "gene1", "cell_type",
                         category_order = character()),
               "cannot be empty")
  expect_error(dotlin_plot(cells, "gene1", "cell_type",
                         category_order = c("A", "Z")),
               "'Z'")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", palette = c(A = "red")),
               "missing colours .*'B', 'C'")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", palette = "red"),
               "1 colours, but 3")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", palette = 1:3),
               "character vector")
  expect_error(dotlin_plot(as.matrix(cells), "gene1", "cell_type"),
               "must be a Seurat object")
})
