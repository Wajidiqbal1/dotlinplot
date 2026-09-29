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

test_that("bars are scaled to the maximum expression of all genes", {
  d <- detection_data(dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type"))

  # gene1 reaches 4 and gene2 reaches 6: both use 6.
  expect_equal(d$bar_ymin, rep(-0.15 * 6, 6))
  expect_equal(d$bar_ymax, rep(-0.05 * 6, 6))

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

test_that("the legend follows the category order", {
  # Only category B has enough expressing cells for a violin.
  p <- dotlin_plot(toy_data(), "gene1", "cell_type")
  fill <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")

  expect_equal(fill$get_limits(), c("A", "B", "C"))
})

test_that("x-axis breaks and labels match the categories", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type",
                   category_order = c("C", "A"))
  x_scale <- ggplot2::layer_scales(p)$x

  expect_equal(x_scale$get_breaks(), c(1, 2))
  expect_equal(x_scale$get_labels(), c("C", "A"))
})

test_that("category_order sets the order and filters categories", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type",
                   category_order = c("C", "A"))
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

test_that("y-axis ticks are as fine as on an ordinary plot of the data", {
  cells <- toy_data()
  ticks <- function(p) {
    y <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y
    breaks <- y$get_breaks()
    breaks[!is.na(breaks)]
  }
  ordinary <- ggplot2::ggplot(cells, ggplot2::aes(cell_type, gene1)) +
    ggplot2::geom_point()

  expect_equal(ticks(dotlin_plot(cells, "gene1", "cell_type")),
               ticks(ordinary))
})

test_that("percentages hang below their bars, inside the frame", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type")
  is_text <- vapply(p$layers, function(l) inherits(l$geom, "GeomText"),
                    logical(1))
  labels <- ggplot2::layer_data(p, which(is_text))
  d <- detection_data(p)
  y_range <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$
    continuous_range

  expect_equal(sort(labels$y), sort(d$bar_ymin))
  expect_gt(p$layers[[which(is_text)]]$aes_params$vjust, 1)
  # Room for the percentages inside the panel, above the x-axis line
  expect_lt(y_range[1], min(d$label_bottom))
  expect_s3_class(ggplot2::calc_element("axis.line.x.bottom", p$theme),
                  "element_line")
})

test_that("a single gene is named in a boxed strip above its panel", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type")

  expect_s3_class(p$facet, "FacetWrap")
  expect_false(inherits(p$theme$strip.background, "element_blank"))
})

test_that("the gene box lines up with the y-axis line", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type")
  axis <- ggplot2::calc_element("axis.line.y.left", p$theme)
  box <- ggplot2::calc_element("strip.background", p$theme)

  # Same width and not clipped, so both are centred on the panel edge.
  expect_equal(box$linewidth, axis$linewidth)
  expect_equal(p$theme$strip.clip, "off")
})

test_that("several genes are named to the left of their panels, unboxed", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type")

  expect_s3_class(p$facet, "FacetGrid")
  expect_equal(p$theme$strip.placement, "outside")
  expect_equal(p$theme$strip.text.y.left$angle, 0)
  expect_s3_class(p$theme$strip.background, "element_blank")
})

test_that("category names are angled only when they are long", {
  short <- dotlin_plot(toy_data(), "gene1", "cell_type")
  cells <- toy_data()
  cells$cell_type <- paste(cells$cell_type, "cells")
  long <- dotlin_plot(cells, "gene1", "cell_type")

  expect_equal(long$theme$axis.text.x$angle, 45)
  expect_null(short$theme$axis.text.x$angle)
})

test_that("a single gene has the original geometry below zero", {
  # gene1 reaches 4: bar from -0.2 to -0.6, percentage centred at -0.9, and
  # the panel ends at the usual 5% margin below it.
  p <- dotlin_plot(toy_data(), "gene1", "cell_type")
  is_text <- vapply(p$layers, function(l) inherits(l$geom, "GeomText"),
                    logical(1))
  labels <- ggplot2::layer_data(p, which(is_text))
  y_range <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$
    continuous_range

  expect_equal(unique(labels$y), -0.9)
  expect_equal(unique(labels$vjust), 0.5)
  expect_equal(y_range[1], -0.9 - 0.05 * (4 + 0.9))
})

test_that("all gene panels share the same y-axis", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type")
  panels <- ggplot2::ggplot_build(p)$layout$panel_params
  y_range <- lapply(panels, function(panel) panel$y$continuous_range)
  y_breaks <- lapply(panels, function(panel) panel$y$get_breaks())

  expect_length(panels, 2)
  expect_equal(y_range[[1]], y_range[[2]])
  expect_equal(y_breaks[[1]], y_breaks[[2]])
})

test_that("shared_y = FALSE gives each gene its own y-axis and bars", {
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type",
                   shared_y = FALSE)
  panels <- ggplot2::ggplot_build(p)$layout$panel_params
  y_top <- vapply(panels, function(panel) panel$y$continuous_range[2],
                  numeric(1))
  d <- detection_data(p)

  # gene1 reaches 4 and gene2 reaches 6.
  expect_lt(y_top[1], y_top[2])
  expect_equal(d$bar_ymin, rep(-0.15 * c(4, 6), each = 3))
})

test_that("violins have the same width whatever the share of expression", {
  # With min_nonzero = 4, violins are drawn for shares of 40% (gene1 A),
  # 75% (gene1 B) and 100% (gene2 A and C).
  p <- dotlin_plot(toy_data(), c("gene1", "gene2"), "cell_type",
                   min_nonzero = 4)
  is_violin <- vapply(p$layers, function(l) inherits(l$geom, "GeomViolin"),
                      logical(1))
  violins <- ggplot2::layer_data(p, which(is_violin))
  widths <- tapply(violins$xmax - violins$xmin,
                   paste(violins$PANEL, violins$x), max)

  expect_equal(as.vector(widths), rep(0.9, 4))
})

test_that("points have the default size of 0.7, as in the original", {
  p <- dotlin_plot(toy_data(), "gene1", "cell_type")

  expect_true(all(ggplot2::layer_data(p, 2)$size == 0.7))
})

test_that("labels never hide a few expressing or non-expressing cells", {
  cells <- data.frame(
    group = rep(c("one", "most", "none", "all"), each = 1000),
    gene = c(1, rep(0, 999), rep(1, 999), 0, rep(0, 1000), rep(1, 1000))
  )

  d <- detection_data(dotlin_plot(cells, "gene", "group",
                                  category_order = c("one", "most", "none",
                                                     "all")))

  expect_equal(d$label, c("<1%", ">99%", "0%", "100%"))
})

test_that("negative and missing values give warnings", {
  cells <- toy_data()
  cells$gene1[1] <- -1
  expect_warning(dotlin_plot(cells, "gene1", "cell_type"), "negative values")

  cells <- toy_data()
  cells$gene1[1] <- NA
  expect_warning(p <- dotlin_plot(cells, "gene1", "cell_type"),
                 "treated as zero")
  expect_equal(detection_data(p)$label, c("40%", "70%", "20%"))
})

test_that("invalid arguments give informative errors", {
  cells <- toy_data()

  expect_error(dotlin_plot(cells, "gene1", "cell_type", min_nonzero = 0),
               "positive integer")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", min_nonzero = 2.5),
               "positive integer")
  expect_error(dotlin_plot(cells, "gene1", "cell_type", shared_y = NA),
               "TRUE or FALSE")
  expect_error(dotlin_plot(cells, character(), "cell_type"),
               "at least one gene")
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
