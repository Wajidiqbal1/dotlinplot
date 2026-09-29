# Dotlin plots: the share of cells expressing a gene as a bar below zero, and
# the distribution of the non-zero values as a violin above it. Private
# helpers (validation, data extraction, categories, palette, plot data and
# plot construction) come first, followed by the public `dotlin_plot()`.


# General utilities ------------------------------------------------------------

# Generate y-axis breaks for the expression values only. They are computed
# from 0 to the top of the axis, so that the space taken by the bars below
# zero does not make them coarser than on an ordinary plot of the same data.
nonnegative_breaks <- function(limits) {
  breaks <- scales::breaks_extended(n = 5)(c(0, limits[2]))
  breaks[breaks >= 0]
}

# Format values for error messages: 'a', 'b', 'c'.
format_values <- function(x) {
  paste0("'", x, "'", collapse = ", ")
}

check_installed <- function(pkg, version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this type of object.",
         call. = FALSE)
  }
  if (!is.null(version) && utils::packageVersion(pkg) < version) {
    stop("Package '", pkg, "' >= ", version, " is required, but version ",
         utils::packageVersion(pkg), " is installed.", call. = FALSE)
  }
}


# Input validation -------------------------------------------------------------

is_string <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x)
}

# Validate the arguments that do not depend on the type of `object`. Returns
# the genes to plot, without duplicates.
validate_inputs <- function(genes, category_col, palette, min_nonzero,
                            shared_y, layer, assay) {
  is_count <- is.numeric(min_nonzero) && length(min_nonzero) == 1 &&
    is.finite(min_nonzero) && min_nonzero >= 1 && min_nonzero %% 1 == 0

  if (!is_count) {
    stop("`min_nonzero` must be a positive integer.", call. = FALSE)
  }

  if (!isTRUE(shared_y) && !isFALSE(shared_y)) {
    stop("`shared_y` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(category_col) && !is_string(category_col)) {
    stop("`category_col` must be a single column name.", call. = FALSE)
  }

  if (!is.null(layer) && !is_string(layer)) {
    stop("`layer` must be a single layer name.", call. = FALSE)
  }

  if (!is.null(assay) && !is_string(assay)) {
    stop("`assay` must be a single assay name.", call. = FALSE)
  }

  genes <- unique(as.character(genes))

  if (length(genes) == 0) {
    stop("`genes` must contain at least one gene.", call. = FALSE)
  }

  if (!is.null(palette) && !is.character(palette) && !is.list(palette)) {
    stop("`palette` must be a character vector (or list) of colours, or NULL.",
         call. = FALSE)
  }

  genes
}


# Expression extraction --------------------------------------------------------

# Extract expression and categories from a supported object. Returns a list
# with `expression` (a numeric matrix with one row per cell and one column per
# gene), `category` (the category of each cell) and `category_label` (used for
# the axis and legend titles).
get_cell_data <- function(object, genes, category_col, layer, assay) {
  if (inherits(object, "Seurat")) {
    data <- get_seurat_data(object, genes, category_col, layer, assay)
  } else if (inherits(object, "SummarizedExperiment")) {
    data <- get_sce_data(object, genes, category_col, layer, assay)
  } else if (is.data.frame(object)) {
    data <- get_data_frame_data(object, genes, category_col, layer, assay)
  } else {
    stop("`object` must be a Seurat object, a SingleCellExperiment (or other ",
         "SummarizedExperiment) or a data frame.", call. = FALSE)
  }

  if (anyNA(data$expression)) {
    warning("Missing (NA) expression values were treated as zero.",
            call. = FALSE)
    data$expression[is.na(data$expression)] <- 0
  }

  if (any(data$expression < 0)) {
    warning("The expression data contain negative values (are they scaled?). ",
            "Only values > 0 count as expressed, so the proportions shown may ",
            "be misleading. Use non-negative data such as log-normalised ",
            "expression.", call. = FALSE)
  }

  data
}

get_seurat_data <- function(object, genes, category_col, layer, assay) {
  check_installed("SeuratObject", "5.0.0")

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  if (!assay %in% SeuratObject::Assays(object)) {
    stop("Assay '", assay, "' was not found. Available assays: ",
         format_values(SeuratObject::Assays(object)), ".", call. = FALSE)
  }

  assay_data <- object[[assay]]

  # Match the layer exactly, as Seurat's own lookup falls back to partial
  # matching. Layers split by sample (e.g. "data.1", "data.2") are combined.
  if (is.null(layer)) {
    layer <- "data"
  }

  available_layers <- SeuratObject::Layers(assay_data)

  if (layer %in% available_layers) {
    layers <- layer
  } else {
    layers <- available_layers[
      startsWith(available_layers, paste0(layer, "."))
    ]
  }

  if (length(layers) == 0) {
    stop("Layer '", layer, "' was not found in assay '", assay, "'. ",
         "Available layers: ", format_values(available_layers), ".",
         call. = FALSE)
  }

  available_genes <- unique(unlist(lapply(
    layers,
    function(l) SeuratObject::Features(assay_data, layer = l)
  )))

  missing_genes <- setdiff(genes, available_genes)

  if (length(missing_genes) > 0) {
    stop("The following genes were not found in layer '", layer, "' of ",
         "assay '", assay, "': ", format_values(missing_genes), call. = FALSE)
  }

  # Fetching from the assay rather than the Seurat object avoids clashes
  # between gene names and metadata columns.
  expression <- SeuratObject::FetchData(
    assay_data,
    vars = genes,
    layer = layers,
    clean = FALSE
  )
  expression <- as.matrix(expression[genes])

  if (is.null(category_col)) {
    category <- SeuratObject::Idents(object)
    category_label <- "Identity"
  } else {
    if (!category_col %in% colnames(object[[]])) {
      stop("'", category_col, "' was not found in the cell metadata of ",
           "`object`.", call. = FALSE)
    }

    category <- object[[category_col, drop = TRUE]]
    category_label <- category_col
  }

  list(
    expression = expression,
    category = category[match(rownames(expression), colnames(object))],
    category_label = category_label
  )
}

get_sce_data <- function(object, genes, category_col, layer, assay) {
  check_installed("SummarizedExperiment")

  if (!is.null(assay)) {
    stop("`assay` only applies to Seurat objects. For a ",
         "SingleCellExperiment, choose the assay with `layer`, e.g. ",
         "`layer = \"counts\"`.", call. = FALSE)
  }

  available_assays <- SummarizedExperiment::assayNames(object)

  if (length(available_assays) == 0) {
    stop("`object` has no named assays.", call. = FALSE)
  }

  if (is.null(layer)) {
    layer <- if ("logcounts" %in% available_assays) {
      "logcounts"
    } else {
      available_assays[1]
    }
  }

  if (!layer %in% available_assays) {
    stop("Assay '", layer, "' was not found. Available assays: ",
         format_values(available_assays), ".", call. = FALSE)
  }

  missing_genes <- setdiff(genes, rownames(object))

  if (length(missing_genes) > 0) {
    stop("The following genes were not found in rownames(object): ",
         format_values(missing_genes), call. = FALSE)
  }

  expression <- SummarizedExperiment::assay(object, layer)[genes, ,
                                                           drop = FALSE]
  expression <- t(as.matrix(expression))

  col_data <- SummarizedExperiment::colData(object)

  if (is.null(category_col)) {
    # The column used by SingleCellExperiment::colLabels().
    if (!"label" %in% colnames(col_data)) {
      stop("`category_col` must be supplied, as `object` has no cell labels ",
           "(colLabels).", call. = FALSE)
    }

    category_col <- "label"
  }

  if (!category_col %in% colnames(col_data)) {
    stop("'", category_col, "' was not found in colData(object).",
         call. = FALSE)
  }

  list(
    expression = expression,
    category = col_data[[category_col]],
    category_label = category_col
  )
}

get_data_frame_data <- function(object, genes, category_col, layer, assay) {
  if (!is.null(layer) || !is.null(assay)) {
    stop("`layer` and `assay` do not apply when `object` is a data frame.",
         call. = FALSE)
  }

  if (is.null(category_col)) {
    stop("`category_col` must be supplied when `object` is a data frame.",
         call. = FALSE)
  }

  if (!category_col %in% names(object)) {
    stop("'", category_col, "' is not a column of `object`.", call. = FALSE)
  }

  missing_genes <- setdiff(genes, names(object))

  if (length(missing_genes) > 0) {
    stop("The following genes were not found in the columns of `object`: ",
         format_values(missing_genes), call. = FALSE)
  }

  is_numeric <- vapply(genes, function(g) is.numeric(object[[g]]), logical(1))

  if (!all(is_numeric)) {
    stop("The following gene columns are not numeric: ",
         format_values(genes[!is_numeric]), call. = FALSE)
  }

  columns <- lapply(genes, function(g) as.double(object[[g]]))

  expression <- matrix(
    unlist(columns, use.names = FALSE),
    ncol = length(genes),
    dimnames = list(NULL, genes)
  )

  list(
    expression = expression,
    category = object[[category_col]],
    category_label = category_col
  )
}


# Category handling ------------------------------------------------------------

# Determine the categories and their plotting order: `category_order` if
# given, otherwise the factor levels (or sorted values) that occur in the data.
get_categories <- function(category, category_label, category_order = NULL) {
  observed <- levels(factor(category))

  # Numbers stored as text, such as cluster numbers, are sorted as numbers
  # ("2" before "10"), not as text.
  if (!is.factor(category)) {
    numbers <- suppressWarnings(as.numeric(observed))

    if (!anyNA(numbers)) {
      observed <- observed[order(numbers)]
    }
  }

  if (length(observed) == 0) {
    stop("'", category_label, "' has no non-missing categories.",
         call. = FALSE)
  }

  if (is.null(category_order)) {
    return(observed)
  }

  category_order <- unique(as.character(category_order))

  if (length(category_order) == 0) {
    stop("`category_order` cannot be empty.", call. = FALSE)
  }

  missing_categories <- setdiff(category_order, observed)

  if (length(missing_categories) > 0) {
    stop("The following categories in `category_order` were not found in '",
         category_label, "': ", format_values(missing_categories),
         call. = FALSE)
  }

  category_order
}


# Palette handling -------------------------------------------------------------

# Convert the palette specification into a named vector: category -> colour.
get_palette <- function(palette, categories) {
  if (is.null(palette)) {
    return(stats::setNames(scales::hue_pal()(length(categories)), categories))
  }

  palette <- unlist(palette)

  if (is.null(names(palette))) {
    if (length(palette) < length(categories)) {
      stop("`palette` has ", length(palette), " colours, but ",
           length(categories), " categories are shown.", call. = FALSE)
    }

    return(stats::setNames(palette[seq_along(categories)], categories))
  }

  missing_categories <- setdiff(categories, names(palette))

  if (length(missing_categories) > 0) {
    stop("The supplied palette is missing colours for the following ",
         "categories: ", format_values(missing_categories), call. = FALSE)
  }

  stats::setNames(palette[match(categories, names(palette))], categories)
}


# Plot data preparation --------------------------------------------------------

# Prepare long-format data and non-zero proportion information.
prepare_plot_data <- function(expression, category, categories, min_nonzero,
                              shared_y) {
  genes <- colnames(expression)

  keep <- as.character(category) %in% categories
  expression <- expression[keep, , drop = FALSE]
  category <- factor(as.character(category[keep]), levels = categories)

  # One row per cell and gene; genes keep the order supplied by the user.
  expression_long <- data.frame(
    category = rep(category, times = length(genes)),
    gene = factor(rep(genes, each = nrow(expression)), levels = genes),
    expression = as.vector(expression)
  )

  expression_long$x <- as.numeric(expression_long$category)

  # Categories in rows, genes in columns. `match()` rather than `[` so that
  # a category named "" is found too.
  nonzero_count <- rowsum((expression > 0) * 1, category)
  nonzero_count <- nonzero_count[match(categories, rownames(nonzero_count)), ,
                                 drop = FALSE]
  nonzero_prop <- nonzero_count / as.vector(table(category))

  cell_group <- cbind(
    as.integer(expression_long$category),
    as.integer(expression_long$gene)
  )
  expression_long$nonzero_count <- nonzero_count[cell_group]

  detection_data <- data.frame(
    category = factor(rep(categories, times = length(genes)),
                      levels = categories),
    gene = factor(rep(genes, each = length(categories)), levels = genes),
    nonzero_count = as.vector(nonzero_count),
    nonzero_prop = as.vector(nonzero_prop)
  )

  detection_data$x <- as.numeric(detection_data$category)

  # The bars are scaled to the highest value of all genes, so that every panel
  # has the same y-axis and the same bars, or to the highest value of each
  # gene when each gene has its own y-axis. Where nothing is expressed, 1 keeps
  # the bars (all at 0%) visible.
  if (shared_y) {
    max_expression <- rep(max(expression), length(genes))
  } else {
    max_expression <- apply(expression, 2, max)
  }

  max_expression[max_expression <= 0] <- 1

  detection_data$max_expression <- rep(max_expression,
                                       each = length(categories))

  bar_width <- 0.65

  detection_data$bar_height <- detection_data$max_expression * 0.15
  detection_data$bar_ymin <- -detection_data$bar_height
  detection_data$bar_ymax <- -detection_data$max_expression * 0.05
  detection_data$bar_xmin <- detection_data$x - bar_width / 2
  detection_data$bar_xmax <- detection_data$x + bar_width / 2
  detection_data$detection_xmax <- detection_data$bar_xmin +
    bar_width * detection_data$nonzero_prop
  # Rounding must not hide a few expressing cells ("0%") or a few cells
  # without expression ("100%").
  percent <- round(detection_data$nonzero_prop * 100)
  detection_data$label <- paste0(percent, "%")
  detection_data$label[percent == 0 & detection_data$nonzero_prop > 0] <- "<1%"
  detection_data$label[percent == 100 & detection_data$nonzero_prop < 1] <-
    ">99%"
  # Centre of the percentage below each bar in a single-gene plot, as in the
  # original plot. With several genes the percentage hangs from its bar and
  # room is kept down to `label_bottom`.
  detection_data$label_y <- -detection_data$bar_height * 1.5
  detection_data$label_bottom <- -detection_data$max_expression * 0.3

  strip_data <- expression_long[expression_long$expression > 0, ,
                                drop = FALSE]

  # A violin needs at least two points; ggplot2 would drop it with a warning.
  violin_data <- strip_data[strip_data$nonzero_count >= max(min_nonzero, 2), ,
                            drop = FALSE]

  list(
    expression_long = expression_long,
    strip_data = strip_data,
    violin_data = violin_data,
    detection_data = detection_data
  )
}


# Percentages ------------------------------------------------------------------

# The percentages are drawn at 7 pt when they fit. How much room they have is
# only known when the plot is drawn (in the RStudio pane, after resizing it,
# or by ggsave()), so they are sized then: if two neighbouring percentages
# would come closer than half the height of the text, or a percentage would
# reach past its bar or the edge of its panel, all percentages are made
# smaller by the same factor, just enough to fit.
# Named like ggplot2's own geoms.
GeomPercentage <- ggplot2::ggproto( # nolint: object_name_linter.
  "GeomPercentage", ggplot2::Geom,
  required_aes = c("x", "y", "ymin", "label"),
  default_aes = ggplot2::aes(colour = "#333333", size = 7 / ggplot2::.pt,
                             vjust = 0.5),
  draw_key = ggplot2::draw_key_blank,
  # `ymin` is the bottom of the bar above each percentage. `rows` holds the
  # percentages of every panel in x-axis order, so that all panels get the
  # same size.
  draw_panel = function(data, panel_params, coord, rows) {
    categories <- coord$transform(data.frame(x = c(1, 2), y = c(0, 0)),
                                  panel_params)
    grid::gTree(
      data = coord$transform(data, panel_params),
      rows = rows,
      category_spacing = abs(diff(categories$x)),
      cl = "dotlin_percentages"
    )
  }
)

#' @exportS3Method grid::makeContent
makeContent.dotlin_percentages <- function(x) {
  data <- x$data
  labels <- unique(unlist(x$rows))
  in_points <- function(npc, convert) {
    convert(grid::unit(npc, "npc"), "pt", valueOnly = TRUE)
  }

  # Room between two categories, and from each percentage down to the panel
  # edge and up to its bar.
  spacing <- in_points(x$category_spacing, grid::convertWidth)
  y <- in_points(data$y, grid::convertHeight)
  bar_bottom <- in_points(data$ymin, grid::convertHeight)

  # The factor by which text of the given size must shrink to fit (1 if it
  # fits): neighbouring percentages need a gap of half the text height, and
  # half a point is kept to the bar and to the panel edge.
  shrink_needed <- function(font_size) {
    gp <- grid::gpar(fontsize = font_size)
    widths <- vapply(labels, function(label) {
      grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = gp)),
                         "pt", valueOnly = TRUE)
    }, numeric(1))
    height <- grid::convertHeight(
      grid::grobHeight(grid::textGrob(labels, gp = gp)), "pt",
      valueOnly = TRUE
    )
    # Space needed between the centres of two neighbouring categories.
    needed <- vapply(x$rows, function(row) {
      w <- widths[row]
      if (length(w) < 2) 0 else max(w[-1] + w[-length(w)]) / 2
    }, numeric(1))

    min(
      1,
      spacing / (max(needed) + height / 2),
      (y - 0.5) / (data$vjust * height),
      ((bar_bottom - y - 0.5) / ((1 - data$vjust) * height))[data$vjust < 1]
    )
  }

  # Some devices, such as pdf(), draw text in whole points, so the text is
  # measured again at the smaller size until it fits.
  font_size <- data$size[1] * ggplot2::.pt
  for (i in 1:5) {
    shrink <- shrink_needed(font_size)
    if (shrink >= 1) {
      break
    }
    font_size <- max(font_size * shrink, 0.1)
  }

  grid::setChildren(x, grid::gList(grid::textGrob(
    data$label,
    x = grid::unit(data$x, "npc"),
    y = grid::unit(data$y, "npc"),
    vjust = data$vjust,
    gp = grid::gpar(col = data$colour, fontsize = font_size),
    name = "percentages"
  )))
}

# The percentages below the bars, centred on (`vjust = 0.5`) or hanging from
# (`vjust > 1`) the y position in column `y` of the detection data.
percentage_layer <- function(detection_data, y, vjust) {
  ggplot2::layer(
    geom = GeomPercentage,
    stat = "identity",
    position = "identity",
    data = detection_data,
    mapping = ggplot2::aes(
      x = .data$x,
      y = .data[[y]],
      ymin = .data$bar_ymin,
      label = .data$label
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    params = list(
      vjust = vjust,
      rows = unname(split(detection_data$label[order(detection_data$x)],
                          detection_data$gene[order(detection_data$x)]))
    )
  )
}


# Plot construction ------------------------------------------------------------

# Construct the ggplot.
build_plot <- function(plot_data, categories, category_label, palette, title,
                       shared_y, point_size, point_alpha, jitter_width) {
  single_gene <- nlevels(plot_data$detection_data$gene) == 1

  # A single gene gets a boxed title strip above its panel. Several genes are
  # stacked in rows named on the left, as a strip above each panel would take
  # height from every gene.
  if (single_gene) {
    gene_layout <- list(
      ggplot2::facet_wrap(~gene),
      # The box has the axis line's width and is not clipped, so both lines
      # are centred on the same edge and the y-axis continues into the box
      # without a step.
      ggplot2::theme(
        strip.background = ggplot2::element_rect(colour = "black",
                                                 linewidth = ggplot2::rel(1)),
        strip.clip = "off"
      )
    )
  } else {
    gene_layout <- list(
      ggplot2::facet_grid(
        gene ~ .,
        scales = if (shared_y) "fixed" else "free_y",
        switch = "y"
      ),
      ggplot2::theme(
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1)
      )
    )
  }

  # With one gene, each percentage is centred 1.5 bar heights below zero and
  # the panel ends just below it, as in the original plot. With several
  # genes, panels can be short, so the percentage hangs from its bar
  # (vjust > 1) to never overlap it, and a blank layer keeps room for it.
  if (single_gene) {
    label_layers <- list(
      percentage_layer(plot_data$detection_data, y = "label_y", vjust = 0.5)
    )
  } else {
    label_layers <- list(
      percentage_layer(plot_data$detection_data, y = "bar_ymin", vjust = 1.2),
      ggplot2::geom_blank(
        data = plot_data$detection_data,
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$label_bottom
        ),
        inherit.aes = FALSE
      )
    )
  }

  # Category names are angled only when some are long enough to run into
  # each other.
  if (max(nchar(categories)) > 5) {
    category_text <- ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
  } else {
    category_text <- ggplot2::element_text()
  }

  ggplot2::ggplot(
    plot_data$expression_long,
    ggplot2::aes(x = .data$x, y = .data$expression, fill = .data$category)
  ) +
    # Every violin has the same width, so that its shape shows only how the
    # expressing cells are distributed; how many cells express the gene is
    # shown by the bar.
    ggplot2::geom_violin(
      data = plot_data$violin_data,
      trim = TRUE,
      scale = "width",
      width = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_jitter(
      data = plot_data$strip_data,
      shape = 21,
      size = point_size,
      alpha = point_alpha,
      width = jitter_width,
      height = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_rect(
      data = plot_data$detection_data,
      mapping = ggplot2::aes(
        xmin = .data$bar_xmin,
        xmax = .data$bar_xmax,
        ymin = .data$bar_ymin,
        ymax = .data$bar_ymax
      ),
      fill = "#bdbdbd",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = plot_data$detection_data,
      mapping = ggplot2::aes(
        xmin = .data$bar_xmin,
        xmax = .data$detection_xmax,
        ymin = .data$bar_ymin,
        ymax = .data$bar_ymax,
        fill = .data$category
      ),
      inherit.aes = FALSE
    ) +
    label_layers +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "#444444",
      linewidth = 0.3
    ) +
    ggplot2::labs(
      title = title,
      fill = category_label,
      x = category_label,
      y = "Expression"
    ) +
    # Without explicit limits the legend would list categories in the order
    # the layers first meet them, putting those without a violin last.
    ggplot2::scale_fill_manual(
      values = palette,
      limits = categories
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(categories),
      labels = categories,
      expand = c(0.1, 0.1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = nonnegative_breaks
    ) +
    # In a very short panel a percentage may reach past the panel: better
    # than cutting it off.
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = category_text
    ) +
    gene_layout
}


# Public API -------------------------------------------------------------------

#' Dotlin plot of zero-inflated expression
#'
#' Shows the expression of each gene across categories of cells (clusters,
#' cell types, conditions, ...) in two parts that are read separately:
#'
#' * **Below zero**, a bar for each category stands for all of its cells
#'   (100%): the coloured part is the share with non-zero expression, the grey
#'   part the share without, and the label gives the percentage.
#' * **Above zero**, a violin and jittered points show the expression levels
#'   of the expressing cells only; zeros are left out. All violins have the
#'   same width, so their shapes can be compared even where few cells express
#'   the gene.
#'
#' Categories with fewer than `min_nonzero` expressing cells get points but no
#' violin, as a density estimated from a handful of values is not meaningful.
#' By default all gene panels share one y-axis, so genes can be compared
#' directly.
#'
#' The percentages are drawn at 7 pt. When the figure is too small for that
#' (many genes, many categories or a narrow figure), they are drawn just small
#' enough not to overlap each other, their bars or the panel edges; a larger
#' figure keeps them at full size.
#'
#' Only values above zero count as expressed, so use non-negative data such as
#' log-normalised expression (the default for Seurat and SingleCellExperiment
#' objects), not scaled data. A warning is given if negative values are found.
#'
#' @param object A Seurat object, a SingleCellExperiment (or other
#'   SummarizedExperiment) or a data frame with one row per cell, one numeric
#'   column per gene and a column of categories.
#' @param genes Names of the genes to show, one panel each, in the order given.
#' @param category_col Name of the cell metadata column that holds the
#'   categories. Defaults to the cell identities (`Idents()`) of a Seurat
#'   object or the cell labels (`colLabels()`) of a SingleCellExperiment.
#'   Required for data frames.
#' @param palette Colours of the categories: a named vector (or list) mapping
#'   each category to a colour, or an unnamed vector used in category order.
#'   If `NULL`, the default ggplot2 hue palette is used.
#' @param category_order Categories to show, in the order they should appear
#'   on the x-axis; categories that are not listed are left out. By default all
#'   categories are shown, in factor-level order (or sorted, if the column is
#'   not a factor; numbers stored as text are sorted as numbers).
#' @param min_nonzero Smallest number of expressing cells for which a violin
#'   is drawn. This is a number of cells, not a percentage: in a small
#'   category, a high percentage can still be too few cells for a violin.
#' @param shared_y If `TRUE` (the default), all genes share one y-axis. Set it
#'   to `FALSE` to give each gene its own y-axis, e.g. for raw counts, where
#'   genes can differ a lot in range.
#' @param layer Expression matrix to use. For a Seurat object, a layer of the
#'   assay (default `"data"`, the log-normalised expression). For a
#'   SingleCellExperiment, an assay (default `"logcounts"` if present,
#'   otherwise the first assay). Not used for data frames.
#' @param assay Seurat objects only: the assay to use. Defaults to
#'   `DefaultAssay(object)`.
#' @param title Plot title.
#' @param point_size Size of the points.
#' @param point_alpha Opacity of the points, between 0 and 1.
#' @param jitter_width Amount of horizontal jitter of the points, in each
#'   direction.
#'
#' @return A ggplot object, which can be customised further with the usual
#'   ggplot2 functions.
#'
#' @examples
#' # Simulated zero-inflated expression of one gene in three cell types
#' set.seed(1)
#' cells <- data.frame(
#'   cell_type = rep(c("A", "B", "C"), each = 200),
#'   gene_x = rbinom(600, 1, rep(c(0.05, 0.6, 0.3), each = 200)) * rlnorm(600)
#' )
#'
#' dotlin_plot(cells, "gene_x", category_col = "cell_type")
#'
#' \dontrun{
#' # Seurat object: cells are grouped by their identities by default
#' dotlin_plot(pbmc, genes = c("CST3", "NKG7", "PPBP"))
#'
#' # SingleCellExperiment: uses the "logcounts" assay by default
#' dotlin_plot(sce, genes = c("CST3", "NKG7", "PPBP"),
#'             category_col = "cell_type")
#' }
#'
#' @export
dotlin_plot <- function(object,
                        genes,
                        category_col = NULL,
                        palette = NULL,
                        category_order = NULL,
                        min_nonzero = 10,
                        shared_y = TRUE,
                        layer = NULL,
                        assay = NULL,
                        title = NULL,
                        point_size = 0.7,
                        point_alpha = 0.6,
                        jitter_width = 0.1) {
  genes <- validate_inputs(
    genes = genes,
    category_col = category_col,
    palette = palette,
    min_nonzero = min_nonzero,
    shared_y = shared_y,
    layer = layer,
    assay = assay
  )

  cell_data <- get_cell_data(
    object = object,
    genes = genes,
    category_col = category_col,
    layer = layer,
    assay = assay
  )

  categories <- get_categories(
    category = cell_data$category,
    category_label = cell_data$category_label,
    category_order = category_order
  )

  resolved_palette <- get_palette(
    palette = palette,
    categories = categories
  )

  plot_data <- prepare_plot_data(
    expression = cell_data$expression,
    category = cell_data$category,
    categories = categories,
    min_nonzero = min_nonzero,
    shared_y = shared_y
  )

  build_plot(
    plot_data = plot_data,
    categories = categories,
    category_label = cell_data$category_label,
    palette = resolved_palette,
    title = title,
    shared_y = shared_y,
    point_size = point_size,
    point_alpha = point_alpha,
    jitter_width = jitter_width
  )
}
