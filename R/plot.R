#' Spatial scatter plot of inferred subclones
#'
#' @description
#' Draws a spatial scatter plot from \code{cnv_obj@cluster}, coloring points by
#' clone label.
#' Optionally flips the axes with \code{coord_flip()}.
#'
#' @param cnv_obj highSpaClone object.
#' @param colors Character vector of hex/color names used for clusters.
#' @param point_size Numeric point size passed to \code{geom_point()}.
#' @param use_coord_flip Logical; if \code{TRUE}, flips the coordinate axes (useful when the
#'   image orientation requires rotation). Default is \code{FALSE}.
#' @param use_x_reverse Logical; if \code{TRUE}, reverses the x-axis direction
#'   (i.e., flips the plot horizontally). This option is useful when the
#'   spatial coordinates are mirrored relative to the original histology image
#'   or when you want to align the orientation with external annotations.
#'   Default is \code{FALSE}.
#' @param title Character plot title.
#'
#' @examples
#' \dontrun{
#' p <- spatialplot(cnv_obj,
#'                  point_size = 0.01,
#'                  use_x_reverse = T,
#'                  use_coord_flip = T,
#'                  title = "Tumor subclones")
#' print(p)
#' }
#'
#' @export
spatialplot <- function(cnv_obj,
                        colors = c("#ebe5c2","#D57358","#8a508f","#023047","#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85"),
                        point_size = 0.01,
                        use_x_reverse = FALSE,
                        use_coord_flip = FALSE,
                        title = "") {
  # --- input ---
  df <- cnv_obj@cluster
  req <- c("x", "y", "cell.label")

  # factorize cluster and freeze level order for stable color mapping
  df$cluster <- as.factor(df$cell.label)
  clv <- levels(df$cluster)
  n_clust <- length(clv)

  # --- plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$cluster)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(name = "Cluster", values = colors[seq_len(n_clust)]) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.text   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_text()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))

  if (isTRUE(use_coord_flip)) {
    p <- p + ggplot2::coord_flip()
  }

  if (isTRUE(use_x_reverse)) {
    p <- p + scale_x_reverse()
  }

  return(p)
}

#' Chromosome-binned CNV heatmap with group annotations
#'
#' @description
#' Draw a CNV heatmap split by chromosomes and annotated by a sample/cell grouping
#' variable. Designed for objects that provide `cnv.data` (cells × bins),
#' `annotation` (with `cell.id` and a grouping column), and `chr_pos` (per-chromosome
#' bin start indices). Supports on-screen display and saving to PDF/PNG with
#' optional auto-sizing based on data dimensions.
#'
#' @param cnv.obj A list-like object with components:
#'   \itemize{
#'     \item \code{cnv.data}: numeric matrix, rows = cells/spots, cols = bins.
#'     \item \code{annotation}: \code{data.frame} with at least columns
#'       \code{cell.id} and the column referenced by \code{groupby}.
#'     \item \code{chr_pos}: named numeric/integer vector or list mapping
#'       chromosome names to **0-based** start indices of bins
#'       (i.e., first bin index for chromosome \emph{k} is \code{start+1}).
#'   }
#' @param groupby Character scalar; column name in \code{cnv.obj$annotation}
#'   used to color rows (default \code{"cell.label"}).
#' @param color_panel Character vector of colors used for group annotation
#'   (recycled to the number of groups).
#' @param min_value,max_value Numeric; lower/upper limits for the heatmap
#'   color mapping (default \code{-3} and \code{3}).
#' @param show Logical; draw to current device (default \code{TRUE}).
#' @param save Logical; if \code{TRUE}, write a file to \code{outfile}.
#' @param outfile Character; output filename. If it ends with \code{.pdf} or
#'   \code{.png}, that format is used. Otherwise both PDF and PNG are written
#'   with sensible names.
#' @param auto_size Logical; if \code{TRUE} (default) choose device size based on
#'   matrix dimensions.
#' @param width_in,height_in Numeric; device width/height in inches used when
#'   \code{auto_size = FALSE}. Ignored otherwise.
#' @param res Integer; PNG resolution (dpi) when saving PNG (default \code{300}).
#' @param ... Passed to \code{ComplexHeatmap::Heatmap()}.
#'
#' @return (Invisibly) a \code{HeatmapList} object.
#'
#' @examples
#' \dontrun{
#' ht <- cnv_heatmap(
#'   cnv.obj,
#'   groupby = "cell.label",
#'   show = TRUE,
#'   save = TRUE,
#'   outfile = "cnv_heatmap.pdf"
#' )
#' }
#'
#'
#' @importFrom grid unit gpar
#' @importFrom grDevices pdf png dev.off
#' @export
#'
cnv_heatmap <- function(
    cnv.obj,
    groupby     = "cell.label",
    color_panel = c("#ffadad","#ffd6a5","#fdffb6","#caffbf",
                    "#9bf6ff","#a0c4ff","#bdb2ff","#ffc6ff"),
    min_value   = -3,
    max_value   =  3,
    show        = TRUE,
    save        = FALSE,
    outfile     = "cnv_heatmap.pdf",
    auto_size   = TRUE,
    width_in    = 10,
    height_in   = 6,
    res         = 300,
    ...
){
  # ---- checks & extraction ----
  if (is.null(cnv.obj$cnv.data) || is.null(cnv.obj$annotation) || is.null(cnv.obj$chr_pos)) {
    stop("`cnv.obj` must contain `cnv.data`, `annotation`, and `chr_pos`.")
  }
  cnv_mat <- as.matrix(cnv.obj$cnv.data)
  if (is.null(rownames(cnv_mat))) stop("`cnv.obj$cnv.data` must have rownames = cell IDs.")
  annotation_df <- cnv.obj$annotation
  if (!all(c("cell.id", groupby) %in% colnames(annotation_df))) {
    stop("`annotation` must contain columns `cell.id` and `", groupby, "`.")
  }

  # align rows to annotation cell.id
  row_idx <- match(annotation_df$cell.id, rownames(cnv_mat))
  valid   <- !is.na(row_idx)
  if (!any(valid)) stop("No annotation cell.id matched CNV rownames.")
  row_idx <- row_idx[valid]
  cnv_mat <- cnv_mat[row_idx, , drop = FALSE]
  group_vec <- annotation_df[[groupby]][valid]

  # order rows by group
  ord <- order(group_vec)
  cnv_mat   <- cnv_mat[ord, , drop = FALSE]
  sorted_grp <- group_vec[ord]

  # ---- chromosome splitting from chr_pos (0-based starts) ----
  chr_pos <- cnv.obj$chr_pos
  # accept list or numeric vector; coerce to named numeric vector
  if (is.list(chr_pos)) chr_pos <- unlist(chr_pos, use.names = TRUE)
  if (is.null(names(chr_pos))) stop("`chr_pos` must be named with chromosome names.")
  all_chr <- names(chr_pos)
  # offsets: starts then final end at ncol
  offsets <- c(unname(chr_pos[all_chr]), ncol(cnv_mat))
  if (any(is.na(offsets))) stop("`chr_pos` contains NA or missing names.")

  # map each bin (column) to a chromosome
  bin_to_chr <- character(ncol(cnv_mat))
  for (i in seq_along(all_chr)) {
    start0 <- offsets[i]
    end0   <- offsets[i + 1]
    cols   <- (start0 + 1):end0
    cols   <- cols[cols <= ncol(cnv_mat)]
    if (length(cols)) bin_to_chr[cols] <- all_chr[i]
  }
  chr_fac <- factor(bin_to_chr, levels = all_chr)

  # ---- colors & annotations ----
  col_fun <- circlize::colorRamp2(c(min_value, 0, max_value),
                                  c("blue", "white", "red"))

  group_fac   <- factor(sorted_grp)
  grp_levels  <- levels(group_fac)
  grp_colors  <- structure(rep(color_panel, length.out = length(grp_levels)),
                           names = grp_levels)

  ha_row <- ComplexHeatmap::rowAnnotation(
    cluster = group_fac,
    col     = list(cluster = grp_colors),
    show_annotation_name = TRUE,
    annotation_legend_param = list(cluster = list(title = groupby))
  )

  # ---- build heatmap ----
  ht <- ComplexHeatmap::Heatmap(
    cnv_mat,
    name               = "CNV",
    col                = col_fun,
    cluster_rows       = FALSE,
    cluster_columns    = FALSE,
    show_row_names     = FALSE,
    show_column_names  = FALSE,
    column_split       = chr_fac,
    row_split          = group_fac,
    row_title          = NULL,
    border             = TRUE,
    rect_gp            = grid::gpar(col = NA),
    row_gap            = grid::unit(0.1, "pt"),
    column_gap         = grid::unit(0.1, "pt"),
    column_title_gp    = grid::gpar(fontsize = 9, fontface = "bold"),
    column_title_rot   = 90,
    heatmap_legend_param = list(
      direction      = "horizontal",
      legend_width   = grid::unit(2, "cm"),
      legend_height  = grid::unit(0.4, "cm"),
      at             = c(min_value, max_value),
      labels         = c("Loss", "Gain"),
      title_position = "leftcenter"
    ),
    ...
  )

  ht_list <- ha_row + ht

  # helper to draw separators between groups
  .decorate_groups <- function() {
    rle_grp   <- rle(as.character(sorted_grp))
    boundaries <- cumsum(rle_grp$lengths)
    if (length(boundaries) > 1L) {
      for (b in boundaries[-length(boundaries)]) {
        grid::grid.lines(x = c(0, 1),
                         y = grid::unit(b, "native"),
                         gp = grid::gpar(lwd = 1, col = "black"))
      }
    }
  }

  # ---- auto size (inches) if requested ----
  if (auto_size) {
    nr <- nrow(cnv_mat); nc <- ncol(cnv_mat)
    width_in  <- max(6, min(20, 3 + nc * 0.01))   # 0.01 in per column + margins
    height_in <- max(4, min(20, 2 + nr * 0.006))  # 0.006 in per row + margins
  }

  # ---- save to file(s) ----
  if (isTRUE(save)) {
    ext <- tolower(tools::file_ext(outfile))
    write_pdf <- isTRUE(ext == "pdf")
    write_png <- isTRUE(ext == "png")

    if (!write_pdf && !write_png) {
      # write both if no recognized extension
      pdf_file <- sub("\\.png$", ".pdf", outfile)
      png_file <- sub("\\.pdf$", ".png", outfile)
      # PDF
      grDevices::pdf(pdf_file, width = width_in, height = height_in)
      ComplexHeatmap::draw(ht_list, merge_legend = TRUE,
                           annotation_legend_side = "bottom",
                           heatmap_legend_side = "bottom")
      ComplexHeatmap::decorate_heatmap_body("CNV", .decorate_groups())
      grDevices::dev.off()
      # PNG
      grDevices::png(png_file,
                     width = round(width_in * res),
                     height = round(height_in * res),
                     res = res)
      ComplexHeatmap::draw(ht_list, merge_legend = TRUE,
                           annotation_legend_side = "bottom",
                           heatmap_legend_side = "bottom")
      ComplexHeatmap::decorate_heatmap_body("CNV", .decorate_groups())
      grDevices::dev.off()
    } else {
      if (write_pdf) {
        grDevices::pdf(outfile, width = width_in, height = height_in)
        ComplexHeatmap::draw(ht_list, merge_legend = TRUE,
                             annotation_legend_side = "bottom",
                             heatmap_legend_side = "bottom")
        ComplexHeatmap::decorate_heatmap_body("CNV", .decorate_groups())
        grDevices::dev.off()
      }
      if (write_png) {
        grDevices::png(outfile,
                       width = round(width_in * res),
                       height = round(height_in * res),
                       res = res)
        ComplexHeatmap::draw(ht_list, merge_legend = TRUE,
                             annotation_legend_side = "bottom",
                             heatmap_legend_side = "bottom")
        ComplexHeatmap::decorate_heatmap_body("CNV", .decorate_groups())
        grDevices::dev.off()
      }
    }
  }

  # ---- on-screen draw ----
  if (isTRUE(show)) {
    ComplexHeatmap::draw(ht_list, merge_legend = TRUE,
                         annotation_legend_side = "bottom",
                         heatmap_legend_side = "bottom")
    ComplexHeatmap::decorate_heatmap_body("CNV", .decorate_groups())
  }

  invisible(ht_list)
}



