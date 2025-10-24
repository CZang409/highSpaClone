#################################################################################################################
# Package: highSpaClone
# Date : 2025-10-10
# Title : Copy Number Variation Inference and Tumor Subclone Analysis for High-Resolution Spatial Transcriptomics
# Authors: Chenxuan Zang
# Contacts: czang@mdanderson.org
#          The University of Texas MD Anderson, Department of Biostatistics
#################################################################################################################

#' The highSpaClone Class
#'
#' @slot counts.data Gene-by-cell raw count matrix.
#' @slot location The spatial coordinates of each cell.
#' @slot chr_pos The chromosome position corresponding to each bin, used for subsequent CNV heatmap drawing.
#' @slot smoothed.data Gene-by-cell smoothed count matrix by gene order.
#' @slot cnv.data Final CNV matrix. Dimension: cell-by-bin.
#' @slot cluster Tumor/Subclone label for each cell.
#' @slot annotation Cell type annotation information.
#' @slot gene_order Gene ordering information, contains chromosome, start, and stop position for each gene.
#' @slot project The name of the project.
#'
#'

setClass("highSpaClone", slots=list(
  counts.data = "ANY",
  location = "ANY",
  chr_pos = "ANY",
  smoothed.data = "ANY",
  cnv.data = "ANY",
  cluster = "ANY",
  annotation = "data.frame",
  gene_order= "data.frame",
  project = "character"
))

#' Create a highSpaClone object.
#'
#' @param counts Gene-by-cell raw count matrix.
#' @param location A data.frame or matrix of spatial coordinates for each cell, with three columns (cell.id, x, y). The row names of location must match with the column names of counts.
#' @param min_avg_expression Minimum average expression per gene.
#' @param min_gene_counts Minimum total counts per cell.
#' @param gene_order_file The gene order file, contains chromosome, start, and stop position for each gene.
#' @param annotations_file Cell type annotation information. Should contain two columns: cell.id and cell.label. Must be a data.frame, .csv or .txt file.
#' @param project The name of the project.
#'
#' @return Return a preliminary highSpaClone object.
#' @importFrom tools file_ext
#' @importFrom utils read.csv read.table
#' @importFrom dplyr filter mutate select arrange left_join
#' @importFrom magrittr `%>%`
#' @export
#'
createObject <- function(counts,
                         location,
                         min_avg_expression=0.01,
                         min_gene_counts=100,
                         gene_order_file=NULL,
                         annotations_file=NULL,
                         project=""){

  if (any(duplicated(rownames(counts)))) {
    stop("Ensure your gene names are unique!")
  }

  ## check dimension
  if(ncol(counts)!=nrow(location)){
    stop("The number of spots in counts and location should be consistent.")
  }

  ## check data order
  if(!identical(colnames(counts), rownames(location))){
    stop("The column names of counts and row names of location should be matched.")
  }

  colnames(location) <- c('cell.id', 'x', 'y')

  ## Filter genes by average expression
  if(!is.null(min_avg_expression)) {

    avg_expression <- Matrix::rowMeans(counts)
    genes_to_keep <- avg_expression >= min_avg_expression
    remained_gene_count <- sum(genes_to_keep)
    removed_gene_count <- sum(!genes_to_keep)
    counts <- counts[genes_to_keep, ]

    message(sprintf("Filtered genes based on average expression threshold: %d genes removed.", removed_gene_count))
    message(sprintf("Remaining genes after filtering: %d genes.", remained_gene_count))

  }
  ## Filter cells by minimum gene counts
  if(!is.null(min_gene_counts)) {
    cell_counts <- Matrix::colSums(counts)
    cells_to_keep <- cell_counts >= min_gene_counts
    remained_cells_count <- sum(cells_to_keep)
    removed_cells_count <- sum(!cells_to_keep)
    counts <- counts[, cells_to_keep]

    message(sprintf("Filtered cells based on minimum gene counts threshold: %d cells removed.", removed_cells_count))
    message(sprintf("Remaining cells after filtering: %d cells.", remained_cells_count))

  }

  #### get gene order info
  if(is.null(gene_order_file)){

    data(hg38_annotation)
    gene_order <- hg38_annotation

  }else if(is.matrix(gene_order_file)||is.data.frame(gene_order_file)){
    gene_order <- gene_order_file
  }else{
    extension <- tools::file_ext(gene_order_file)
    if (extension == "txt") {
      gene_order <- read.table(gene_order_file, header=FALSE, sep="\t")
    } else if (extension == "csv") {
      gene_order <- read.csv(gene_order_file, header=TRUE)
    } else {
      stop("Unsupported file format. Only txt and csv files are supported.")
    }
  }

  colnames(gene_order) <- c("gene","chromosome", "start", "end")
  rownames(gene_order) <- gene_order$gene

  ## extract common genes of count matrix and gene file
  common_genes <- intersect(rownames(counts), gene_order$gene)
  annotation_genes <- gene_order$gene
  left_genes <- intersect(common_genes, annotation_genes)
  removed_genes_count <- length(common_genes) - length(left_genes)
  remaining_genes_count <- length(left_genes)
  message(sprintf("Removed %d genes that did not match the gene order file.", removed_genes_count))
  message(sprintf("Remaining %d genes after matching with the gene order file.", remaining_genes_count))

  gene_order <- gene_order[gene_order$gene %in% common_genes,]
  counts <- counts[rownames(counts) %in% common_genes,]

  #### read annotations file
  if(is.null(annotations_file)){
    stop("Please input an annotation file. It can be a data.frame, .csv or .txt file.")

  } else if(is.data.frame(annotations_file)){
    annotations <- annotations_file

  } else {
    extension <- tools::file_ext(annotations_file)
    if (extension == "txt") {
      annotations <- read.table(annotations_file, header=FALSE, sep="\t",
                                stringsAsFactors=FALSE,
                                colClasses = c('character', 'character'))
    } else if (extension == "csv") {
      annotations <- read.table(annotations_file, header=TRUE, sep=",",
                                stringsAsFactors=FALSE,
                                colClasses = c('character', 'character'))
    } else {
      stop("Unsupported file format. Only dataframe, csv and txt files are supported.")
    }
  }

  colnames(annotations) <- c("cell.id", "cell.label")

  ## sort
  gene_order2 <- gene_order
  chr_num <- gsub("chr", "", gene_order2$chromosome)
  chr_num[chr_num=='X'] <- 23
  chr_num[chr_num=='Y'] <- 24
  chr_num[chr_num=='M'] <- 25
  gene_order2$chr_num <- as.numeric(chr_num)
  sorted_indices <- with(gene_order2, order(chr_num, start, end))
  gene_order <- gene_order[sorted_indices,,drop=FALSE]
  chr_levels <- unique(gene_order$chromosome)
  gene_order$chromosome <- factor(gene_order$chromosome,levels = chr_levels)
  counts <- counts[match(rownames(gene_order), rownames(counts)),,drop=FALSE]
  annotations <- annotations %>% dplyr::filter(cell.id %in% colnames(counts))
  annotations <- annotations[match(colnames(counts),annotations$cell.id),]

  object <- new(

    Class = "highSpaClone",
    counts.data = counts,
    location = location,
    annotation = annotations,
    gene_order= gene_order,
    project = ""

  )

  return(object)
}


#' Smooth expression along chromosomes (windowed running mean) with optional chunking/parallelism.
#'
#' @param obj highSpaClone object.
#' @param window_size Integer. Size of the running-mean window (in genes) used within each chromosome (default: 101).
#' @param step Integer. Number of genes in a bin (default: 50).
#' @param exclude_chromosomes Chromosomes that are excluded from analysis (default: c('chrX', 'chrY', 'chrM')).
#' @param smooth_with_ends Logical. Whether to apply an end-aware smoother (default: \code{FALSE}).
#' @param use_chunk Logical. If \code{TRUE}, split cells into chunks of size \code{chunk_size}
#' @param chunk_size Integer. Number of cells per chunk when \code{use_chunk=TRUE}.
#' @param parallel Logical. If \code{TRUE}, process chunks in parallel.
#' @param n_cores Integer or \code{NULL}. Number of CPU cores to use when \code{parallel=TRUE}.
#'   If \code{NULL}, uses \code{parallel::detectCores()-1}. The effective cores are capped
#'   by \code{max_cores} and the number of chunks (default \code{NULL}).
#' @param max_cores Integer. Upper bound on cores even if more are available (default \code{4}).
#'
#' @return Return highSpaClone object containing smoothed gene expression data
#' @importFrom Matrix rBind
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster mclapply
#' @export
#'
smooth_expr <- function(
    obj,
    window_size = 101,
    step = 50,
    exclude_chromosomes = c("chrX", "chrY", "chrM"),
    smooth_with_ends = FALSE,
    use_chunk = TRUE,
    chunk_size = 5000,
    parallel = FALSE,
    n_cores = NULL,
    max_cores = 4
){

  start_time <- Sys.time()

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  ## ====== Extract expression matrix and normalize ======
  expr <- obj@counts.data
  expr <- normalize_data(expr)

  step1_time <- Sys.time()
  cat("Normalization time: ", difftime(step1_time, start_time, units="secs"), " seconds\n")

  ## Preserve sparsity
  is_sparse <- inherits(expr, "sparseMatrix")
  if (is_sparse) {
    expr <- as(expr, "CsparseMatrix")
    cat("Working with sparse matrix format (keeping sparse)\n")
  }

  var <- obj@gene_order
  n_cells <- ncol(expr)
  cat("Total cells: ", n_cells, "\n")

  ## ====== Decide whether to chunk by cells ======
  if (use_chunk && n_cells > chunk_size) {
    cat("\n=== Chunking mode enabled ===\n")
    cat("Processing in chunks of size ", chunk_size, "...\n")

    n_chunks <- ceiling(n_cells / chunk_size)
    chunk_indices <- split(seq_len(n_cells), ceiling(seq_along(seq_len(n_cells)) / chunk_size))
    cat("Total chunks: ", n_chunks, "\n")

    ## ====== Parallel or sequential processing ======
    if (parallel) {
      cat("=== Parallel processing enabled ===\n")

      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package 'parallel' is required for parallel processing.")
      }

      is_unix <- .Platform$OS.type == "unix"

      if (is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1L
        n_cores <- max(1L, n_cores)
      }
      n_cores <- min(n_cores, max_cores, n_chunks)
      cat("Requested cores: ", n_cores, " (max_cores limit: ", max_cores, ")\n")

      expr_size_gb <- as.numeric(object.size(expr)) / (1024^3)
      estimated_memory_gb <- expr_size_gb * n_cores * 1.5
      cat("Estimated memory usage: ~", round(estimated_memory_gb, 2), " GB\n")
      if (estimated_memory_gb > 50) {
        warning("High memory usage estimated (", round(estimated_memory_gb, 2),
                " GB). Consider reducing n_cores or max_cores, or using parallel=FALSE.")
      }

      if (is_unix) {
        cat("Unix system detected: using mclapply\n")
        cat("Using ", n_cores, " cores\n")

        chunk_results <- parallel::mclapply(seq_along(chunk_indices), function(i) {
          tryCatch({
            if (i %% 5 == 0) cat("Processing chunk", i, "of", n_chunks, "\n")
            idx <- chunk_indices[[i]]
            expr_chunk <- expr[, idx, drop = FALSE]
            infercnv_chunk(expr_chunk, var, exclude_chromosomes,
                           window_size, step, smooth_with_ends)
          }, error = function(e) {
            structure(list(error = conditionMessage(e)), class = "try-error")
          })
        }, mc.cores = n_cores)

      } else {
        cat("Windows system detected: using parLapply\n")
        cat("Using ", n_cores, " cores\n")

        cl <- parallel::makeCluster(n_cores)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

        parallel::clusterExport(
          cl,
          c("infercnv_chunk", "running_mean_by_chromosome",
            "running_mean_for_chromosome", "running_mean",
            "smooth_with_ends_internal", "sort_chromosomes",
            "var", "exclude_chromosomes", "window_size",
            "step", "smooth_with_ends"),
          envir = environment()
        )
        parallel::clusterEvalQ(cl, { base::loadNamespace("Matrix"); NULL })
        parallel::clusterEvalQ(cl, { base::loadNamespace("stats"); NULL })

        chunk_results <- parallel::parLapply(cl, seq_along(chunk_indices), function(i) {
          tryCatch({
            if (i %% 5 == 0) cat("Processing chunk", i, "of", length(chunk_indices), "\n")
            idx <- chunk_indices[[i]]
            expr_chunk <- expr[, idx, drop = FALSE]
            infercnv_chunk(expr_chunk, var, exclude_chromosomes,
                           window_size, step, smooth_with_ends)
          }, error = function(e) {
            structure(list(error = conditionMessage(e)), class = "try-error")
          })
        })
      }

    } else {
      cat("=== Sequential chunk processing ===\n")
      cat("Processing ", n_chunks, " chunks sequentially...\n")

      chunk_results <- lapply(seq_along(chunk_indices), function(i) {
        tryCatch({
          if (i %% 5 == 0 || i == 1 || i == n_chunks) {
            cat("Processing chunk", i, "of", n_chunks,
                sprintf("(%.1f%% complete)\n", 100*i/n_chunks))
          }
          idx <- chunk_indices[[i]]
          expr_chunk <- expr[, idx, drop = FALSE]
          res <- infercnv_chunk(expr_chunk, var, exclude_chromosomes,
                                window_size, step, smooth_with_ends)
          rm(expr_chunk)
          if (i %% 10 == 0) gc(verbose = FALSE)
          res
        }, error = function(e) {
          structure(list(error = conditionMessage(e)), class = "try-error")
        })
      })

      cat("Sequential chunk processing completed.\n")
    }

    ## ====== Validate chunk results and extract ======
    bad <- vapply(chunk_results, function(x) inherits(x, "try-error") || is.null(x$error) == FALSE, logical(1))
    if (any(bad)) {
      msgs <- vapply(chunk_results[bad], function(x) x$error, character(1))
      stop("Some chunks failed: ", paste(sprintf("#%d: %s", which(bad), msgs), collapse = " | "))
    }

    x_list <- lapply(chunk_results, `[[`, "x_smoothed")
    chr_pos_list <- lapply(chunk_results, `[[`, "chr_pos")

    ## Column-name consistency assertion (all chunks must share the same set/order)
    coln0 <- colnames(x_list[[1]])
    same_cols <- vapply(x_list, function(m) identical(colnames(m), coln0), logical(1))
    if (!all(same_cols)) {
      stop("Column names (smoothed positions) differ across chunks. ",
           "This indicates inconsistent per-chromosome outputs. Check running_mean/running_mean_by_chromosome.")
    }

    ## Combine (stack cells by rows)
    cat("Combining results...\n")
    if (inherits(x_list[[1]], "sparseMatrix")) {
      cat("Combining sparse matrices...\n")
      res <- do.call(Matrix::rBind, x_list)
    } else {
      cat("Combining dense matrices...\n")
      res <- do.call(rbind, x_list)
    }

    cat("Result class: ", class(res)[1], "\n")
    cat("Result dimensions: ", dim(res)[1], " x ", dim(res)[2], "\n")

    ## Use chr_pos from the first chunk; optionally check consistency
    chr_pos <- chr_pos_list[[1]]
    same_chrpos <- vapply(chr_pos_list, function(cp) identical(cp, chr_pos), logical(1))
    if (!all(same_chrpos)) {
      warning("chr_pos differs across chunks; using the first chunk's chr_pos.")
    }

    ## Set row names to the actual concatenation order
    cell_order <- unlist(chunk_indices, use.names = FALSE)
    rownames(res) <- colnames(obj@counts.data)[cell_order]

  } else {
    ## ====== No chunking ======
    if (!use_chunk) {
      cat("\n=== No chunking mode (use_chunk=FALSE) ===\n")
    } else {
      cat("\n=== Cell count below chunk size ===\n")
    }
    cat("Processing all ", n_cells, " cells at once...\n")

    expr_size_gb <- as.numeric(object.size(expr)) / (1024^3)
    cat("Expression matrix size: ~", round(expr_size_gb, 2), " GB\n")
    if (expr_size_gb > 10) {
      warning("Large matrix detected (", round(expr_size_gb, 2),
              " GB). Consider using use_chunk=TRUE to reduce memory usage.")
    }

    cat("Running smoothing on all cells...\n")
    results <- infercnv_chunk(expr, var, exclude_chromosomes,
                              window_size, step, smooth_with_ends)
    chr_pos <- results$chr_pos
    res <- results$x_smoothed

    cat("Result class: ", class(res)[1], "\n")
    cat("Result dimensions: ", dim(res)[1], " x ", dim(res)[2], "\n")

    ## Row names follow the original order
    rownames(res) <- colnames(obj@counts.data)

    cat("No-chunk processing completed.\n")
  }

  ## ====== Store into object ======
  cat("\nStoring results in object...\n")
  obj@smoothed.data <- res
  obj@chr_pos <- chr_pos

  cat("Results successfully stored.\n")

  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units="mins")
  cat("\n========================================\n")
  cat("Total running time: ", round(total_time, 2), " minutes\n")
  cat("========================================\n")

  return(obj)
}

## ====== Helper functions ======

infercnv_chunk <- function(expr, var, exclude_chromosomes, window_size, step, smooth_with_ends) {
  running_means <- running_mean_by_chromosome(expr, var, exclude_chromosomes,
                                              window_size, step, smooth_with_ends)
  list(
    chr_pos    = running_means$chr_pos,
    x_smoothed = running_means$x_smoothed
  )
}

sort_chromosomes <- function(l) {
  extract_number <- function(text) {
    num <- as.integer(gsub("[^0-9]", "", text))
    return(num)
  }
  l[order(sapply(l, extract_number))]
}

normalize_data <- function(expr, size_factor = NULL) {
  if (is.null(size_factor)) {
    if (inherits(expr, "sparseMatrix")) {
      cs <- Matrix::colSums(expr, na.rm = TRUE)
    } else {
      expr[is.na(expr)] <- 0
      cs <- Matrix::colSums(expr, na.rm = TRUE)
    }
    cs[cs == 0] <- 1
    size_factor <- stats::median(cs) / cs
  }
  if (inherits(expr, "sparseMatrix")) {
    normalized_data <- expr %*% Matrix::Diagonal(x = as.numeric(size_factor))
  } else {
    normalized_data <- sweep(as.matrix(expr), MARGIN = 2, STATS = size_factor, FUN = "*")
  }
  normalized_data
}

running_mean_for_chromosome <- function(chr, expr, var, window_size, step, smooth_with_ends) {
  idx <- which(var$chromosome == chr)
  if (length(idx) == 0L) {
    return(list(smoothed_x = matrix(numeric(0), nrow = ncol(expr), ncol = 0)))
  }

  ord <- order(var$start[idx], na.last = TRUE)
  genes_chr <- rownames(var)[idx][ord]

  keep <- genes_chr[genes_chr %in% rownames(expr)]
  if (length(keep) == 0L) {
    return(list(smoothed_x = matrix(numeric(0), nrow = ncol(expr), ncol = 0)))
  }

  expr_chr <- expr[keep, , drop = FALSE]  # genes x cells

  running_mean_res <- running_mean(expr_chr, window_size, step, smooth_with_ends)
  smoothed_x <- running_mean_res$smoothed_x   # cells x positions

  if (!is.matrix(smoothed_x)) smoothed_x <- as.matrix(smoothed_x)
  ## Explicit: rows = cells; row names = cell names
  stopifnot(nrow(smoothed_x) == ncol(expr_chr))
  if (is.null(rownames(smoothed_x))) rownames(smoothed_x) <- colnames(expr_chr)

  list(smoothed_x = smoothed_x)
}

running_mean_by_chromosome <- function(expr, var, exclude_chromosomes, window_size, step, smooth_with_ends) {
  if (!is.null(exclude_chromosomes)) {
    var <- var[!(var$chromosome %in% exclude_chromosomes), , drop = FALSE]
  }

  chromosomes <- unique(var$chromosome)
  chromosomes <- sort_chromosomes(chromosomes)

  running_means <- lapply(chromosomes, function(chr) {
    running_mean_for_chromosome(chr, expr, var, window_size, step, smooth_with_ends)
  })
  names(running_means) <- chromosomes

  running_means_list <- lapply(running_means, `[[`, 1)
  non_empty <- vapply(running_means_list, function(x) ncol(x) > 0, logical(1))
  running_means_list <- running_means_list[non_empty]
  chromosomes <- chromosomes[non_empty]
  if (length(running_means_list) == 0) stop("No valid chromosome data after filtering")

  ## Each chromosome must yield the same number of rows (cells)
  rn_ok <- length(unique(vapply(running_means_list, nrow, integer(1)))) == 1
  if (!rn_ok) stop("Inconsistent number of rows (cells) across chromosomes.")

  x_smoothed <- do.call(cbind, running_means_list)

  ## Fill row names (cells) if missing
  if (is.null(rownames(x_smoothed)) || anyNA(rownames(x_smoothed))) {
    rn <- rownames(running_means_list[[1]])
    if (!is.null(rn)) rownames(x_smoothed) <- rn
  }

  ## Build chr_pos (offset of starting column for each chromosome)
  chr_start_pos <- list()
  chr_pos_offset <- 0L
  for (i in seq_along(chromosomes)) {
    chr_start_pos[[chromosomes[i]]] <- chr_pos_offset
    chr_pos_offset <- chr_pos_offset + ncol(running_means_list[[i]])
  }
  names(chr_start_pos) <- chromosomes

  ## Column names: chr_i_j
  name_vectors <- unlist(lapply(seq_along(running_means_list), function(i){
    n <- ncol(running_means_list[[i]])
    paste(chromosomes[i], seq_len(n), sep = "_")
  }), use.names = FALSE)
  colnames(x_smoothed) <- name_vectors

  list(chr_pos = chr_start_pos, x_smoothed = x_smoothed)
}

smooth_with_ends_internal <- function(x, pyramid) {
  n <- length(pyramid)
  tail_length <- floor(n / 2)
  obs_length <- length(x)

  smoothed <- stats::filter(x, pyramid/sum(pyramid), sides = 2)
  smoothed <- as.vector(smoothed)

  for (tail_end in seq_len(tail_length)) {
    end_tail <- obs_length - tail_end + 1
    d_left <- tail_end - 1
    d_right <- obs_length - tail_end
    d_right <- ifelse(d_right > tail_length, tail_length, d_right)
    r_left <- tail_length - d_left
    r_right <- tail_length - d_right
    denominator <- (n^2 + n) / 2 - ((r_left * (r_left + 1)) / 2) - ((r_right * (r_right + 1)) / 2)
    numerator_counts_vector <- c(seq_len(tail_length), tail_length + 1, rev(seq_len(tail_length)))

    left_chunk <- x[seq_len(tail_end + d_right)]
    numerator_range_left <- numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
    smoothed[tail_end] <- sum(left_chunk * numerator_range_left) / denominator

    right_chunk <- x[(end_tail - d_right):obs_length]
    numerator_range_right <- numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
    smoothed[end_tail] <- sum(right_chunk * rev(numerator_range_right)) / denominator
  }
  smoothed
}

running_mean <- function(x, n = 51, step = 10, smooth_with_ends = FALSE) {
  ## x: genes x cells
  x <- Matrix::t(x)  # -> cells x genes
  n_cells <- nrow(x)
  n_genes <- ncol(x)

  if (n < 2) {
    smoothed_x <- Matrix::t(x)  # Return structure consistent with other branches
    return(list(smoothed_x = smoothed_x))
  }

  if (n < n_genes) {
    pyramid <- pmin(seq_len(n), rev(seq_len(n)))
    smoothed_x <- if (smooth_with_ends) {
      t(apply(x, 1, function(row) smooth_with_ends_internal(row, pyramid)))
    } else {
      t(apply(x, 1, function(row) stats::convolve(row, pyramid, type = "filter"))) / sum(pyramid)
    }
    ## Downsample
    smoothed_x <- smoothed_x[, seq(1, ncol(smoothed_x), by = step), drop = FALSE]
  } else {
    ## Explicit preallocation to ensure stable orientation/dimensions (cells x genes)
    pyramid <- rep(1, n_genes)
    out <- matrix(NA_real_, nrow = n_cells, ncol = n_genes)
    for (i in seq_len(n_cells)) {
      out[i, ] <- stats::convolve(x[i, ], pyramid, type = "filter") / sum(pyramid)
    }
    ## If you want downsampling consistent with the other branch, uncomment next line:
    # out <- out[, seq(1, ncol(out), by = step), drop = FALSE]
    smoothed_x <- out
  }

  ## Add row names (cell names)
  if (!is.null(rownames(x))) rownames(smoothed_x) <- rownames(x)

  list(smoothed_x = smoothed_x)
}

