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
  
  # Check if sparse and keep it sparse
  is_sparse <- inherits(expr, "sparseMatrix")
  if (is_sparse) {
    expr <- as(expr, "CsparseMatrix")
    cat("Working with sparse matrix format (keeping sparse)\n")
  }
  
  var <- obj@gene_order
  n_cells <- ncol(expr)
  cat("Total cells: ", n_cells, "\n")
  
  # Decide whether to use chunking
  if (use_chunk && n_cells > chunk_size) {
    # ===== CHUNKING MODE =====
    cat("\n=== Chunking mode enabled ===\n")
    cat("Processing in chunks of size ", chunk_size, "...\n")
    
    # Split cells into chunks
    n_chunks <- ceiling(n_cells / chunk_size)
    chunk_indices <- split(1:n_cells, ceiling(seq_along(1:n_cells) / chunk_size))
    
    cat("Total chunks: ", n_chunks, "\n")
    
    # Choose processing method based on parallel parameter
    if (parallel) {
      cat("=== Parallel processing enabled ===\n")
      
      # Check platform
      is_unix <- .Platform$OS.type == "unix"
      
      # Determine number of cores to use
      if (is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1
        n_cores <- max(1, n_cores)
      }
      
      # Apply max_cores limit
      n_cores <- min(n_cores, max_cores, n_chunks)
      
      cat("Requested cores: ", n_cores, " (max_cores limit: ", max_cores, ")\n")
      
      if (is_unix) {
        # Unix/Linux/Mac: use mclapply
        cat("Unix system detected: using mclapply\n")
        cat("Using ", n_cores, " cores\n")
        
        # Estimate memory
        expr_size_gb <- as.numeric(object.size(expr)) / (1024^3)
        estimated_memory_gb <- expr_size_gb * n_cores * 1.5
        cat("Estimated memory usage: ~", round(estimated_memory_gb, 2), " GB\n")
        
        if (estimated_memory_gb > 50) {
          warning("High memory usage estimated (", round(estimated_memory_gb, 2), 
                  " GB). Consider reducing n_cores or max_cores, or using parallel=FALSE.")
        }
        
        if (!requireNamespace("parallel", quietly = TRUE)) {
          stop("Package 'parallel' is required for parallel processing.")
        }
        
        library(parallel)
        chunk_results <- mclapply(seq_along(chunk_indices), function(i) {
          if (i %% 5 == 0) cat("Processing chunk", i, "of", n_chunks, "\n")
          indices <- chunk_indices[[i]]
          expr_chunk <- expr[, indices, drop = FALSE]
          results <- infercnv_chunk(expr_chunk, var, exclude_chromosomes, 
                                    window_size, step, smooth_with_ends)
          return(results$x_smoothed)
        }, mc.cores = n_cores)
        
      } else {
        # Windows: use parLapply
        cat("Windows system detected: using parLapply\n")
        cat("Using ", n_cores, " cores\n")
        
        # Estimate memory
        expr_size_gb <- as.numeric(object.size(expr)) / (1024^3)
        estimated_memory_gb <- expr_size_gb * n_cores * 1.5
        cat("Estimated memory usage: ~", round(estimated_memory_gb, 2), " GB\n")
        
        if (estimated_memory_gb > 50) {
          warning("High memory usage estimated (", round(estimated_memory_gb, 2), 
                  " GB). Consider reducing n_cores or max_cores, or using parallel=FALSE.")
        }
        
        if (!requireNamespace("parallel", quietly = TRUE)) {
          stop("Package 'parallel' is required for parallel processing.")
        }
        
        # Create cluster
        cl <- parallel::makeCluster(n_cores)
        
        # Export necessary objects
        parallel::clusterExport(cl, 
                                c("infercnv_chunk", "running_mean_by_chromosome", 
                                  "running_mean_for_chromosome", "running_mean",
                                  "smooth_with_ends_internal", "sort_chromosomes",
                                  "var", "exclude_chromosomes", "window_size", 
                                  "step", "smooth_with_ends", "chunk_indices"),
                                envir = environment())
        
        # Load required packages on each worker
        parallel::clusterEvalQ(cl, {
          library(dplyr)
          library(Matrix)
        })
        
        # Process chunks
        chunk_results <- tryCatch({
          parallel::parLapply(cl, seq_along(chunk_indices), function(i) {
            if (i %% 5 == 0) cat("Processing chunk", i, "of", length(chunk_indices), "\n")
            indices <- chunk_indices[[i]]
            expr_chunk <- expr[, indices, drop = FALSE]
            results <- infercnv_chunk(expr_chunk, var, exclude_chromosomes, 
                                      window_size, step, smooth_with_ends)
            return(results$x_smoothed)
          })
        }, finally = {
          parallel::stopCluster(cl)
        })
      }
      
      cat("Parallel processing completed.\n")
      
    } else {
      # Sequential processing with chunks
      cat("=== Sequential chunk processing ===\n")
      cat("Processing ", n_chunks, " chunks sequentially...\n")
      
      chunk_results <- lapply(seq_along(chunk_indices), function(i) {
        if (i %% 5 == 0 || i == 1 || i == n_chunks) {
          cat("Processing chunk", i, "of", n_chunks, 
              sprintf("(%.1f%% complete)\n", 100*i/n_chunks))
        }
        indices <- chunk_indices[[i]]
        expr_chunk <- expr[, indices, drop = FALSE]
        results <- infercnv_chunk(expr_chunk, var, exclude_chromosomes, 
                                  window_size, step, smooth_with_ends)
        
        # Free memory
        rm(expr_chunk)
        if (i %% 10 == 0) gc(verbose = FALSE)
        
        return(results$x_smoothed)
      })
      
      cat("Sequential chunk processing completed.\n")
    }
    
    # Combine results - use sparse-friendly method if applicable
    cat("Combining results...\n")
    
    # Check if results are sparse
    if (inherits(chunk_results[[1]], "sparseMatrix")) {
      cat("Combining sparse matrices...\n")
      # Use Matrix::rBind to keep sparse format
      res <- do.call(Matrix::rBind, chunk_results)
    } else {
      cat("Combining dense matrices...\n")
      res <- do.call(rbind, chunk_results)
    }
    
    cat("Result class: ", class(res)[1], "\n")
    cat("Result dimensions: ", dim(res)[1], " x ", dim(res)[2], "\n")
    
    # Free memory
    rm(chunk_results)
    gc(verbose = FALSE)
    
    # Get chr_pos from first chunk
    cat("Calculating chromosome positions...\n")
    chr_pos_result <- infercnv_chunk(expr[, 1:min(100, n_cells), drop = FALSE], 
                                     var, exclude_chromosomes, 
                                     window_size, step, smooth_with_ends)
    chr_pos <- chr_pos_result$chr_pos
    
  } else {
    # ===== NO CHUNKING MODE =====
    if (!use_chunk) {
      cat("\n=== No chunking mode (use_chunk=FALSE) ===\n")
    } else {
      cat("\n=== Cell count below chunk size ===\n")
    }
    cat("Processing all ", n_cells, " cells at once...\n")
    
    # Estimate memory for full processing
    expr_size_gb <- as.numeric(object.size(expr)) / (1024^3)
    cat("Expression matrix size: ~", round(expr_size_gb, 2), " GB\n")
    
    if (expr_size_gb > 10) {
      warning("Large matrix detected (", round(expr_size_gb, 2), 
              " GB). Consider using use_chunk=TRUE to reduce memory usage.")
    }
    
    # Process all cells at once
    cat("Running smoothing on all cells...\n")
    results <- infercnv_chunk(expr, var, exclude_chromosomes, 
                              window_size, step, smooth_with_ends)
    chr_pos <- results$chr_pos
    res <- results$x_smoothed
    
    cat("Result class: ", class(res)[1], "\n")
    cat("Result dimensions: ", dim(res)[1], " x ", dim(res)[2], "\n")
    
    cat("No-chunk processing completed.\n")
  }
  
  # Store results - keep the original format
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

running_mean_by_chromosome <- function(expr, var, exclude_chromosomes, window_size, step, smooth_with_ends) {
  
  if (!is.null(exclude_chromosomes)) {
    var <- var[!(var$chromosome %in% exclude_chromosomes),]
  }
  
  chromosomes <- unique(var$chromosome)
  chromosomes <- sort_chromosomes(chromosomes)
  
  
  running_means <- lapply(chromosomes, function(chr) {
    running_mean_for_chromosome(chr, expr, var, window_size, step, smooth_with_ends)
  })
  
  names(running_means) <- chromosomes
  running_means_list <- lapply(running_means, `[[`, 1)
  
  chr_start_pos <- list()
  chr_pos_offset <- 0
  for (i in seq_along(chromosomes)) {
    chr_start_pos[[chromosomes[i]]] <- chr_pos_offset
    chr_pos_offset <- chr_pos_offset + ncol(running_means_list[[i]])
  }
  names(chr_start_pos) <- chromosomes
  
  x_smoothed <- do.call(cbind, running_means_list)
  
  name_vectors <- c()
  for (list_name in names(running_means_list)) {
    current_matrix <- running_means_list[[list_name]]
    n <- ncol(current_matrix)
    name_vector <- paste(list_name, 1:n, sep = "_")
    name_vectors <- c(name_vectors,name_vector)
  }
  
  colnames(x_smoothed) <- name_vectors
  
  return(list(chr_pos = chr_start_pos, x_smoothed = x_smoothed))
}

running_mean_for_chromosome <- function(chr, expr, var, window_size, step, smooth_with_ends) {
  
  genes <- var %>% filter(chromosome == chr) %>% arrange(start) %>% rownames()
  expr_chr <- expr[which(rownames(expr) %in% genes),]
  running_mean_res <- running_mean(expr_chr, window_size, step, smooth_with_ends)
  
  return(running_mean_res)
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
  
  x <- t(x)
  
  if (n < 2) {
    print("window length < 2, returning original unmodified data")
    return(x)
    
  } else if (n < ncol(x)) {
    
    pyramid <- pmin(1:n, rev(1:n))
    
    smoothed_x <- if (smooth_with_ends) {
      Matrix::t(apply(x, 1, function(row) {
        smooth_with_ends_internal(row, pyramid)
      }))
    } else {
      Matrix::t(apply(x, 1, function(row) {
        stats::convolve(row, pyramid, type = "filter")
      })) / sum(pyramid)
    }
    
    smoothed_x <- smoothed_x[, seq(1, ncol(smoothed_x), by = step), drop = FALSE]
    return(list(smoothed_x = smoothed_x))
    
  } else {
    
    pyramid <- rep(1, ncol(x))
    smoothed_x <- apply(x, 1, function(row) {
      stats::convolve(row, pyramid, type = "filter")
    }) / sum(pyramid)
    
    smoothed_x <- as.matrix(smoothed_x)
    return(list(smoothed_x = smoothed_x))
  }
}


