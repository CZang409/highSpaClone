#' Tumor region annotation using reference (normal) cells
#'
#' @param obj A highSpaClone object
#' @param ref Character vector of reference (normal) cell-type labels matching `obj@annotation$cell.label`.
#' @param lambda Numeric; spatial regularization strength passed to `run_iter()`.
#' @param K K=2.
#' @param max_iter Integer; maximum number of outer iterations.
#' @param min_iter Integer; minimum number of iterations before allowing early stop.
#' @param epsilon Numeric; convergence threshold on relative objective change.
#' @param seed Integer; random seed used in k-means initialization/updates.
#'
#' @return
#' The input `obj` with updated fields:
#' \itemize{
#'   \item `obj@cnv.data`: optimized CNV matrix (cells × bins).
#'   \item `obj@cluster`: `data.frame(cell.id, x, y, cell.label)` with `cell.label` ∈ {`"Tumor"`, `"Other"`}.
#' }
#'
#' @section Convergence and stopping:
#' Early stop if relative change < `epsilon` after `min_iter`, or if the objective increases.
#' If `max_iter` is reached without meeting the criterion, a message is printed.
#'
#' @examples
#' \dontrun{
#' obj <- FindTumor(
#'   obj,
#'   ref = c("B cells", "T cells"),
#'   lambda = 1/1000,
#'   K = 2,
#'   max_iter = 300,
#'   min_iter = 30
#' )
#' head(obj@cluster)
#' }
#'
#' @importFrom dplyr filter select pull mutate
#' @importFrom Matrix Matrix
#' @export

FindTumor <- function(
    obj,
    ref,
    lambda = 1/1000,
    K = 2,
    max_iter = 500,
    min_iter = 30,
    epsilon = 5e-04,
    seed = 12345678
){

  start_time <- Sys.time()

  cat("========================================\n")
  cat("Start Tumor Annotation...\n")
  cat("========================================\n")

  # Extract smoothed data
  cnv <- t(obj@smoothed.data)
  cnv <- as.matrix(cnv)
  cnv <- cnv+1e-6

  cat("CNV matrix dimensions: ", nrow(cnv), " x ", ncol(cnv), "\n")

  # Get annotation
  label <- obj@annotation

  # Filter normal reference cells
  cat("Reference cell types: ", paste(ref, collapse = ", "), "\n")
  norm <- label %>% filter(cell.label %in% ref)
  cat("Number of reference cells: ", nrow(norm), "\n")

  if(nrow(norm) == 0){
    stop("No reference cells found. Please check the ref parameter.")
  }

  # Calculate reference mean
  norm_mat <- cnv[, which(colnames(cnv) %in% norm$cell.id)]
  row_means <- rowMeans(norm_mat)

  # Get spatial location
  spatial_location <- obj@location
  spatial_location <- spatial_location[match(colnames(cnv), spatial_location$cell.id), ]

  # Construct B matrix
  B <- diag(x = row_means)
  colnames(B) <- row.names(cnv)
  row.names(B) <- row.names(cnv)

  # Calculate V matrix
  V <- sweep(cnv, 1, row_means, "/")
  V <- t(V)

  # Create adjacency matrix
  AMatrix <- createA(spatial_location)
  Bins <- colnames(V)
  cnv.df <- getConcated(V, spatial_location)

  # Initialize k-means clustering
  kmeans_result <- kmeansFunc_Initialize(cnv.df[, Bins], K)
  numCluster <- length(unique(kmeans_result$kmeans))
  cnv.df$kmeans_cluster <- as.numeric(kmeans_result$kmeans)
  TotalClusters <- unique(cnv.df$kmeans_cluster)[order(unique(cnv.df$kmeans_cluster))]

  # Convert to sparse matrices for efficiency
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  countMatrix <- Matrix::Matrix(cnv, sparse = TRUE)
  AMatrix <- Matrix::Matrix(AMatrix, sparse = TRUE)

  # ===== Step 3: First iteration =====
  ResList <- run_iter(
    Y = countMatrix,
    BIn = B,
    A = AMatrix,
    V = V,
    vecStr = cnv.df$kmeans_cluster,
    lambda = lambda
  )

  obj_old <- ResList$Obj

  # ===== Step 4: Iterative optimization =====
  centers_old <- NULL
  converged <- FALSE
  final_iter <- 0

  for(iter in 1:max_iter){

    if(iter %% 50 == 0 || iter == 1){
      cat("Iteration ", iter, " / ", max_iter, "\n")
    }

    # Update V matrix
    rownames(ResList$V) <- colnames(countMatrix)
    colnames(ResList$V) <- Bins
    V <- ResList$V
    rm(ResList)
    gc(verbose = FALSE)

    # Update concatenated matrix
    cnv.df <- getConcated(V, spatial_location)

    # Update k-means clustering
    if(iter == 1){
      centers_old <- NULL
    }

    kmeans_result <- kmeansFunc_Iter(
      cnv.df[, Bins],
      K,
      centers_old,
      seed = seed
    )

    numCluster <- length(unique(kmeans_result$kmeans))
    cnv.df$kmeans_cluster <- as.numeric(kmeans_result$kmeans)
    rm(kmeans_result)
    TotalClusters <- unique(cnv.df$kmeans_cluster)[order(unique(cnv.df$kmeans_cluster))]

    # Update centers
    centers_old <- getMu(Bins, TotalClusters, cnv.df)

    # Run iteration
    ResList <- run_iter(
      Y = countMatrix,
      BIn = B,
      A = AMatrix,
      V = V,
      vecStr = cnv.df$kmeans_cluster,
      lambda = lambda
    )

    obj_new <- ResList$Obj

    logicalObjective <- (obj_old - obj_new) * 2.0 / abs(obj_new + obj_old) < epsilon

    if(is.na(obj_new) || logicalObjective){
      if(iter >= min_iter){
        cat("\nConverged at iteration ", iter, "\n")
        converged <- TRUE
        final_iter <- iter
        break
      } else {
        obj_old <- obj_new
      }
    } else if(obj_new >= obj_old){
      cat("\nObjective increased. Stopping at iteration ", iter, "\n")
      final_iter <- iter
      break
    } else {
      obj_old <- obj_new
    }

    final_iter <- iter
  }

  if(!converged && final_iter == max_iter){
    cat("\nReached maximum iterations (", max_iter, ") without convergence.\n")
  }

  step1_time <- Sys.time()
  iteration_time <- difftime(step1_time, start_time, units = "mins")
  cat("\nIteration time: ", round(iteration_time, 2), " minutes\n")

  # ===== Step 5: Store results in object =====
  # Store CNV data
  rownames(ResList$V) <- colnames(countMatrix)
  colnames(ResList$V) <- Bins
  obj@cnv.data <- ResList$V

  # Store cluster assignments

  loc <- obj@location
  loc <- loc[match(rownames(cnv.df),loc$cell.id),]

  # assign tumor labels
  cnv_mat <- obj@cnv.data
  norm_mat=cnv_mat[which(rownames(cnv_mat) %in% norm$cell.id),]

  id1 <- dplyr::pull(dplyr::select(dplyr::filter(cnv.df, kmeans_cluster == 1), cellID), 1)
  id2 <- dplyr::pull(dplyr::select(dplyr::filter(cnv.df, kmeans_cluster == 2), cellID), 1)
  mat1 <- cnv_mat[id1,]
  mat2 <- cnv_mat[id2,]

  ##assign labels to two clusters
  #calculate the distance between the centroids of two clusters and that of the normal reference.
  #the cluster that is further away from the norm is annotated as tumor.
  centroids <- lapply(list(norm_mat, mat1, mat2), function(cluster) {
    colMeans(cluster, na.rm = TRUE)
  })
  centroids_matrix <- do.call(rbind, centroids)
  rownames(centroids_matrix) <- c("Norm", "Mat1", "Mat2")
  dist_matrix <- dist(centroids_matrix)
  tumor <- which.max(dist_matrix[1:2])

  cnv.df$tumor <- ifelse(cnv.df$kmeans_cluster == tumor, "Tumor", "Other")

  cluster_df <- data.frame(
    cell.id = rownames(cnv.df),
    x = loc$x,
    y = loc$y,
    cell.label = cnv.df$tumor,
    row.names = rownames(cnv.df)
  )

  obj@cluster <- cluster_df

  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")

  cat("\n========================================\n")
  cat("Tumor annotation completed!\n")
  cat("Total time: ", round(total_time, 2), " minutes\n")
  cat("========================================\n")

  return(obj)
}

#' Tumor subclone identification with spatial regularization
#'
#' @param obj a highSpaClone object
#' @param ref Character vector of reference (normal) labels.
#' @param ref.id Optional character vector of reference cell IDs.
#' @param tumor Optional character vector of tumor labels; if `NULL`, all non-reference cells are used.
#' @param tumor.id Optional vector of tumor cell IDs; if both `tumor` and `tumor.id` are `NULL`,
#'   all non-reference cells are used.
#' @param lambda Numeric; spatial regularization strength passed to `run_iter()`.
#' @param K Integer; number of subclones (clusters). If `NULL`, it should be set by the caller.
#' @param max_iter Integer; maximum number of iterations.
#' @param min_iter Integer; minimum iterations before early stop.
#' @param epsilon Numeric; convergence threshold on relative objective change.
#' @param seed Integer; random seed.
#'
#' @return
#' The input `obj` with:
#' \itemize{
#'   \item `@cnv.data`: optimized CNV matrix (tumor cells × bins).
#'   \item `@cluster`: `data.frame(cell.id, x, y, cell.label)` where `cell.label` ∈ {`"Clone 1"`, …, `"Clone K"`}.
#' }
#'
#' @section Convergence and stopping:
#' Same as `FindTumor()`: early stop on relative improvement < `epsilon` after `min_iter`,
#' stop on objective increase, or stop at `max_iter`.
#'
#' @examples
#' \dontrun{
#' obj <- FindClone(
#'   obj,
#'   ref = c("B cells", "T cells"),
#'   K = 3,
#'   lambda = 1,
#'   max_iter = 500,
#'   min_iter = 30
#' )
#' table(obj@cluster$cluster)
#' }
#'
#' @importFrom dplyr filter select pull mutate
#' @importFrom Matrix Matrix
#' @export
FindClone <- function(
    obj,
    ref,
    ref.id = NULL,
    tumor = NULL,
    tumor.id = NULL,
    lambda = NULL,
    K = NULL,
    max_iter = 500,
    min_iter = 30,
    epsilon = 5e-04,
    seed = 12345678
){
  start_time <- Sys.time()
  cat("========================================\n")
  cat("Running highSpaClone...\n")
  cat("========================================\n")

  # ===== Step 1: Data preparation =====
  cat("\n[Step 1] Preparing data...\n")

  # 1) CNV matrix: transpose to (bins x cells) and add small epsilon to avoid zeros
  cnv <- Matrix::t(obj@smoothed.data)
  cnv <- as.matrix(cnv)
  eps <- 1e-6
  cnv <- cnv + eps
  all_cells <- colnames(cnv)
  cat("CNV matrix dimensions (bins x cells): ", nrow(cnv), " x ", ncol(cnv), "\n")

  # 2) Annotation must contain cell.id and cell.label
  label <- obj@annotation
  if (!all(c("cell.id","cell.label") %in% colnames(label))) {
    stop("`obj@annotation` must contain columns: cell.id, cell.label")
  }

  # Reference set (by types and/or IDs), then intersect with existing cells
  ref_ids_from_label <- if (!missing(ref) && !is.null(ref)) {
    label$cell.id[label$cell.label %in% ref]
  } else character(0)

  ref_ids <- unique(c(ref_ids_from_label, ref.id))
  ref_ids <- intersect(ref_ids, all_cells)

  cat("Reference cell types (ref): ",
      if(length(ref)) paste(ref, collapse=", ") else "NULL", "\n")
  cat("Reference cell number: ", length(ref_ids), "\n")

  if (length(ref_ids) == 0) {
    stop("No reference cells found. Check `ref` / `ref.id` against annotation and smoothed.data colnames.")
  }

  # 3) Reference mean per gene
  norm_mat <- cnv[, ref_ids, drop = FALSE]
  row_means <- rowMeans(norm_mat)

  # 4) Tumor set (by types and/or IDs); default to all non-reference cells if not provided
  tumor_ids_from_label <- if (!is.null(tumor)) {
    label$cell.id[label$cell.label %in% tumor]
  } else character(0)

  tumor_ids <- unique(c(tumor_ids_from_label, tumor.id))
  tumor_ids <- intersect(tumor_ids, all_cells)

  if (length(tumor_ids) == 0) {
    default_tumor <- setdiff(all_cells, ref_ids)
    if (length(default_tumor) == 0) {
      stop("After excluding reference cells, no cells remain for tumor set.")
    }
    tumor_ids <- default_tumor
    cat("Tumor set not provided. Using all non-reference cells: ", length(tumor_ids), "\n")
  } else {
    cat("Tumor cell number: ", length(tumor_ids), "\n")
  }

  # 5) Keep tumor columns only for downstream clustering/optimization
  cnv <- cnv[, tumor_ids, drop = FALSE]

  # 6) Align spatial coordinates to tumor cells
  spatial_location <- obj@location
  if (!is.null(rownames(spatial_location)) && all(tumor_ids %in% rownames(spatial_location))) {
    spatial_location <- spatial_location[tumor_ids, , drop = FALSE]
  } else if ("id" %in% colnames(spatial_location)) {
    idx <- match(tumor_ids, spatial_location$id)
    if (anyNA(idx)) stop("Some tumor_ids are not found in `obj@location$id`.")
    spatial_location <- spatial_location[idx, , drop = FALSE]
    rownames(spatial_location) <- spatial_location$id
  } else {
    stop("`obj@location` must have rownames as cell IDs or an `id` column.")
  }

  # ===== Step 2: Initialize matrices =====
  cat("\n[Step 2] Initializing matrices...\n")

  # Construct diagonal B with reference means (bins x bins)
  B <- diag(x = row_means)
  colnames(B) <- rownames(cnv)  # gene names
  rownames(B) <- rownames(cnv)

  # Compute V = CNV / ref_mean (cells x bins)
  V <- sweep(cnv, 1, row_means, "/")
  V <- t(V)  # -> cells x bins (tumor subset)

  # Spatial adjacency
  AMatrix <- createA(spatial_location)

  # Feature table for clustering
  Bins <- colnames(V)  # gene names
  cnv.df <- getConcated(V, spatial_location)

  # k-means initialization
  kmeans_result <- kmeansFunc_Initialize(cnv.df[, Bins, drop = FALSE], K)
  cnv.df$kmeans_cluster <- as.numeric(kmeans_result$kmeans)
  TotalClusters <- sort(unique(cnv.df$kmeans_cluster))

  # Sparse matrices for efficiency
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  countMatrix <- Matrix::Matrix(cnv, sparse = TRUE)     # bins x tumor_cells
  AMatrix     <- Matrix::Matrix(AMatrix, sparse = TRUE)

  # Initial IRIS iteration
  ResList <- run_iter(
    Y = countMatrix,
    BIn = B,
    A = AMatrix,
    V = V,
    vecStr = cnv.df$kmeans_cluster,
    lambda = lambda
  )
  obj_old <- ResList$Obj

  # ===== Step 3: Iterative optimization =====
  cat("\n[Step 3] Starting iterative optimization...\n")
  centers_old <- NULL
  converged <- FALSE
  final_iter <- 0

  for (iter in 1:max_iter) {
    if (iter %% 25 == 0 || iter == 1) cat("Iteration ", iter, " / ", max_iter, "\n")

    # Update V from ResList
    rownames(ResList$V) <- colnames(countMatrix)  # cell IDs
    colnames(ResList$V) <- Bins                   # gene names
    V <- ResList$V
    rm(ResList); gc(verbose = FALSE)

    # Update feature table
    cnv.df <- getConcated(V, spatial_location)

    # k-means update
    if (iter == 1) centers_old <- NULL
    kmeans_result <- kmeansFunc_Iter(
      cnv.df[, Bins, drop = FALSE],
      K,
      centers_old,
      seed = seed
    )
    cnv.df$kmeans_cluster <- as.numeric(kmeans_result$kmeans)
    TotalClusters <- sort(unique(cnv.df$kmeans_cluster))

    # Update cluster centers
    centers_old <- getMu(Bins, TotalClusters, cnv.df)

    # iteration
    ResList <- run_iter(
      Y = countMatrix,
      BIn = B,
      A = AMatrix,
      V = V,
      vecStr = cnv.df$kmeans_cluster,
      lambda = lambda
    )
    obj_new <- ResList$Obj

    # Relative improvement check
    logicalObjective <- (obj_old - obj_new) * 2.0 / abs(obj_new + obj_old) < epsilon

    if (is.na(obj_new) || logicalObjective) {
      if (iter >= min_iter) {
        cat("\nConverged at iteration ", iter, "\n")
        converged <- TRUE
        final_iter <- iter
        break
      } else {
        obj_old <- obj_new
      }
    } else if (obj_new >= obj_old) {
      cat("\nObjective increased. Stopping at iteration ", iter, "\n")
      final_iter <- iter
      break
    } else {
      obj_old <- obj_new
    }
    final_iter <- iter
  }

  if (!converged && final_iter == max_iter) {
    cat("\nReached maximum iterations (", max_iter, ") without convergence.\n")
  }

  step1_time <- Sys.time()
  iteration_time <- difftime(step1_time, start_time, units = "mins")
  cat("\nIteration time: ", round(iteration_time, 2), " minutes\n")

  # ===== Step 4: Store results in object =====
  cat("\n[Step 4] Storing results...\n")

  # Write back CNV (cells x bins; tumor cells only)
  rownames(ResList$V) <- colnames(countMatrix)  # tumor cell IDs
  colnames(ResList$V) <- Bins
  obj@cnv.data <- ResList$V

  # Write back clustering (tumor cells only)
  loc <- obj@location
  if (!is.null(rownames(loc)) && all(rownames(ResList$V) %in% rownames(loc))) {
    loc <- loc[rownames(ResList$V), , drop = FALSE]
  } else if ("id" %in% colnames(loc)) {
    loc <- loc[match(rownames(ResList$V), loc$cell.id), , drop = FALSE]
    rownames(loc) <- loc$cell.id
  } else {
    stop("`obj@location` must have rownames as cell IDs or an `id` column.")
  }

  cluster_df <- data.frame(
    cell.id = rownames(cnv.df),
    x = loc$x,
    y = loc$y,
    cell.label = paste("Clone", cnv.df$kmeans_cluster),
    row.names = rownames(cnv.df),
    stringsAsFactors = FALSE
  )
  obj@cluster <- cluster_df

  cat("Results stored in object.\n")

  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")

  cat("\n========================================\n")
  cat("CNV inference and tumor subclone detection completed!\n")
  cat("Total time: ", round(total_time, 2), " minutes\n")
  cat("========================================\n")

  return(obj)
}

#' Suggest the number of subclones (K) using average silhouette
#'
#' @description
#' Computes a recommendation for the number of subclones \code{K}
#' using the average silhouette on an MDS embedding derived from CNV data.
#' Optionally fuses spatial distances from \code{obj@location} with CNV distances
#' (weight \code{alpha}) before embedding. Saves two figures and prints the
#' silhouette-vs-K plot to the device. Pick \code{K} that obtaions the largest average silhoutte score.
#'
#' @param obj a highSpaClone object
#' @param ref Character vector of **reference (normal)** cell labels.
#' @param ref.id Character vector of reference **cell IDs**.
#' @param tumor Character vector of **tumor** cell labels.
#' @param tumor.id Character vector of tumor **cell IDs**.
#'   If both tumor label/IDs are absent, all non-reference cells are used.
#' @param k_range Integer vector of K values to evaluate.
#' @param n_sub Integer. Max number of cells to subsample for evaluation (default: 10000).
#' @param alpha \code{NULL} for CNV-only. If a number in \code{[0,1]},
#'   fuse CNV and spatial distances: \eqn{d = alpha * d_{CNV} + (1-alpha) * d_{spatial}} (default: 0.7).
#' @param seed Integer seed for reproducibility.
#' @param out_dir Output directory for figures.
#'
#' @details
#' Two figures are saved to \code{out_dir}:
#' \itemize{
#'   \item \code{silhouette_vs_k.png}: Average silhouette vs \code{K}.
#'   \item \code{silhouette_bars_k_recommended.png}: Silhouette bars at the recommended \code{K}.
#' }
#'
#'
#' @seealso \code{\link[cluster]{silhouette}}, \code{\link[stats]{cmdscale}}, \code{\link[stats]{kmeans}}
#'
#' @importFrom stats dist cmdscale kmeans
#' @export
suggest_k <- function(
    obj,
    ref,
    ref.id = NULL,
    tumor,
    tumor.id = NULL,
    k_range = 2:8,
    n_sub = 10000,
    alpha = 0.7,
    seed = 123,
    out_dir = "figs_k"
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("cluster", quietly = TRUE)) stop("Package 'cluster' is required.")
  set.seed(seed)

  ## ---------- 0) Pull matrices ----------
  cnv_full <- t(obj@smoothed.data)
  cnv_full <- as.matrix(cnv_full) + 1e-6
  all_cells <- colnames(cnv_full)

  label <- obj@annotation
  if (!all(c("cell.id","cell.label") %in% colnames(label))) {
    stop("`obj@annotation` must contain columns: cell.id, cell.label")
  }

  ## ---------- 1) Build ref / tumor sets ----------
  ref_ids_from_label <- if (!missing(ref) && !is.null(ref)) label$cell.id[label$cell.label %in% ref] else character(0)
  ref_ids <- unique(c(ref_ids_from_label, ref.id))
  ref_ids <- intersect(ref_ids, all_cells)
  if (length(ref_ids) == 0) stop("No reference cells found. Check `ref` / `ref.id`.")

  tumor_ids_from_label <- if (!is.null(tumor)) label$cell.id[label$cell.label %in% tumor] else character(0)
  tumor_ids <- unique(c(tumor_ids_from_label, tumor.id))
  tumor_ids <- intersect(tumor_ids, all_cells)
  if (length(tumor_ids) == 0) {
    tumor_ids <- setdiff(all_cells, ref_ids)
    if (length(tumor_ids) == 0) stop("After excluding reference cells, no cells remain for tumor set.")
  }

  ## ---------- 2) Reference mean from FULL CNV (not tumor-only) ----------
  ref_mat  <- cnv_full[, ref_ids, drop = FALSE]    # bins x ref
  row_means <- rowMeans(ref_mat)                   # length = bins

  ## ---------- 3) Restrict to tumor cells, then compute V (cells x bins) ----------
  cnv_tumor <- cnv_full[, tumor_ids, drop = FALSE] # bins x tumor
  V <- t(sweep(cnv_tumor, 1, row_means, "/"))      # cells x bins

  ## ---------- 4) (Optional) Align location if alpha is provided ----------
  loc_use <- NULL
  if (!is.null(alpha)) {
    loc <- obj@location
    need_cols <- c("cell.id", "x", "y")
    if (!all(need_cols %in% colnames(loc))) stop("`obj@location` must contain columns: cell.id, x, y")
    idx_loc <- match(tumor_ids, loc$cell.id)
    if (anyNA(idx_loc)) {
      miss <- tumor_ids[is.na(idx_loc)]
      stop(sprintf("Some tumor IDs missing in obj@location$id (n=%d). e.g., %s",
                   length(miss), paste(head(miss, 5), collapse = ", ")))
    }
    loc_use <- loc[idx_loc, c("cell.id","x","y"), drop = FALSE]
    rownames(loc_use) <- loc_use$cell.id
  }

  ## ---------- 5) Drop columns with no variation / non-finite ----------
  finite_cols <- colSums(is.finite(V)) > 0
  if (!all(finite_cols)) V <- V[, finite_cols, drop = FALSE]

  has_variation <- apply(V, 2, function(x) {
    r <- range(x, na.rm = TRUE)
    is.finite(r[1]) && is.finite(r[2]) && (r[2] > r[1])
  })
  has_variation[is.na(has_variation)] <- FALSE

  if (!all(has_variation)) V <- V[, has_variation, drop = FALSE]
  if (ncol(V) == 0) stop("All columns became constant/NA after filtering.")

  ## ---------- 6) Subsample cells for evaluation ----------
  n_sub <- min(n_sub, nrow(V))
  idx_sub <- sample(seq_len(nrow(V)), n_sub)
  X_cnv  <- V[idx_sub, , drop = FALSE]
  if (!is.null(alpha)) {
    loc_sub <- loc_use[idx_sub, , drop = FALSE]
  } else {
    loc_sub <- NULL
  }

  ## ---------- 7) Distances: CNV-only or fused (CNV+spatial) ----------
  d_cnv <- stats::dist(X_cnv, method = "euclidean")
  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1) stop("`alpha` must be in [0,1] or NULL.")
    d_sp  <- stats::dist(loc_sub[, c("x","y")], method = "euclidean")
    rng1 <- range(d_cnv); d_cnv_n <- (d_cnv - rng1[1]) / (rng1[2] - rng1[1] + 1e-12)
    rng2 <- range(d_sp);  d_sp_n  <- (d_sp  - rng2[1]) / (rng2[2] - rng2[1]  + 1e-12)
    dmat <- alpha * d_cnv_n + (1 - alpha) * d_sp_n
  } else {
    dmat <- d_cnv
  }

  ## ---------- 8) Embed + evaluate silhouette ----------
  mds_embed <- stats::cmdscale(as.matrix(dmat), k = 10)
  ks <- as.integer(k_range)
  sil_mean <- numeric(length(ks))
  for (i in seq_along(ks)) {
    k  <- ks[i]
    km <- stats::kmeans(mds_embed, centers = k, nstart = 25, iter.max = 100)
    sil <- cluster::silhouette(km$cluster, dmat)
    sil_mean[i] <- mean(sil[, "sil_width"])
  }
  res_df <- data.frame(k = ks, mean_silhouette = sil_mean)

  ## ---------- 9) Choose K and print ----------
  k_reco <- ks[which.max(sil_mean)]
  cat(sprintf("Suggested K by average silhouette = %d (mean=%.3f)\n",
              k_reco, max(sil_mean, na.rm = TRUE)))

  ## ---------- 10) Save exactly two figures ----------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # 1) Average Silhouette vs k
  p_sil <- ggplot2::ggplot(res_df, ggplot2::aes(x = k, y = mean_silhouette)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3, shape = 21, fill = "white") +
    ggplot2::geom_vline(xintercept = k_reco, linetype = "dashed", linewidth = 0.8, color = "red") +
    ggplot2::labs(title = "Average Silhouette vs k", x = "k", y = "Mean silhouette") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
    )

  print(p_sil)
  ggplot2::ggsave(file.path(out_dir, "silhouette_vs_k.png"),
                  plot = p_sil, width = 7, height = 5, dpi = 300)

  # 2) Bars at recommended K
  km_reco <- stats::kmeans(mds_embed, centers = k_reco, nstart = 50, iter.max = 200)
  sil_rec <- cluster::silhouette(km_reco$cluster, dmat)
  grDevices::png(file.path(out_dir, "silhouette_bars_k_recommended.png"),
                 width = 2000, height = 1500, res = 250)
  graphics::plot(
    sil_rec,
    main = sprintf("Silhouette bars (k=%d), mean=%.3f", k_reco, mean(sil_rec[, 3])),
    col = "gray", border = NA
  )
  grDevices::dev.off()

  invisible(list(k_recommended = k_reco,
                 summary_df    = res_df,
                 out_dir       = out_dir))
}

################################################################################################
###internal functions
createA <- function(location) {
  location <- as.data.frame(location)
  norm_cords <- location[, c("x", "y")]

  # normalize to [0, 1] range
  norm_cords$x <- norm_cords$x - min(norm_cords$x)
  norm_cords$y <- norm_cords$y - min(norm_cords$y)
  scaleFactor <- max(norm_cords$x, norm_cords$y)
  norm_cords$x <- norm_cords$x / scaleFactor
  norm_cords$y <- norm_cords$y / scaleFactor
  rownames(norm_cords) <- rownames(location)

  # find 10 nearest neighbors (excluding self)
  ineibor <- 11
  near_data <- RANN::nn2(norm_cords[, 1:2], k = ineibor)
  neighbors <- near_data$nn.idx[, -1]

  # build sparse adjacency matrix
  Nmat <- Matrix::Matrix(0, nrow = nrow(neighbors), ncol = nrow(neighbors), sparse = TRUE)
  for (icol in 1:ncol(neighbors)) {
    edges <- data.frame(i = 1:nrow(neighbors), j = neighbors[, icol])
    adjacency <- Matrix::sparseMatrix(
      i = as.integer(edges$i),
      j = as.integer(edges$j),
      x = 1,
      dims = rep(nrow(neighbors), 2),
      use.last.ij = TRUE
    )
    Nmat <- Nmat + adjacency
  }

  # keep only mutual neighbors
  Nmat <- Nmat * Matrix::t(Nmat)
  rownames(Nmat) <- colnames(Nmat) <- rownames(norm_cords)

  return(Nmat)
}

getConcated <- function(V, location) {
  V <- as.data.frame(V)
  V$cellID <- rownames(V)
  V$x <- location$x
  V$y <- location$y

  return(V)
}

getMu <- function(Bins, TotalClusters, data) {

  mu <- matrix(0, nrow = length(Bins), ncol = length(TotalClusters))
  colnames(mu) <- TotalClusters
  rownames(mu) <- Bins

  Clusters <- unique(data$kmeans_cluster)
  Clusters <- Clusters[order(Clusters)]

  for (i in Clusters) {
    mu[, colnames(mu) == i] <- colMeans(data[data$kmeans_cluster == i, Bins])
  }

  return(mu)
}


kmeansFunc_Initialize <- function(data,k,seed=12345678){
  set.seed(seed)
  if(nrow(data)  < 300000){
    numStart = 100
  }else{
    numStart = 1
  }
  cl <- suppressWarnings(try(kmeans(data, k,nstart = numStart,iter.max=100),silent = TRUE))

  return(list(kmeans =cl$cluster))
}


kmeansFunc_Iter <- function(data,k,centers_old = NULL,seed=12345678){
  set.seed(seed)
  if(nrow(data)  < 300000){
    numStart = 10
  }else{
    numStart = 1
  }
  if(is.null(centers_old)){
    cl <- suppressWarnings(try(kmeans(data, k,nstart = 1,iter.max=numStart),silent = TRUE))
  }else{
    cl <- suppressWarnings(try(kmeans(data, centers_old,nstart = 2,iter.max=numStart),silent = TRUE))
  }
  if(class(cl) == "try-error"){
    cl <- suppressWarnings(try(kmeans(data, k,nstart = 1,iter.max=numStart),silent = TRUE))
  }

  if (inherits(cl, "kmeans")) {
    resList <- list(kmeans = cl$cluster, centers = cl$centers)
  } else {
    resList <- list(kmeans = NULL, centers = NULL)
    warning("K-means clustering failed. Returning NULL values.")
  }
  return(resList)
}

