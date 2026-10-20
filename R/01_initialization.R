# =========================
# Helpers
# =========================

# helper: rescale to [0,1]
.scale01 <- function(x) {
  r <- range(x, finite = TRUE)
  if (!is.finite(r[1]) || r[1] == r[2]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

# helper: orient sign using a reference vector (e.g., PC1)
.orient_by_ref <- function(t, ref) {
  if (length(t) != length(ref)) return(t)
  if (stats::cor(t, ref, use = "pairwise.complete.obs") < 0) -t else t
}

# helper: get PC1 quickly (prcomp is fine for now)
.pc1 <- function(X, center = TRUE, scale. = FALSE) {
  stats::prcomp(X, center = center, scale. = scale.)$x[, 1]
}

# =========================
# Ordering methods
# =========================

# Internal numerical implementation of fiedler_ordering.
.fiedler_ordering <- function(
    X,
    k = 15,
    weight = c("binary", "rbf", "inv"),
    sigma = NULL,
    keep = c("giant", "all"),
    return_full = TRUE
) {
  X <- as.matrix(X)
  n0 <- nrow(X)
  k_missing <- missing(k)
  if (n0 < 3) stop("Need n >= 3.")
  if (k_missing) {
    k <- min(as.integer(k), n0 - 1L)
  }
  if (k < 2 || k >= n0) stop("Require 2 <= k < n.")

  weight <- match.arg(weight)
  keep <- match.arg(keep)

  build_knn_graph <- function(k_use) {
    nn <- RANN::nn2(data = X, query = X, k = k_use + 1L, treetype = "kd")
    idx <- nn$nn.idx[, -1, drop = FALSE]
    dist <- nn$nn.dists[, -1, drop = FALSE]

    i <- rep(seq_len(n0), each = k_use)
    j <- as.vector(t(idx))
    d <- as.vector(t(dist))

    sigma_use <- sigma
    w <- switch(
      weight,
      binary = rep(1, length(d)),
      inv = 1 / pmax(d, 1e-12),
      rbf = {
        if (is.null(sigma_use)) sigma_use <- stats::median(d)
        exp(-(d^2) / (2 * sigma_use^2))
      }
    )

    W <- Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(n0, n0))
    W <- pmax(W, Matrix::t(W))
    g <- igraph::graph_from_adjacency_matrix(W > 0, mode = "undirected", diag = FALSE)

    list(
      W = W,
      comp = igraph::components(g)
    )
  }

  graph_fit <- build_knn_graph(k)

  if (k_missing && graph_fit$comp$no > 1L) {
    lower_k <- k
    upper_k <- k

    while (graph_fit$comp$no > 1L && upper_k < (n0 - 1L)) {
      lower_k <- upper_k
      upper_k <- min(
        n0 - 1L,
        max(upper_k + 1L, as.integer(ceiling(1.5 * upper_k)))
      )
      graph_fit <- build_knn_graph(upper_k)
    }

    if (graph_fit$comp$no == 1L && upper_k > lower_k) {
      best_fit <- graph_fit
      left <- lower_k
      right <- upper_k
      while ((right - left) > 1L) {
        mid <- left + (right - left) %/% 2L
        candidate <- build_knn_graph(mid)
        if (candidate$comp$no == 1L) {
          right <- mid
          best_fit <- candidate
        } else {
          left <- mid
        }
      }
      k <- right
      graph_fit <- best_fit
    } else {
      k <- upper_k
    }
  }

  W <- graph_fit$W
  comp <- graph_fit$comp
  keep_idx <- seq_len(n0)

  if (keep == "giant" && comp$no > 1) {
    giant <- which.max(comp$csize)
    keep_idx <- which(comp$membership == giant)
    Wk <- W[keep_idx, keep_idx, drop = FALSE]
  } else {
    Wk <- W
  }

  if (nrow(Wk) < 3) {
    if (return_full) {
      return(list(
        t = rep(NA_real_, n0),
        keep_idx = keep_idx,
        n_components = comp$no,
        k_used = k
      ))
    }
    stop("Largest connected component has < 3 nodes; increase k or filter outliers.")
  }

  deg <- Matrix::rowSums(Wk)
  D_inv_sqrt <- Matrix::Diagonal(nrow(Wk), x = 1 / sqrt(pmax(deg, 1e-12)))
  S <- D_inv_sqrt %*% Wk %*% D_inv_sqrt
  S <- (S + Matrix::t(S)) / 2

  e <- RSpectra::eigs_sym(
    S, k = 2, which = "LA",
    opts = list(tol = 1e-10, maxitr = 20000, retvec = TRUE)
  )

  u2 <- e$vectors[, 2]
  t_kept <- .scale01(u2)

  if (return_full) {
    t_full <- rep(NA_real_, n0)
    t_full[keep_idx] <- t_kept
    t_out <- t_full
  } else {
    t_out <- t_kept
  }

  # if comp$no is larger than 1, warn the user
  if (comp$no > 1) {
    warning("The kNN graph has ", comp$no, " connected components; ",
            "consider increasing k or filtering outliers.")
  }

  list(
    t = t_out,
    keep_idx = keep_idx,
    n_components = comp$no,
    k_used = k
  )
}


# Internal numerical implementation of tSNE_ordering.
.tSNE_ordering <- function(
    X,
    tSNE_dims = 1,
    component = 1,
    perplexity = 10,
    max_iter = 500,
    seed = NULL,
    scale01 = TRUE,
    orient_by_pc1 = TRUE
) {
  X <- as.matrix(X)
  n <- nrow(X)

  if (component < 1 || component > tSNE_dims) {
    stop("component must be between 1 and tSNE_dims.")
  }
  if (perplexity >= (n - 1) / 3) {
    stop("perplexity must be < (n-1)/3 for Rtsne (reduce perplexity or increase n).")
  }
  if (!is.null(seed)) set.seed(seed)

  tsne_result <- Rtsne::Rtsne(
    X = X, dims = tSNE_dims, perplexity = perplexity,
    verbose = FALSE, max_iter = max_iter
  )

  t <- tsne_result$Y[, component]

  if (orient_by_pc1) t <- .orient_by_ref(t, .pc1(X))
  if (scale01) t <- .scale01(t)

  list(t = t, keep_idx = seq_len(n))
}


#' Principal-component ordering
#'
#' @param X Numeric matrix (n x D). Missing values are filled with each
#'   feature's observed mean, with a warning.
#' @param center logical.
#' @param scale logical.
#' @param scale01 logical; if TRUE, rescale returned t to [0,1].
#' @param component Integer principal-component index to use as the ordering
#'   score. Defaults to 1.
#' @return list(t, keep_idx)
#' @export
PCA_ordering <- function(
    X,
    center = TRUE,
    scale = FALSE,
    scale01 = TRUE,
    component = 1L
) {
  X <- .mpcurve_ordering_input(X, "PCA_ordering")
  component <- .mpcurve_integer_setting(component, "component")
  max_comp <- min(nrow(X), ncol(X))
  if (component > max_comp)
    stop(sprintf("component %d exceeds available PCs (%d).", component, max_comp))

  if (component == 1L) {
    t <- .pc1(X, center = center, scale. = scale)
  } else {
    pca <- stats::prcomp(X, center = center, scale. = scale)
    t <- pca$x[, component]
  }
  if (scale01) t <- .scale01(t)
  list(t = t, keep_idx = seq_len(nrow(X)))
}


# Internal numerical implementation of pcurve_ordering.
.pcurve_ordering <- function(
    X,
    smoother = c("smooth_spline", "lowess", "periodic_lowess"),
    thresh = 0.001,
    maxit = 10,
    stretch = 2,
    approx_points = FALSE,
    scale01 = TRUE,
    orient_by_pc1 = TRUE
) {
  X <- as.matrix(X)
  smoother <- match.arg(smoother)
  pc_fit <- princurve::principal_curve(
    X,
    smoother = smoother,
    thresh = thresh,
    stretch = stretch,
    approx_points = approx_points,
    maxit = maxit
  )

  t <- pc_fit$lambda
  if (orient_by_pc1) t <- .orient_by_ref(t, .pc1(X))
  if (scale01) t <- .scale01(t)

  list(t = t, keep_idx = seq_len(nrow(X)), fit = pc_fit)
}

# Exact neighbors excluding the query sample. Break distance ties by row index
# so graphs are nested as k increases, including when observations coincide.
.isomap_neighbors <- function(X, k) {
  n <- nrow(X)
  query_count <- min(n, k + 2L)
  repeat {
    nn <- RANN::nn2(X, X, k = query_count, treetype = "kd", eps = 0)
    idx <- matrix(NA_integer_, n, k)
    distances <- matrix(NA_real_, n, k)
    boundary_tie <- FALSE
    for (row in seq_len(n)) {
      candidates <- which(nn$nn.idx[row, ] != row)
      ranked <- candidates[order(nn$nn.dists[row, candidates],
                                 nn$nn.idx[row, candidates])]
      chosen <- ranked[seq_len(k)]
      idx[row, ] <- nn$nn.idx[row, chosen]
      distances[row, ] <- nn$nn.dists[row, chosen]
      if (query_count < n && nn$nn.dists[row, query_count] <= distances[row, k]) {
        boundary_tie <- TRUE
      }
    }
    if (!boundary_tie) return(list(idx = idx, distances = distances))
    # Include every tie at the kth distance before applying the row-index rule.
    query_count <- min(n, 2L * query_count)
  }
}

.isomap_graph <- function(X, k) {
  nn <- .isomap_neighbors(X, k)
  edges <- data.frame(from = rep(seq_len(nrow(X)), each = k),
                      to = as.vector(t(nn$idx)),
                      weight = as.vector(t(nn$distances)))
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                        vertices = seq_len(nrow(X)))
  list(graph = graph, components = igraph::components(graph), k = k)
}

# Connectivity is monotone for nested union-kNN graphs. Bracket the first
# connected graph by doubling k, then find its exact minimum by binary search.
# Only the selected graph is embedded; no candidate MPCurve fits are needed.
.isomap_connected_graph <- function(X) {
  lower <- 0L
  upper <- 1L
  best <- .isomap_graph(X, upper)
  while (best$components$no > 1L) {
    lower <- upper
    upper <- min(nrow(X) - 1L, 2L * upper)
    best <- .isomap_graph(X, upper)
  }
  while (upper - lower > 1L) {
    mid <- lower + (upper - lower) %/% 2L
    candidate <- .isomap_graph(X, mid)
    if (candidate$components$no == 1L) {
      upper <- mid
      best <- candidate
    } else {
      lower <- mid
    }
  }
  best
}

# Internal numerical implementation of isomap_ordering.
.isomap_ordering <- function(
    X,
    k = NULL,
    ndim = 1,
    component = 1,
    landmark = NULL,
    landmark_method = c("random", "kmeans"),
    seed = NULL,
    keep = c("giant", "all"),
    return_full = TRUE,
    orient_by_pc1 = TRUE,
    scale01 = TRUE,
    eps = 1e-8
) {
  X <- as.matrix(X)
  n0 <- nrow(X)

  if (n0 < 3) stop("Need n >= 3.")
  if (!is.null(k) && (k < 1 || k >= n0)) stop("Require 1 <= k < n.")
  if (ndim < 1) stop("ndim must be >= 1.")
  if (component < 1 || component > ndim) stop("component must be between 1 and ndim.")

  keep <- match.arg(keep)
  landmark_method <- match.arg(landmark_method)

  if (!is.null(seed)) set.seed(seed)

  # --- kNN graph (weighted by Euclidean distances) ---
  graph_fit <- if (is.null(k)) .isomap_connected_graph(X) else .isomap_graph(X, k)
  k <- graph_fit$k
  g <- graph_fit$graph

  # --- connected components handling ---
  comp <- graph_fit$components
  keep_idx <- seq_len(n0)

  if (keep == "giant" && comp$no > 1) {
    giant <- which.max(comp$csize)
    keep_idx <- which(comp$membership == giant)
  }

  if (length(keep_idx) < 3) {
    if (return_full) {
      return(list(
        t = rep(NA_real_, n0),
        keep_idx = keep_idx,
        n_components = comp$no,
        k_used = k,
        embed = NULL,
        landmark_idx = integer(0),
        geodesic_to_landmark = NULL
      ))
    }
    stop("Largest connected component has < 3 nodes; increase k or filter outliers.")
  }

  # induce subgraph if needed
  if (length(keep_idx) < n0) {
    gk <- igraph::induced_subgraph(g, vids = keep_idx)
    Xk <- X[keep_idx, , drop = FALSE]
  } else {
    gk <- g
    Xk <- X
  }

  n <- nrow(Xk)

  # --- landmark count (always safe) ---
  if (is.null(landmark)) landmark <- min(1000L, n)
  landmark <- as.integer(landmark)
  landmark <- max(3L, min(landmark, n))  # critical: cap to n

  # --- choose landmarks (default random to avoid kmeans-centers pitfalls) ---
  if (landmark_method == "random") {
    Lk <- sort(sample.int(n, landmark))
  } else {
    # extra safety: centers must be in [1, n]
    centers <- max(1L, min(landmark, n))

    # kmeans can still error if X has NA/Inf; let it fall back cleanly
    Lk <- tryCatch({
      km <- stats::kmeans(Xk, centers = centers, iter.max = 30)

      # nearest observed point to each center
      nnc <- RANN::nn2(data = Xk, query = km$centers, k = 1, treetype = "kd")
      L0 <- unique(as.integer(nnc$nn.idx[, 1]))

      # if duplicates reduced count, top up randomly
      if (length(L0) < 3) L0 <- sort(sample.int(n, min(centers, n)))
      if (length(L0) < centers) {
        extra <- setdiff(seq_len(n), L0)
        if (length(extra) > 0) {
          add <- sample(extra, min(centers - length(L0), length(extra)))
          L0 <- sort(c(L0, add))
        }
      }
      sort(L0)
    }, error = function(e) {
      # fall back to random landmarks
      sort(sample.int(n, landmark))
    })
  }

  # --- geodesic distances ---
  D_to_L <- igraph::distances(
    gk,
    v = seq_len(n),
    to = Lk,
    weights = igraph::E(gk)$weight
  )
  D_LL <- D_to_L[Lk, , drop = FALSE]

  if (any(!is.finite(D_LL))) {
    if (keep == "all") {
      warning("kNN graph appears disconnected: Inf geodesic distances found. ",
              "Consider keep='giant' or increase k.")
    }
    # For embedding, replace Inf with a large finite value.
    finite_max <- suppressWarnings(max(D_LL[is.finite(D_LL)]))
    if (!is.finite(finite_max)) finite_max <- 1
    D_LL2 <- D_LL
    D_LL2[!is.finite(D_LL2)] <- 2 * finite_max
    YL <- stats::cmdscale(as.dist(D_LL2), k = ndim, eig = FALSE)
  } else {
    YL <- stats::cmdscale(as.dist(D_LL), k = ndim, eig = FALSE)
  }

  YL <- as.matrix(YL)
  if (ncol(YL) < ndim) {
    YL <- cbind(YL, matrix(0, nrow = nrow(YL), ncol = ndim - ncol(YL)))
  }

  # --- out-of-sample projection by inverse-distance weighting ---
  W <- 1 / (pmax(D_to_L, 0) + eps)
  W[!is.finite(W)] <- 0
  denom <- rowSums(W)
  denom[denom <= 0] <- 1

  Y <- (W %*% YL) / denom

  t_kept <- Y[, component]

  if (orient_by_pc1) t_kept <- .orient_by_ref(t_kept, .pc1(Xk))
  if (scale01) t_kept <- .scale01(t_kept)

  if (return_full) {
    t_full <- rep(NA_real_, n0)
    t_full[keep_idx] <- t_kept

    Y_full <- matrix(NA_real_, nrow = n0, ncol = ndim)
    Y_full[keep_idx, ] <- Y

    landmark_idx_full <- keep_idx[Lk]

    t_out <- t_full
    embed_out <- Y_full
    landmark_idx_out <- landmark_idx_full
  } else {
    t_out <- t_kept
    embed_out <- Y
    landmark_idx_out <- Lk
  }

  if (comp$no > 1) {
    warning("The kNN graph has ", comp$no, " connected components; ",
            "Isomap computed on ", length(keep_idx), " nodes (keep='", keep, "').")
  }

  list(
    t = t_out,
    keep_idx = keep_idx,
    n_components = comp$no,
    k_used = k,
    embed = embed_out,
    landmark_idx = landmark_idx_out,
    geodesic_to_landmark = D_to_L
  )
}


# Random ordering start: a uniformly random permutation of the samples, scaled
# to [0, 1]. Binned like any other ordering score, it assigns samples to
# position bins at random; quantile bins then hold floor(n/K) or ceil(n/K)
# samples each. Reproducible through set.seed().
.random_ordering <- function(X) {
  n <- nrow(as.matrix(X))
  list(t = .scale01(sample.int(n)), keep_idx = seq_len(n))
}


# =========================
# Ordering initialization
# =========================

#' Cut an ordering score into position bins
#'
#' @description
#' Discretizes \code{ordering_vec} into \code{K} ordered position bins and
#' summarizes the samples in each bin. These summaries start every
#' ordering-based CAVI fit.
#'
#' @details
#' \code{"quantile"} and \code{"kmeans"} bins are non-empty. Quantile cuts that
#' coincide because of tied scores reduce the number of bins, with a warning.
#' \code{"equal"} uses equal-width bins, which may be empty; an empty bin takes
#' its mean by linear interpolation between the nearest non-empty bins, or
#' from the single nearest non-empty bin at either end.
#'
#' @param X Numeric matrix (\code{n x d}).
#' @param ordering_vec Numeric ordering score of length \code{n}; \code{NA}
#'   scores are handled by \code{na_action}.
#' @param K Integer >= 2; requested number of bins.
#' @param discretization One of \code{"quantile"}, \code{"equal"}, \code{"kmeans"}.
#' @param na_action \code{"drop"} removes samples with \code{NA} scores;
#'   \code{"error"} stops.
#' @param eps Positive floor for the pooled variances.
#'
#' @return A list with \code{pi} (bin proportions), \code{mu} (length-\code{K}
#'   list of bin mean vectors), \code{sigma2} (length-\code{d} pooled
#'   within-bin variances), \code{keep_idx} (rows kept after \code{NA}
#'   handling), and \code{cluster_rank} (bin index of each kept row, in
#'   increasing score order).
#' @noRd
.cavi_bin_ordering <- function(
    X, ordering_vec, K,
    discretization = c("quantile", "equal", "kmeans"),
    na_action = c("drop", "error"),
    eps = 1e-12
) {
  X <- as.matrix(X)
  discretization <- match.arg(discretization)
  na_action <- match.arg(na_action)

  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be an integer >= 2.")
  if (length(ordering_vec) != nrow(X)) stop("ordering_vec must have length nrow(X).")

  keep_idx <- seq_len(nrow(X))
  if (anyNA(ordering_vec)) {
    if (na_action == "error") stop("ordering_vec contains NA.")
    keep_idx <- which(!is.na(ordering_vec))
    ordering_vec <- ordering_vec[keep_idx]
    X <- X[keep_idx, , drop = FALSE]
  }

  n <- nrow(X)
  d <- ncol(X)
  if (n < K) stop("After NA handling, n must be >= K.")

  # ---- discretize ordering_vec into K bins / clusters ----
  if (discretization == "kmeans") {
    cl <- stats::kmeans(ordering_vec, centers = K)
    centers <- as.numeric(cl$centers)
    cluster_ids <- cl$cluster

  } else if (discretization == "quantile") {
    probs <- seq(0, 1, length.out = K + 1)
    brks <- as.numeric(stats::quantile(ordering_vec, probs = probs, na.rm = TRUE))
    brks <- unique(brks)
    if (length(brks) < 2) stop("ordering_vec has too few unique values for quantile discretization.")

    cluster_ids <- cut(ordering_vec, breaks = brks, include.lowest = TRUE, labels = FALSE)

    if (all(is.na(cluster_ids))) stop("Quantile discretization failed (all clusters NA).")
    if (anyNA(cluster_ids)) {
      warning("Some values were not assigned to any quantile bin; dropping them.")
      ok <- which(!is.na(cluster_ids))
      cluster_ids <- cluster_ids[ok]
      ordering_vec <- ordering_vec[ok]
      X <- X[ok, , drop = FALSE]
      keep_idx <- keep_idx[ok]
      n <- nrow(X)
    }

    K_eff <- max(cluster_ids, na.rm = TRUE)
    if (K_eff != K) {
      warning("Quantile discretization produced ", K_eff, " groups (requested K=", K, ").")
      K <- K_eff
    }

    centers <- tapply(ordering_vec, cluster_ids, mean)
    centers <- as.numeric(centers)

  } else { # "equal"
    brks <- seq(min(ordering_vec), max(ordering_vec), length.out = K + 1)
    cluster_ids <- cut(ordering_vec, breaks = brks, include.lowest = TRUE, labels = FALSE)
    if (all(is.na(cluster_ids))) stop("Equal discretization failed (all clusters NA).")

    centers <- tapply(ordering_vec, cluster_ids, mean)
    centers <- as.numeric(centers)
  }

  # Assign cluster_rank: rank in ascending-t order (1 = smallest t, K = largest t).
  #
  # For "kmeans" and "quantile", all K bins are non-empty, so
  #   cluster_order = order(centers)  is a permutation of 1..K, and
  #   match(cluster_ids, cluster_order) correctly maps each cluster ID to its rank.
  #
  # For "equal", cut() may produce empty bins, so centers has length < K and
  #   cluster_order = 1..length(centers) != 1..K.  Calling match() against this
  #   shorter sequence returns NA for any bin ID > length(centers), which forces
  #   those points to pseudotime = 0 - a silent, non-monotone bug.
  #   Since equal-width bins from cut() are already ordered by t (bin 1 = lowest,
  #   bin K = highest), the bin ID itself is the correct rank.
  if (discretization == "equal") {
    cluster_rank <- cluster_ids          # already 1..K in ascending-t order
  } else {
    cluster_order <- order(centers)
    cluster_rank  <- match(cluster_ids, cluster_order)
  }

  # ---- mu ----
  # First pass: compute means for non-empty bins (NULL for empty ones).
  mu_raw <- lapply(seq_len(K), function(k) {
    idx <- which(cluster_rank == k)
    if (length(idx) == 0L) NULL else colMeans(X[idx, , drop = FALSE])
  })

  # Second pass: fill empty bins by linear interpolation between the nearest
  # non-empty neighbours on each side.  Edge cases (no left / no right
  # neighbour) fall back to the single available neighbour (constant extrapolation).
  non_empty <- which(!vapply(mu_raw, is.null, logical(1L)))
  if (length(non_empty) == 0L) stop("All K position bins are empty.")

  mu_list <- vector("list", K)
  for (k in seq_len(K)) {
    if (!is.null(mu_raw[[k]])) {
      mu_list[[k]] <- mu_raw[[k]]
    } else {
      left_nb  <- non_empty[non_empty < k]
      right_nb <- non_empty[non_empty > k]
      has_left  <- length(left_nb)  > 0L
      has_right <- length(right_nb) > 0L

      if (has_left && has_right) {
        l <- max(left_nb);  r <- min(right_nb)
        w_r <- (k - l) / (r - l)          # weight toward right neighbour
        mu_list[[k]] <- (1 - w_r) * mu_raw[[l]] + w_r * mu_raw[[r]]
      } else if (has_left) {
        mu_list[[k]] <- mu_raw[[max(left_nb)]]
      } else {
        mu_list[[k]] <- mu_raw[[min(right_nb)]]
      }
    }
  }

  # ---- pi ----
  Nk <- tabulate(cluster_rank, nbins = K)
  pi_vec <- Nk / n

  # ---- sigma2: pooled within-bin residuals, one variance per feature ----
  diffs <- do.call(rbind, lapply(seq_len(K), function(k) {
    idx <- which(cluster_rank == k)
    if (length(idx) == 0) return(NULL)
    sweep(X[idx, , drop = FALSE], 2, mu_list[[k]])
  }))
  sigma2 <- if (is.null(diffs) || nrow(diffs) == 0) rep(1, d) else colMeans(diffs^2)
  sigma2 <- pmax(as.numeric(sigma2), eps)

  list(
    pi = pi_vec,
    mu = mu_list,
    sigma2 = sigma2,
    keep_idx = keep_idx,
    cluster_rank = cluster_rank
  )
}

#' Initialize one ordering for a CAVI fit
#'
#' @description
#' Computes an ordering score with \code{method} and cuts it into \code{K}
#' position bins with \code{.cavi_bin_ordering()}. \code{"random"} uses a
#' random permutation of the samples and follows the same path. A failed
#' non-PCA ordering falls back to PCA component 1 when
#' \code{init_on_failure = "pca"}.
#'
#' @param X Numeric matrix (\code{n x d}) without missing values.
#' @param K Integer >= 2; requested number of bins.
#' @param method Ordering method.
#' @param discretization Passed to \code{.cavi_bin_ordering()}.
#' @param init_on_failure \code{"pca"} or \code{"error"}.
#' @param method_args Named list of arguments for the ordering helper.
#'
#' @return A list with \code{params} (bin proportions \code{pi}, bin means
#'   \code{mu}, pooled variances \code{sigma2}), \code{keep_idx},
#'   \code{cluster_rank}, and \code{ordering} (the ordering result, whose
#'   \code{init_info} attribute records the method used and any fallback).
#' @noRd
.cavi_initialize_ordering <- function(
    X, K,
    method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
    discretization = c("quantile", "equal", "kmeans"),
    init_on_failure = "pca",
    method_args = list()
) {
  method <- match.arg(method)
  discretization <- match.arg(discretization)
  X <- as.matrix(X)
  if (nrow(X) < 1L) stop("X must have at least one row.")
  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be a single integer >= 2.")

  ordering_result <- .cavi_get_ordering_result(
    X, method, method_args = method_args, init_on_failure = init_on_failure)
  init <- .cavi_bin_ordering(
    X = X,
    ordering_vec = ordering_result$t,
    K = K,
    discretization = discretization,
    na_action = "drop"
  )

  ordering <- c(ordering_result,
    list(method = attr(ordering_result, "init_info")$method_used))
  attr(ordering, "init_info") <- attr(ordering_result, "init_info")
  list(
    params       = list(pi = init$pi, mu = init$mu, sigma2 = init$sigma2),
    keep_idx     = init$keep_idx,
    cluster_rank = init$cluster_rank,
    ordering     = ordering
  )
}
