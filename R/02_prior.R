# =========================
# 1D random walk precision
# =========================

#' Prior precision for q-th order random walk on K time points (d-dimensional)
#'
#' Constructs Q = lambda * (D_q' D_q) kronecker I_d, with optional ridge.
#'
#' @param K integer number of time points.
#' @param d integer dimension of the process at each time.
#' @param q integer order of the random walk (q>=1).
#' @param lambda nonnegative scalar precision multiplier.
#' @param ridge nonnegative scalar ridge added to the 1D precision.
#' @param sparse logical; return a sparse Matrix object if TRUE.
#' @return (dK) x (dK) precision matrix.
#' @noRd
make_random_walk_precision <- function(K, d, q = 1, lambda = 1, ridge = 0, sparse = FALSE) {
  if (K < 1 || K != as.integer(K)) stop("K must be a positive integer.")
  if (d < 1 || d != as.integer(d)) stop("d must be a positive integer.")
  if (q < 1 || q != as.integer(q) || q >= K) stop("q must be an integer in {1, ..., K-1}.")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (ridge < 0) stop("ridge must be >= 0.")

  # Dense construction (fine for moderate K); for large K you may want a fully sparse builder later
  D_q <- diag(K)
  for (i in seq_len(q)) D_q <- diff(D_q)

  Q_scalar <- crossprod(D_q) + ridge * diag(K)
  Q_U <- lambda * kronecker(Q_scalar, diag(d))

  if (sparse) {
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Please install.packages('Matrix')")
    Q_U <- Matrix::Matrix(Q_U, sparse = TRUE)
  }

  Q_U
}
