#' MPCurver package
#'
#' CAVI-first trajectory and pseudotime inference with legacy smoothEM
#' backends retained for internal benchmarking and compatibility checks.
#' Single-ordering fits use ordinary CAVI; fixed multi-ordering fits use the
#' structural family \eqn{q(C)\prod_j q(Z_j)q(U_j\mid Z_j)} with one canonical
#' shared state.
#'
#' @name MPCurver
#' @keywords internal
#' @importFrom graphics abline arrows axis barplot box legend lines mtext par plot.new plot.window points rect text title
#' @importFrom stats as.dist cor optimize prcomp setNames var
#' @importFrom utils tail
"_PACKAGE"

#' Unified MPCurve fit object
#'
#' The public return type of \code{\link{fit_mpcurve}} and
#' \code{\link{do_mpcurve}}. An \code{mpcurve} object normalizes the public
#' fields of the recommended CAVI backend while retaining the raw fit in
#' \code{$fit}.
#'
#' @section Single-ordering objects:
#' For \code{intrinsic_dim = 1}, \code{$params$pi} is a length-\code{K}
#' position-prior vector, \code{$params$mu} is a \code{d x K} posterior-mean
#' trajectory matrix, \code{$params$sigma2} is a length-\code{d} estimated
#' variance vector (or \code{NULL} when measurement standard deviations are
#' supplied through \code{S}), and \code{$gamma} is the \code{n x K}
#' responsibility matrix. The full Gaussian trajectory posterior is available
#' from \code{$fit$posterior}. When \code{S} is supplied,
#' \code{$measurement_sd} stores either its length-\code{d} feature-shared
#' form or its \code{n x d} observation-level form.
#'
#' @section Fixed-M structural partition objects:
#' For \code{intrinsic_dim = M >= 2}, the variational family is
#' \deqn{q(C)\prod_{j=1}^d q(Z_j)q(U_j\mid Z_j).}
#' The normalized object has one canonical structural state:
#' \itemize{
#'   \item \code{$params$pi}: named list of \code{M} length-\code{K}
#'     position-prior vectors;
#'   \item \code{$params$mu}: named list of \code{M} \code{d x K}
#'     conditional posterior-mean matrices;
#'   \item \code{$params$sigma2}: one shared length-\code{d} estimated
#'     variance vector, or \code{NULL} for known measurement noise;
#'   \item \code{$gamma}: named list of \code{M} \code{n x K}
#'     position-responsibility matrices;
#'   \item \code{$measurement_sd}: \code{NULL} when noise is estimated, or the
#'     supplied length-\code{d} / \code{n x d} known-standard-deviation
#'     specification;
#'   \item \code{$conditional_posterior}: named conditional Gaussian state for
#'     \eqn{q(U_j\mid Z_j=m)}. Each \code{$mean[[m]]} is \code{d x K};
#'     each \code{$cov[[m]][[j]]} is \code{K x K}; and
#'     \code{$var[[m]]} (also available as \code{$diag[[m]]}) is
#'     \code{d x K};
#'   \item \code{$lambda_mat}: a \code{d x M} matrix of smoothness
#'     parameters;
#'   \item \code{$partition$pi_weights}: a \code{d x M} matrix containing
#'     \eqn{q(Z_j=m)}, with the hard display assignment in
#'     \code{$partition$assign};
#'   \item \code{$objective_history} and \code{$temperature_history}: the
#'     fixed-\code{M} structural objective and its temperature schedule.
#'     Objective values at different temperatures are not directly
#'     comparable; use the exact \eqn{T=1} segment for monotonicity checks.
#' }
#' \code{$fits} contains named, derived single-ordering views for compatibility
#' with existing inspection and plotting code. Those views are not independent
#' fitting states and are never used to continue the structural algorithm.
#' Structural objects always retain their requested fixed \code{M}; the
#' compatibility fields \code{$active_intrinsic_dim} and
#' \code{$displayed_intrinsic_dim} therefore both equal \code{M}.
#'
#' @name mpcurve
#' @keywords models
NULL
