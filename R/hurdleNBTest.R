#' @importFrom stats optim pnorm pchisq plogis dnbinom dnorm qlogis
NULL
#' Negative Log-Likelihood for Hurdle NB with Fixed Zero Probabilities and Dispersion
#'
#' Internal function to compute the negative log-likelihood of a hurdle
#' negative binomial model where both zero probabilities (\code{theta}) and
#' dispersion (\code{phi}) are fixed. Only mean parameters are optimized.
#'
#' @param params Numeric vector:
#'   \describe{
#'     \item{1}{If \code{one_group = TRUE}, log(mu) shared across groups.}
#'     \item{1}{If \code{one_group = FALSE}, log(mu_A) for group A.}
#'     \item{2}{If \code{one_group = FALSE}, delta_rate = log(mu_B) - log(mu_A).}
#'   }
#' @param counts Numeric vector of observed counts for a single gene.
#' @param groups Factor indicating group membership for each observation.
#' @param offset Optional numeric vector of normalization factors. If \code{NULL},
#'   defaults to 1 for all observations.
#' @param fixed_theta_A Numeric, fixed non-zero probability for group A.
#' @param fixed_theta_B Numeric, fixed non-zero probability for group B.
#' @param fixed_phi Numeric, fixed NB dispersion parameter (size).
#' @param one_group Logical; if \code{TRUE}, a single mean is shared across both groups.
#'
#' @return Numeric, the negative log-likelihood value.
#'
#' @details
#' The likelihood consists of:
#' \itemize{
#'   \item Bernoulli component for zero vs non-zero counts using fixed theta
#'   \item Truncated negative binomial component for positive counts using fixed phi
#' }
#'
#' This function is used internally during optimization when only mean parameters
#' are estimated, while theta and phi are kept constant.
#'
#' @keywords internal
#' @noRd
nll_hurdle_fixed_theta_phi <- function(params, counts, groups, offset = NULL,
                                       fixed_theta_A, fixed_theta_B,
                                       fixed_phi, one_group = TRUE) {
  eps <- 1e-12

  # --- 0. Setup ---
  if (is.null(offset)) offset <- rep(1, length(counts))
  levs <- levels(groups)

  # --- 1. Map Parameters ---
  # Only log_rates are in params now
  if (one_group) {
    # params = [log_rate_shared]
    log_rate_A <- params[1]
    log_rate_B <- params[1]
  } else {
    # params = [log_rate_A, delta_rate]
    log_rate_A <- params[1]
    delta_rate <- params[2]
    log_rate_B <- log_rate_A + delta_rate
  }

  # Map fixed thetas and rates to observations
  is_grp_A <- (groups == levs[1])
  thetas <- ifelse(is_grp_A, fixed_theta_A, fixed_theta_B)
  thetas <- pmin(pmax(thetas, eps), 1 - eps)
  rates  <- ifelse(is_grp_A, exp(log_rate_A), exp(log_rate_B))
  mus    <- rates * offset

  # --- 2. Likelihood Calculation ---
  is_zero <- counts == 0

  # A. Binary Part (using fixed thetas)
  ll_binary <- sum(log1p(-thetas[is_zero])) + sum(log(thetas[!is_zero]))

  # B. Truncated Negative Binomial Part
  pos_counts <- counts[!is_zero]
  pos_mus    <- mus[!is_zero]

  # Standard NB log-likelihood
  ll_nb <- sum(dnbinom(pos_counts, mu = pos_mus, size = fixed_phi, log = TRUE))

  # Truncation Adjustment: log(1 / (1 - P(Y=0)))
  prob_zero_nb     <- dnbinom(0, mu = pos_mus, size = fixed_phi)
  prob_zero_nb     <- pmin(pmax(prob_zero_nb, 0), 1 - eps)
  trunc_adjustment <- -sum(log1p(-prob_zero_nb))

  return(-(ll_binary + ll_nb + trunc_adjustment))
}

#' Negative Log-Likelihood for Hurdle NB with Fixed Dispersion
#'
#' Internal function to compute the negative log-likelihood of a hurdle
#' negative binomial model with fixed dispersion (\code{phi}). The nonzero
#' probabilities (\code{theta}) and mean parameters are optimized.
#'
#' @param params Numeric vector:
#'   \describe{
#'     \item{1}{If \code{one_group = TRUE}, shared logit nonzero probability.}
#'     \item{2}{If \code{one_group = TRUE}, shared log mean parameter.}
#'     \item{1--2}{If \code{one_group = FALSE}, group-specific logit nonzero
#'       probabilities for groups A and B.}
#'     \item{3}{If \code{one_group = FALSE}, log mean parameter for group A.}
#'     \item{4}{If \code{one_group = FALSE},
#'       \code{delta_rate = log(mu_B) - log(mu_A)}.}
#'   }
#' @param counts Numeric vector of observed counts for a single gene.
#' @param groups Factor indicating group membership for each observation.
#' @param offset Optional numeric vector of normalization factors. If
#'   \code{NULL}, defaults to 1 for all observations.
#' @param fixed_phi Numeric, fixed negative binomial dispersion parameter
#'   (\code{size}).
#' @param one_group Logical; if \code{TRUE}, the nonzero probability and mean
#'   parameters are shared across both groups.
#'
#' @return Numeric, the negative log-likelihood value.
#'
#' @details
#' The likelihood consists of:
#' \itemize{
#'   \item A Bernoulli component for zero versus nonzero counts.
#'   \item A zero-truncated negative binomial component for positive counts
#'         using the fixed dispersion parameter.
#' }
#'
#' This function is used internally during optimization when the nonzero
#' probabilities and mean parameters are estimated while dispersion is held
#' fixed.
#'
#' @keywords internal
#' @noRd
nll_hurdle_fixed_phi <- function(params, counts, groups, offset = NULL,
                                 fixed_phi, one_group = TRUE) {
  eps <- 1e-12

  # --- 0. Setup ---
  if (is.null(offset)) offset <- rep(1, length(counts))
  levs <- levels(groups)

  # --- 1. Map Parameters ---
  if (one_group) {
    # params = [logit_theta_shared, log_rate_shared]
    logit_theta_A <- params[1]
    logit_theta_B <- params[1]

    log_rate_A <- params[2]
    log_rate_B <- params[2]
  } else {
    # params = [logit_theta_A, logit_theta_B, log_rate_A, delta_rate]
    logit_theta_A <- params[1]
    logit_theta_B <- params[2]

    log_rate_A <- params[3]
    delta_rate <- params[4]
    log_rate_B <- log_rate_A + delta_rate

  }

  theta_A <- plogis(logit_theta_A)
  theta_B <- plogis(logit_theta_B)

  # Map thetas and rates to observations
  is_grp_A <- (groups == levs[1])
  thetas <- ifelse(is_grp_A, theta_A, theta_B)
  thetas <- pmin(pmax(thetas, eps), 1 - eps)
  rates  <- ifelse(is_grp_A, exp(log_rate_A), exp(log_rate_B))
  mus    <- rates * offset

  # --- 2. Likelihood Calculation ---
  is_zero <- counts == 0

  # A. Binary Part
  ll_binary <- sum(log1p(-thetas[is_zero])) + sum(log(thetas[!is_zero]))

  # B. Truncated Negative Binomial Part
  pos_counts <- counts[!is_zero]
  pos_mus    <- mus[!is_zero]

  # Standard NB log-likelihood
  ll_nb <- sum(dnbinom(pos_counts, mu = pos_mus, size = fixed_phi, log = TRUE))

  # Truncation Adjustment: log(1 / (1 - P(Y=0)))
  prob_zero_nb <- dnbinom(0, mu = pos_mus, size = fixed_phi)
  prob_zero_nb <- pmin(pmax(prob_zero_nb, 0), 1 - eps)
  trunc_adjustment <- -sum(log1p(-prob_zero_nb))

  return(-(ll_binary + ll_nb + trunc_adjustment))
}

#' Mean-Based Likelihood Ratio Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise likelihood ratio tests (LRTs) for differential expression
#' in the expression mean using a hurdle negative binomial model with fixed
#' zero probabilities and tag-wise dispersions.
#'
#' @param y A list-like object returned from \code{tagwiseEst()} containing:
#'   \describe{
#'     \item{counts}{Numeric matrix of gene expression counts (genes x samples).}
#'     \item{samples}{Data frame with columns \code{group} (factor) and \code{size.factor} (numeric).}
#'     \item{tagwise.disp}{Numeric vector of estimated tag-wise dispersions.}
#'     \item{zero_prob_matrix}{Matrix of zero probabilities per group per gene.}
#'   }
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{lrt_stats}{Numeric vector of LRT statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' For each gene:
#' \itemize{
#'   \item The null model assumes a single shared mean across groups.
#'   \item The alternative model estimates group-specific means.
#'   \item Zero probabilities and dispersions are fixed at estimates obtained
#'         from \code{\link{tagwiseEst}}.
#'   \item When one group has all zero counts, a one-sided Z-test is applied instead.
#' }
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' y <- hurdle_LRT.mean(y)
#'
#' @export
hurdle_LRT.mean <- function(y){
  n_row <- nrow(y$counts)
  lrt_stats <- numeric(n_row)
  p_values  <- numeric(n_row)

  groups_factor <- factor(y$samples$group)
  offsets       <- y$samples$size.factor
  levs          <- levels(groups_factor)

  for (i in 1:n_row) {
    current_counts <- as.numeric(y$counts[i, ])
    current_phi <- 1 / y$tagwise.disp[i]

    # Retrieve the fixed priors for this gene
    f_theta_A <- 1 - y$zero_prob_matrix[i, "prob_A"]
    f_theta_B <- 1 - y$zero_prob_matrix[i, "prob_B"]

    all_zero_A <- all(current_counts[groups_factor == levs[1]] == 0)
    all_zero_B <- all(current_counts[groups_factor == levs[2]] == 0)

    # Initial estimate for log_rate (shared)
    obs_mean <- mean(current_counts[current_counts > 0] / offsets[current_counts > 0], na.rm=TRUE)
    if(is.na(obs_mean) || obs_mean <= 0) obs_mean <- 0.1

    tryCatch({
      # --- Model 1: Null (1 Shared Mean) ---
      # We estimate 1 parameter: log_rate
      null_model <- optim(
        par = log(obs_mean),
        fn = nll_hurdle_fixed_theta_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_theta_A = f_theta_A,
        fixed_theta_B = f_theta_B,
        fixed_phi = current_phi,
        one_group = TRUE,
        method = "BFGS",
        hessian = TRUE
      )

      if (all_zero_A || all_zero_B) {
        # when one group has all zeros, test the mean of the other group
        se_beta <- sqrt(1 / null_model$hessian[1, 1])

        # Z-test: (Observed log-rate - "Zero" Floor) / SE
        # testing if the shared mean is significantly above a floor
        z_stat <- (null_model$par - 1) / se_beta

        lrt_stats[i] <- z_stat
        # One-sided test (upper tail)
        p_values[i]  <- pnorm(z_stat, lower.tail = FALSE)

      } else {
        # --- Model 2: Alt (Standard LRT) ---
        m_A <- mean(current_counts[groups_factor == levs[1] & current_counts > 0] /
                      offsets[groups_factor == levs[1] & current_counts > 0], na.rm=TRUE)
        m_B <- mean(current_counts[groups_factor == levs[2] & current_counts > 0] /
                      offsets[groups_factor == levs[2] & current_counts > 0], na.rm=TRUE)

        alt_model <- optim(
          par = c(log(m_A), log(m_B) - log(m_A)),
          fn = nll_hurdle_fixed_theta_phi,
          counts = current_counts,
          groups = groups_factor,
          offset = offsets,
          fixed_theta_A = f_theta_A,
          fixed_theta_B = f_theta_B,
          fixed_phi = current_phi,
          one_group = FALSE,
          method = "BFGS"
        )

        lrt <- 2 * (null_model$value - alt_model$value)
        lrt_stats[i] <- max(0, lrt)
        p_values[i]  <- pchisq(lrt_stats[i], df = 1, lower.tail = FALSE)
      }

    }, error = function(e) {
      lrt_stats[i] <<- NA
      p_values[i]  <<- NA
    })
  }

  y$lrt_stats <- lrt_stats
  y$p.values  <- p_values
  return(y)
}

#' Distributional Likelihood Ratio Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise likelihood ratio tests (LRTs) for distributional
#' differential expression using a hurdle negative binomial model. The test
#' jointly evaluates differences in the zero probability and positive-count
#' mean between groups while holding the tag-wise dispersion fixed.
#'
#' @param y A list-like object returned from \code{\link{tagwiseEst}} containing:
#'   \describe{
#'     \item{counts}{Numeric matrix of gene expression counts (genes x samples).}
#'     \item{samples}{Data frame with columns \code{group} and
#'       \code{size.factor}.}
#'     \item{tagwise.disp}{Numeric vector of estimated tag-wise dispersions.}
#'   }
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{lrt_stats}{Numeric vector of LRT statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' For each gene:
#' \itemize{
#'   \item The null model assumes a shared nonzero probability and a shared
#'         positive-count mean across groups.
#'   \item The alternative model estimates group-specific nonzero probabilities
#'         and positive-count means.
#'   \item The dispersion is fixed at its tag-wise estimate obtained from
#'         \code{\link{tagwiseEst}}.
#'   \item The LRT uses two degrees of freedom when both groups contain positive
#'         counts and one degree of freedom when either group contains only
#'         zero counts.
#' }
#'
#' Thus, rejection of the null hypothesis indicates a difference in the overall
#' expression distribution arising from the zero probability, positive-count
#' mean, or both.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' y <- hurdle_LRT.dist(y)
#'
#' @export
hurdle_LRT.dist <- function(y) {
  n_row <- nrow(y$counts)
  lrt_stats <- numeric(n_row)
  p_values  <- numeric(n_row)

  groups_factor <- factor(y$samples$group)
  offsets       <- y$samples$size.factor
  levs          <- levels(groups_factor)

  for (i in 1:n_row) {
    current_counts <- as.numeric(y$counts[i, ])
    current_phi <- 1 / y$tagwise.disp[i]

    all_zero_A <- all(current_counts[groups_factor == levs[1]] == 0)
    all_zero_B <- all(current_counts[groups_factor == levs[2]] == 0)

    # Initial shared positive probability
    p_shared <- mean(current_counts > 0)
    p_shared <- max(0.01, min(0.99, p_shared))

    # Initial shared positive-count mean
    obs_mean <- mean(
      current_counts[current_counts > 0] /
        offsets[current_counts > 0],
      na.rm = TRUE
    )

    if (is.na(obs_mean) || obs_mean <= 0) obs_mean <- 0.1

    tryCatch({
      # --- Model 1: Null ---
      # Shared theta and shared count mean
      null_model <- optim(
        par = c(
          qlogis(p_shared),
          log(obs_mean)
        ),
        fn = nll_hurdle_fixed_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_phi = current_phi,
        one_group = TRUE,
        method = "BFGS"
      )

      df <- 2
      if (all_zero_A || all_zero_B){
        df <- 1
      }

      p_A <- mean(current_counts[groups_factor == levs[1]] > 0)
      p_B <- mean(current_counts[groups_factor == levs[2]] > 0)

      p_A <- max(0.01, min(0.99, p_A))
      p_B <- max(0.01, min(0.99, p_B))

      # Initial group-specific positive-count means
      m_A <- mean(current_counts[groups_factor == levs[1] & current_counts > 0] /
                    offsets[groups_factor == levs[1] & current_counts > 0], na.rm = TRUE)

      m_B <- mean(current_counts[groups_factor == levs[2] & current_counts > 0] /
                    offsets[groups_factor == levs[2] & current_counts > 0], na.rm = TRUE)

      if (!is.finite(m_A) || m_A <= 0) m_A <- obs_mean
      if (!is.finite(m_B) || m_B <= 0) m_B <- obs_mean

      # --- Model 2: Alternative ---
      # Group-specific theta and count mean
      alt_model <- optim(
        par = c(
          qlogis(p_A),
          qlogis(p_B),
          log(m_A),
          log(m_B) - log(m_A)
        ),
        fn = nll_hurdle_fixed_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_phi = current_phi,
        one_group = FALSE,
        method = "BFGS"
      )

      lrt <- 2 * (null_model$value - alt_model$value)
      if (!is.finite(lrt) || lrt < -1e-6) {
        stop("Invalid likelihood-ratio statistic.")
      }
      lrt_stats[i] <- max(0, lrt)
      p_values[i] <- pchisq(lrt_stats[i],df = df,lower.tail = FALSE)

    }, error = function(e) {
      lrt_stats[i] <<- NA
      p_values[i]  <<- NA
    })
  }

  y$lrt_stats <- lrt_stats
  y$p.values  <- p_values

  return(y)
}

#' Likelihood Ratio Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise differential expression testing using a hurdle negative
#' binomial model. The function can test differences in either the expression
#' mean or the overall expression distribution between groups.
#'
#' @param y A list-like object returned from \code{\link{tagwiseEst}} containing
#'   the count matrix, sample information, and estimated tag-wise dispersions.
#' @param test Character string specifying the hypothesis to test. Use
#'   \code{"mean"} to test differences in the positive-count mean while fixing
#'   zero probabilities and dispersions, or \code{"distribution"} to jointly
#'   test differences in the zero probability and positive-count mean.
#'   The default is \code{"mean"}.
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{lrt_stats}{Numeric vector of test statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' When \code{test = "mean"}, the function calls
#' \code{\link{hurdle_LRT.mean}}. When \code{test = "distribution"}, it calls
#' \code{\link{hurdle_LRT.dist}}.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#'
#' mean_test <- hurdle.LRT(y, test = "mean")
#' distribution_test <- hurdle.LRT(y, test = "distribution")
#'
#' @export
hurdle.LRT <- function(y, test = c("mean", "distribution")) {
  test <- match.arg(test)

  if (test == "mean") {
    return(hurdle_LRT.mean(y))
  }
  return(hurdle_LRT.dist(y))
}


#' Mean-Based Wald Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise Wald tests for differential expression in the expression
#' mean using a hurdle negative binomial model with fixed zero probabilities
#' and tag-wise dispersions.
#'
#' @param y A list-like object returned from \code{\link{tagwiseEst}} containing:
#'   \describe{
#'     \item{counts}{Numeric matrix of gene expression counts (genes x samples).}
#'     \item{samples}{Data frame with columns \code{group} and
#'       \code{size.factor}.}
#'     \item{tagwise.disp}{Numeric vector of estimated tag-wise dispersions.}
#'     \item{zero_prob_matrix}{Matrix of estimated zero probabilities for each
#'       gene and group.}
#'   }
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{wald_stats}{Numeric vector of Wald or one-sided Z statistics for
#'       each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' For each gene:
#' \itemize{
#'   \item The model estimates group-specific mean parameters while fixing the
#'         zero probabilities and dispersion at estimates obtained from
#'         \code{\link{tagwiseEst}}.
#'   \item When both groups contain positive counts, a two-sided Wald test is
#'         applied to the log mean difference between groups.
#'   \item When either group contains only zero counts, a one-sided Z-test is
#'         applied to the estimated log mean of the nonzero group using the
#'         corresponding Hessian-based standard error.
#' }
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' y <- hurdle_Wald_Test.mean(y)
#'
#' @export
hurdle_Wald_Test.mean <- function(y) {
  n_row <- nrow(y$counts)

  # Initialize output vectors
  wald_stats <- rep(NA, n_row)
  p_values   <- rep(1, n_row)

  groups_factor <- factor(y$samples$group)
  offsets       <- y$samples$size.factor
  levs          <- levels(groups_factor)

  for (i in 1:n_row) {
    current_counts <- as.numeric(y$counts[i, ])
    current_phi <- 1 / y$tagwise.disp[i]
    f_theta_A <- 1 - y$zero_prob_matrix[i, "prob_A"]
    f_theta_B <- 1 - y$zero_prob_matrix[i, "prob_B"]

    # Check for all-zero groups
    is_zero_A <- all(current_counts[groups_factor == levs[1]] == 0)
    is_zero_B <- all(current_counts[groups_factor == levs[2]] == 0)

    # Warm start
    m_A <- mean((current_counts[groups_factor == levs[1]] + 0.1) /
                  offsets[groups_factor == levs[1]], na.rm = TRUE)
    m_B <- mean((current_counts[groups_factor == levs[2]] + 0.1) /
                  offsets[groups_factor == levs[2]], na.rm = TRUE)

    tryCatch({
      fit <- optim(
        par = c(log(m_A), log(m_B) - log(m_A)),
        fn = nll_hurdle_fixed_theta_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_theta_A = f_theta_A,
        fixed_theta_B = f_theta_B,
        fixed_phi = current_phi,
        one_group = FALSE,
        method = "BFGS",
        hessian = TRUE
      )

      H <- fit$hessian

      # Logic for One Group All Zeros
      if (is_zero_A || is_zero_B) {

        # choose which curvature to use
        idx <- if (is_zero_A) 2 else 1   # A all zero -> use delta direction; B all zero -> use A direction

        if (is.na(H[idx, idx]) || H[idx, idx] <= 1e-10) {
          wald_stats[i] <- NA
          p_values[i]   <- NA
        } else {
          se_beta <- sqrt(1 / H[idx, idx])

          # point estimate on log scale for the *non-zero group mean*
          beta_hat <- if (is_zero_A) (fit$par[1] + fit$par[2]) else fit$par[1]

          z_stat <- (beta_hat - 1) / se_beta
          wald_stats[i] <- z_stat
          p_values[i]   <- pnorm(z_stat, lower.tail = FALSE)  # one-sided
        }

      } else {
        # Logic for Standard Comparison (Neither is all zero)
        if (det(H) <= 1e-10 || any(is.na(H))) {
          wald_stats[i] <- NA
          p_values[i]   <- NA
        } else {
          vcov_mat <- solve(H)
          delta_beta <- fit$par[2]
          se_delta  <- sqrt(max(0, vcov_mat[2,2]))

          if (se_delta > 0) {
            z_stat <- delta_beta / se_delta
            wald_stats[i] <- z_stat
            # Two-sided for standard comparison
            p_values[i]   <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
          } else {
            p_values[i] <- 1
          }
        }
      }

    }, error = function(e) {
      wald_stats[i] <<- NA
      p_values[i]   <<- NA
    })
  }

  y$wald_stats <- wald_stats
  y$p.values   <- p_values
  return(y)
}

#' Distributional Wald Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise Wald tests for distributional differential expression
#' using a hurdle negative binomial model. The test jointly evaluates
#' differences in the zero probability and positive-count mean between groups
#' while holding the tag-wise dispersion fixed.
#'
#' @param y A list-like object returned from \code{\link{tagwiseEst}} containing:
#'   \describe{
#'     \item{counts}{Numeric matrix of gene expression counts (genes x samples).}
#'     \item{samples}{Data frame with columns \code{group} and
#'       \code{size.factor}.}
#'     \item{tagwise.disp}{Numeric vector of estimated tag-wise dispersions.}
#'   }
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{wald_stats}{Numeric vector of Wald statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' For each gene:
#' \itemize{
#'   \item When both groups contain positive counts, the model estimates
#'         group-specific nonzero probabilities and mean parameters while
#'         fixing the dispersion at its tag-wise estimate.
#'   \item A joint two-degree-of-freedom Wald test evaluates the nonzero
#'         probability and positive-count mean differences between groups.
#'   \item When either group contains only zero counts, the test compares only
#'         the nonzero probabilities using corrected log odds, with 0.5 added
#'         to each positive- and zero-count frequency.
#'   \item The resulting squared Z statistic is compared with a chi-squared
#'         distribution with one degree of freedom in the all-zero-group case.
#' }
#'
#' Rejection of the null hypothesis indicates a difference in the overall
#' expression distribution arising from the zero probability, positive-count
#' mean, or both.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' y <- hurdle_Wald_Test.dist(y)
#'
#' @export
hurdle_Wald_Test.dist <- function(y) {
  n_row <- nrow(y$counts)

  # Initialize output vectors
  wald_stats <- rep(NA, n_row)
  p_values   <- rep(1, n_row)

  groups_factor <- factor(y$samples$group)
  offsets       <- y$samples$size.factor
  levs          <- levels(groups_factor)

  for (i in 1:n_row) {
    current_counts <- as.numeric(y$counts[i, ])
    current_phi <- 1 / y$tagwise.disp[i]

    # Check for all-zero groups
    is_zero_A <- all(current_counts[groups_factor == levs[1]] == 0)
    is_zero_B <- all(current_counts[groups_factor == levs[2]] == 0)

    if (is_zero_A || is_zero_B) {
      counts_A <- current_counts[groups_factor == levs[1]]
      counts_B <- current_counts[groups_factor == levs[2]]

      positive_A <- sum(counts_A > 0) + 0.5
      zero_A     <- sum(counts_A == 0) + 0.5
      positive_B <- sum(counts_B > 0) + 0.5
      zero_B     <- sum(counts_B == 0) + 0.5

      logit_A <- log(positive_A / zero_A)
      logit_B <- log(positive_B / zero_B)

      se_difference <- sqrt(
        1 / positive_A +
          1 / zero_A +
          1 / positive_B +
          1 / zero_B
      )

      z_stat <- (logit_B - logit_A) / se_difference

      wald_stats[i] <- z_stat^2
      p_values[i] <- pchisq(
        wald_stats[i],
        df = 1,
        lower.tail = FALSE
      )

      next
    }

    # Initial zero-component parameters
    p_A <- mean(current_counts[groups_factor == levs[1]] > 0)
    p_B <- mean(current_counts[groups_factor == levs[2]] > 0)

    p_A <- max(0.01, min(0.99, p_A))
    p_B <- max(0.01, min(0.99, p_B))

    # Warm start
    m_A <- mean((current_counts[groups_factor == levs[1]] + 0.1) /
                  offsets[groups_factor == levs[1]], na.rm = TRUE)

    m_B <- mean((current_counts[groups_factor == levs[2]] + 0.1) /
                  offsets[groups_factor == levs[2]], na.rm = TRUE)

    tryCatch({
      fit <- optim(
        par = c(
          qlogis(p_A),
          qlogis(p_B),
          log(m_A),
          log(m_B) - log(m_A)
        ),
        fn = nll_hurdle_fixed_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_phi = current_phi,
        one_group = FALSE,
        method = "BFGS",
        hessian = TRUE
      )

      H <- fit$hessian

      # Standard joint distribution comparison
      if (det(H) <= 1e-10 || any(is.na(H))) {
        wald_stats[i] <- NA
        p_values[i]   <- NA
      } else {
        vcov_mat <- solve(H)

        contrast <- matrix(
          c(
            -1, 1, 0, 0,
            0, 0, 0, 1
          ),
          nrow = 2,
          byrow = TRUE
        )

        difference <- contrast %*% fit$par
        difference_vcov <- contrast %*% vcov_mat %*% t(contrast)

        wald_stat <- drop(
          t(difference) %*%
            solve(difference_vcov) %*%
            difference
        )

        wald_stats[i] <- max(0, wald_stat)
        p_values[i] <- pchisq(
          wald_stats[i],
          df = 2,
          lower.tail = FALSE
        )
      }

    }, error = function(e) {
      wald_stats[i] <<- NA
      p_values[i]   <<- NA
    })
  }

  y$wald_stats <- wald_stats
  y$p.values   <- p_values

  return(y)
}

#' Wald Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise differential expression testing using a hurdle negative
#' binomial model. The function can test differences in either the expression
#' mean or the overall expression distribution between groups.
#'
#' @param y A list-like object returned from \code{\link{tagwiseEst}} containing
#'   the count matrix, sample information, and estimated tag-wise dispersions.
#' @param test Character string specifying the hypothesis to test. Use
#'   \code{"mean"} to test differences in the positive-count mean while fixing
#'   zero probabilities and dispersions, or \code{"distribution"} to jointly
#'   test differences in the zero probability and positive-count mean.
#'   The default is \code{"mean"}.
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{wald_stats}{Numeric vector of Wald statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' When \code{test = "mean"}, the function calls
#' \code{\link{hurdle_Wald_Test.mean}}. When
#' \code{test = "distribution"}, it calls
#' \code{\link{hurdle_Wald_Test.dist}}.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#'
#' mean_test <- hurdle.Wald.Test(y, test = "mean")
#' distribution_test <- hurdle.Wald.Test(y, test = "distribution")
#'
#' @export
hurdle.Wald.Test <- function(y, test = c("mean", "distribution")) {
  test <- match.arg(test)

  if (test == "mean") {
    return(hurdle_Wald_Test.mean(y))
  }
  return(hurdle_Wald_Test.dist(y))
}
