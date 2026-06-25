#' @importFrom stats optim pnorm pchisq plogis dnbinom dnorm qlogis
NULL
#' Negative Log-Likelihood for Hurdle Model
#'
#' Internal function to compute the negative log-likelihood for a hurdle
#' negative binomial model with optional prior on dispersion.
#'
#' @param params Numeric vector of length 5:
#'   \describe{
#'     \item{1}{logit(theta_A): 1 - zero probability for group A}
#'     \item{2}{logit(theta_B): 1 - zero probability for group B}
#'     \item{3}{log(mu_A): mean parameter for group A}
#'     \item{4}{delta_rate: log mean difference between groups}
#'     \item{5}{log(phi): dispersion parameter (log scale)}
#'   }
#' @param counts Numeric vector of counts.
#' @param groups Factor indicating group membership (2 levels).
#' @param offset Optional numeric vector of offsets (default: all ones).
#' @param priors Optional numeric vector of length 2 specifying mean and standard deviation
#'   for a normal prior on log(phi). Default is \code{NULL}, meaning no prior is applied.
#'
#' @return Numeric scalar: negative log-likelihood.
#'
#' @keywords internal
#' @noRd
nll_hurdle <- function(params, counts, groups, offset = NULL, priors = NULL){
  eps <- 1e-12

  # --- 0. input checks ---
  if (length(params) != 5) {
    stop("'params' must have length 5.")
  }

  groups <- as.factor(groups)
  if (length(groups) != length(counts)) {
    stop("'groups' must have the same length as 'counts'.")
  }
  if (nlevels(groups) != 2) {
    stop("'groups' must have exactly 2 levels.")
  }

  if (is.null(offset)) {
    offset <- rep(1, length(counts))
  }
  if (length(offset) != length(counts)) {
    stop("'offset' must have the same length as 'counts'.")
  }

  if (!is.null(priors)) {
    if (length(priors) != 2) {
      stop("'priors' must be NULL or a vector of length 2.")
    }
    if (priors[2] <= 0) {
      stop("prior sd must be positive.")
    }
  }

  # --- 1. Unpack Parameters ---
  logit_theta_A <- params[1]
  logit_theta_B <- params[2]
  log_rate_A    <- params[3]
  delta_rate    <- params[4]
  log_rate_B    <- log_rate_A + delta_rate
  phi           <- exp(params[5]) # phi is defined as 1/alpha where alpha is the overdispersion

  theta_A <- plogis(logit_theta_A)
  theta_B <- plogis(logit_theta_B)

  # --- 2. Assign Parameters to Data Points ---
  levs <- levels(groups)
  thetas <- ifelse(groups == levs[1], theta_A, theta_B)
  thetas <- pmin(pmax(thetas, eps), 1 - eps)

  base_rates <- ifelse(groups == levs[1], exp(log_rate_A), exp(log_rate_B))
  mus <- base_rates * offset

  # --- 3. Likelihood Calculation ---
  is_zero <- counts == 0

  # binary part
  ll_binary <- sum(log1p(-thetas[is_zero])) + sum(log(thetas[!is_zero]))

  # truncated NB part
  pos_counts <- counts[!is_zero]
  pos_mus <- mus[!is_zero]
  ll_nb <- sum(dnbinom(pos_counts, mu = pos_mus, size = phi, log = TRUE))
  prob_zero_nb <- dnbinom(0, mu = pos_mus, size = phi)
  prob_zero_nb <- pmin(pmax(prob_zero_nb, 0), 1 - eps)
  trunc_adjustment <- -sum(log1p(-prob_zero_nb))
  ll_count <- ll_nb + trunc_adjustment

  # prior part
  if (is.null(priors)) {
    prior_phi <- 0
  } else{
    prior_phi <- dnorm(params[5], mean = priors[1], sd = priors[2], log = TRUE)
  }

  total_ll <- ll_binary + ll_count + prior_phi
  return(-total_ll)
}


#' Estimate Bin-wise Priors for Hurdle Model Parameters (Internal)
#'
#' Internal function to group genes by similar mean expression levels and perform
#' bin-wise estimation of prior parameters for the hurdle negative binomial model.
#' Within each bin, genes are pooled and a global model is fitted to estimate
#' group-specific non-zero probabilities and log dispersion.
#'
#' @param y A list object produced by \code{sizeFactorsEst()}, containing
#' \code{counts}, \code{samples}, and \code{baseMean}.
#' @param n_bins Optional integer specifying the number of bins. If \code{NULL},
#' the number of bins is determined automatically.
#'
#' @details
#' Genes are partitioned into bins based on \code{log(baseMean)}. For each bin,
#' a hurdle negative binomial model is fitted using \code{nll_hurdle()} on pooled
#' counts to obtain bin-level estimates of:
#' \itemize{
#'   \item Group-specific non-zero probabilities
#'   \item Log-dispersion parameter (log(phi))
#' }
#'
#' The estimated parameters are mapped back to individual genes to provide
#' stabilized priors for downstream modeling.
#'
#' @return
#' The input object \code{y} with added components:
#' \describe{
#'   \item{prob_matrix}{Matrix of group-specific non-zero probabilities.}
#'   \item{prior_log_phi_gene}{Gene-wise log(phi) estimates.}
#'   \item{prior_bins}{Data.frame of bin-level parameter estimates.}
#' }
#'
#' @keywords internal
priorEst <- function(y, n_bins = NULL) {
  counts <- y$counts
  offset <- y$samples$size.factor
  grp    <- as.factor(y$samples$group)

  # Helper function for inverse logit
  inv_logit <- function(x) 1 / (1 + exp(-x))

  logmeans <- log(y$baseMean)
  finite   <- is.finite(logmeans)

  min_logmean <- min(logmeans[finite])
  max_logmean <- max(logmeans[finite])
  bin_width   <- 0.5

  if (min_logmean < 0.5) {
    first_upper <- 0.5
    if (is.null(n_bins)) {
      n_bins <- ceiling((max_logmean - first_upper) / bin_width) + 1
    }
  } else {
    first_upper <- min_logmean + bin_width
    if (is.null(n_bins)) {
      n_bins <- ceiling((max_logmean - min_logmean) / bin_width)
    }
  }

  # Containers for bin-level estimates
  local_log_phi <- rep(NA_real_, n_bins)
  local_prob_A  <- rep(NA_real_, n_bins)
  local_prob_B  <- rep(NA_real_, n_bins)
  bin_lower     <- rep(NA_real_, n_bins)
  bin_upper     <- rep(NA_real_, n_bins)

  for (j in seq_len(n_bins)) {
    if (j == 1 && min_logmean < 0.5) {
      lower <- min_logmean
      upper <- first_upper
    } else if (min_logmean < 0.5) {
      lower <- first_upper + (j - 2) * bin_width
      upper <- lower + bin_width
    } else {
      lower <- min_logmean + (j - 1) * bin_width
      upper <- lower + bin_width
    }

    if (j == n_bins) {
      upper <- max_logmean + 1e-8
    }

    bin_lower[j] <- lower
    bin_upper[j] <- upper

    idx <- finite & logmeans < upper & logmeans >= lower
    if (!any(idx)) next

    genes <- counts[idx, , drop = FALSE]
    n_row <- nrow(genes)

    df_long <- data.frame(
      group = rep(grp, times = n_row),
      count = as.vector(t(genes))
    )
    df_long$group <- droplevels(df_long$group)
    levs <- levels(df_long$group)

    # Initial params estimation
    pA_init <- max(0.01, min(0.99, mean(df_long$count[df_long$group == levs[1]] > 0)))
    pB_init <- max(0.01, min(0.99, mean(df_long$count[df_long$group == levs[2]] > 0)))

    mean_A <- mean(df_long$count[df_long$group == levs[1] & df_long$count > 0])
    mean_B <- mean(df_long$count[df_long$group == levs[2] & df_long$count > 0])

    if (is.nan(mean_A)) mean_A <- mean(df_long$count)
    if (is.nan(mean_B)) mean_B <- mean(df_long$count)

    log_rate_A_init <- log(mean_A)
    delta_rate_init <- log(mean_B) - log(mean_A)

    init_params <- c(
      qlogis(pA_init),
      qlogis(pB_init),
      log_rate_A_init,
      delta_rate_init,  # <-- instead of log(mean_B)
      0                 # log(phi) init
    )

    offset_long <- rep(offset, times = n_row)
    opt_global <- optim(
      par    = init_params,
      fn     = nll_hurdle,
      counts = df_long$count,
      groups = df_long$group,
      method = "BFGS",
      priors = NULL,
      offset = offset_long
    )

    # Extract and store bin-level results
    local_prob_A[j]  <- inv_logit(opt_global$par[1])
    local_prob_B[j]  <- inv_logit(opt_global$par[2])
    local_log_phi[j] <- opt_global$par[5]
  }

  # Map genes to bins
  gene_bin <- rep(NA_integer_, length(logmeans))
  for (j in seq_len(n_bins)) {
    idx <- finite & logmeans >= bin_lower[j] & logmeans < bin_upper[j]
    gene_bin[idx] <- j
  }

  # Construct the probability matrix (n_genes x 2)
  # Rows correspond to genes, Columns to Group A and Group B priors
  prob_matrix <- matrix(NA_real_, nrow = length(logmeans), ncol = 2)
  colnames(prob_matrix) <- c("prob_A", "prob_B")

  prob_matrix[, 1] <- local_prob_A[gene_bin]
  prob_matrix[, 2] <- local_prob_B[gene_bin]

  # Save to the object
  y$prob_matrix        <- prob_matrix
  y$prior_log_phi_gene <- local_log_phi[gene_bin]

  y$prior_bins <- data.frame(
    lower   = bin_lower,
    upper   = bin_upper,
    prob_A  = local_prob_A,
    prob_B  = local_prob_B,
    log_phi = local_log_phi
  )

  return(y)
}


#' Negative Log-Likelihood for Hurdle Model with Fixed Zero Probabilities
#'
#' Internal function to compute the negative log-likelihood of a hurdle negative
#' binomial model where the zero (dropout) probabilities are fixed. Only the mean
#' and dispersion parameters are optimized.
#'
#' @param params Numeric vector of length 5:
#'   \describe{
#'     \item{1}{logit(theta_A): 1 - zero probability for group A (fixed)}
#'     \item{2}{logit(theta_B): 1 - zero probability for group B (fixed)}
#'     \item{3}{log(mu_A): mean parameter for group A}
#'     \item{4}{delta_rate: log mean difference between groups}
#'     \item{5}{log(phi): dispersion parameter (log scale)}
#'   }
#' @param counts Numeric vector of observed counts for a single gene.
#' @param groups Factor indicating group membership for each observation.
#' @param offset Optional numeric vector of normalization factors. If \code{NULL},
#' defaults to 1 for all observations.
#' @param priors Optional numeric vector of length 2 specifying the mean and
#' standard deviation for a normal prior on \code{log(phi)}. If \code{NULL},
#' no prior is applied.
#'
#' @return
#' A single numeric value: the negative log-likelihood.
#'
#' @details
#' This function is used internally during tag-wise dispersion estimation, where
#' zero probabilities are fixed based on bin-level prior estimates. The likelihood
#' consists of:
#' \itemize{
#'   \item A Bernoulli component modeling zero vs non-zero counts
#'   \item A truncated negative binomial component for positive counts
#'   \item An optional normal prior on \code{log(phi)}
#' }
#'
#' @keywords internal
#' @noRd
nll_hurdle_fixed_P <- function(params, counts, groups, offset = NULL, priors = NULL){
  eps <- 1e-12

  # --- 0. Handle Offset ---
  if (is.null(offset)) {
    offset <- rep(1, length(counts))
  }

  # --- 1. Unpack Parameters ---
  logit_theta_A <- params[1]
  logit_theta_B <- params[2]
  log_rate_A    <- params[3]
  delta_rate    <- params[4]
  log_rate_B    <- log_rate_A + delta_rate
  # params[5] is log_phi. We exponentiate to get phi (size parameter)
  log_phi       <- params[5]
  phi           <- exp(log_phi)

  theta_A <- plogis(logit_theta_A)
  theta_B <- plogis(logit_theta_B)

  # --- 2. Assign Parameters to Data Points ---
  levs <- levels(groups)
  thetas <- ifelse(groups == levs[1], theta_A, theta_B)
  thetas <- pmin(pmax(thetas, eps), 1 - eps)

  base_rates <- ifelse(groups == levs[1], exp(log_rate_A), exp(log_rate_B))
  mus <- base_rates * offset

  # --- 3. Likelihood Calculation ---
  is_zero <- (counts == 0)

  # Binary part: Probability of being non-zero vs zero
  # ll_binary <- sum(log(1 - thetas[is_zero])) + sum(log(thetas[!is_zero]))
  ll_binary <- sum(log1p(-thetas[is_zero])) + sum(log(thetas[!is_zero]))

  # Truncated NB part
  pos_counts <- counts[!is_zero]
  pos_mus    <- mus[!is_zero]

  # Standard NB log-likelihood for positive counts
  ll_nb <- sum(dnbinom(pos_counts, mu = pos_mus, size = phi, log = TRUE))

  # Truncation adjustment: log(1 - P(Y=0))
  prob_zero_nb <- dnbinom(0, mu = pos_mus, size = phi)

  # Numerical safety: ensure prob_zero_nb doesn't exactly hit 1
  prob_zero_nb <- pmin(pmax(prob_zero_nb, 0), 1 - eps)
  trunc_adjustment <- -sum(log1p(-prob_zero_nb))

  ll_count <- ll_nb + trunc_adjustment

  # --- 4. Prior part ---
  # Prior is typically applied to the log-dispersion parameter
  if (is.null(priors)) {
    prior_val <- 0
  } else {
    # priors[1] is the mean from the trend (prior_log_phi_gene)
    # priors[2] is the width (SD) of the prior
    prior_val <- dnorm(log_phi, mean = priors[1], sd = priors[2], log = TRUE)
  }

  total_ll <- ll_binary + ll_count + prior_val

  return(-total_ll)
}

#' Tag-Wise Dispersion Estimation for Small Samples
#'
#' Internal function that estimates gene-specific dispersions by fixing the
#' bin-level zero probabilities obtained from \code{\link{priorEst}} and
#' optimizing the group means and dispersion for each gene.
#'
#' @param y A list object created by \code{\link{prepareDGE}} with estimated
#'   size factors, containing counts and sample information.
#' @param prior Logical. Whether to use prior information for dispersion
#'   shrinkage.
#'
#' @return The input \code{y} object augmented with \code{tagwise.disp} and
#'   \code{zero_prob_matrix}.
#'
#' @keywords internal
#' @noRd
tagwiseEst.smallSample <- function(y, prior=TRUE)  {
  # 1. Get the bin-based priors for each gene
  y2 <- priorEst(y)
  n_row <- nrow(y2$counts)

  tagwise.dispersions <- rep(NA_real_, n_row)

  for (i in seq_len(n_row)) {
    # Prepare data for current gene
    my_counts <- as.vector(t(y2$counts[i, ]))
    my_groups <- factor(y2$samples$group)
    levs <- levels(my_groups)

    # 2. Retrieve FIXED zero probabilities from the priorEst step
    fixed_prob_A <- y2$prob_matrix[i, "prob_A"]
    fixed_prob_B <- y2$prob_matrix[i, "prob_B"]

    data_A_pos <- my_counts[my_groups == levs[1] & my_counts > 0]
    data_B_pos <- my_counts[my_groups == levs[2] & my_counts > 0]

    # 3. Initial values for the REMAINING parameters only
    mean_A <- if(length(data_A_pos) > 0) mean(data_A_pos) else mean(my_counts + 0.1)
    mean_B <- if(length(data_B_pos) > 0) mean(data_B_pos) else mean(my_counts + 0.1)

    # Reduced parameter vector: log_mu_A, log_mu_B, log_phi
    log_rate_A_init <- log(mean_A)
    delta_rate_init <- log(mean_B) - log(mean_A)
    init_params_reduced <- c(
      log_rate_A_init,
      delta_rate_init,
      0
    )

    prior_mean_phi <- y2$prior_log_phi_gene[i]

    # 4. Optimization
    # Note: We wrap nll_hurdle to keep prob_A and prob_B constant
    opt_prime <- optim(
      par = init_params_reduced,
      fn = function(p) {
        # Reconstruct the full 5-parameter vector for the original nll_hurdle
        # Full vector: [logit_theta_A, logit_theta_B, log_mu_A, log_mu_B, log_phi]
        full_params <- c(
          qlogis(fixed_prob_A),
          qlogis(fixed_prob_B),
          p[1], # log_mu_A
          p[2], # delta_rate
          p[3]  # log_phi
        )
        nll_hurdle_fixed_P(
          full_params,
          counts = my_counts,
          groups = my_groups,
          priors = if (prior) c(prior_mean_phi, 3) else NULL,
          offset = y2$samples$size.factor
        )
      },
      method = "BFGS"
    )

    tagwise.dispersions[i] <- exp(-opt_prime$par[3])
  }

  y$tagwise.disp <- tagwise.dispersions
  # Optionally keep the priors for reference
  # y$prob_matrix <- y2$prob_matrix
  y$zero_prob_matrix <- 1 - y2$prob_matrix

  return(y)
}

#' Tag-Wise Dispersion Estimation for Large Samples
#'
#' Internal function that fits a hurdle negative binomial model separately to
#' each gene, jointly estimating the group-specific nonzero probabilities, mean
#' parameters, and dispersion without borrowing information across genes or
#' applying prior-based shrinkage.
#'
#' @param y A list object created by \code{\link{prepareDGE}} with estimated
#'   size factors, containing counts and sample information.
#'
#' @return The input \code{y} object augmented with \code{tagwise.disp} and
#'   \code{zero_prob_matrix}.
#'
#' @keywords internal
#' @noRd
tagwiseEst.largeSample <- function(y){
  counts <- as.matrix(y$counts)
  offset <- y$samples$size.factor
  grp    <- as.factor(y$samples$group)

  inv_logit <- function(x) 1 / (1 + exp(-x))

  n_genes <- nrow(counts)

  local_log_phi <- rep(NA_real_, n_genes)
  local_prob_A  <- rep(NA_real_, n_genes)
  local_prob_B  <- rep(NA_real_, n_genes)

  for (j in seq_len(n_genes)) {

    # Use only the current row
    genes <- counts[j, , drop = FALSE]

    df_long <- data.frame(
      group = grp,
      count = as.vector(t(genes))
    )

    levs <- levels(df_long$group)

    # Initial parameter estimation
    pA_init <- max(
      0.01,
      min(0.99, mean(
        df_long$count[df_long$group == levs[1]] > 0
      ))
    )

    pB_init <- max(
      0.01,
      min(0.99, mean(
        df_long$count[df_long$group == levs[2]] > 0
      ))
    )

    mean_A <- mean(
      df_long$count[
        df_long$group == levs[1] & df_long$count > 0
      ]
    )

    mean_B <- mean(
      df_long$count[
        df_long$group == levs[2] & df_long$count > 0
      ]
    )

    if (is.nan(mean_A)) mean_A <- mean(df_long$count)
    if (is.nan(mean_B)) mean_B <- mean(df_long$count)

    log_rate_A_init <- log(mean_A)
    delta_rate_init <- log(mean_B) - log(mean_A)

    init_params <- c(
      qlogis(pA_init),
      qlogis(pB_init),
      log_rate_A_init,
      delta_rate_init,
      0
    )

    opt_global <- optim(
      par     = init_params,
      fn      = nll_hurdle,
      counts  = df_long$count,
      groups  = df_long$group,
      method  = "BFGS",
      priors  = NULL,
      offset  = offset
    )

    # Store row-specific estimates
    local_prob_A[j]  <- inv_logit(opt_global$par[1])
    local_prob_B[j]  <- inv_logit(opt_global$par[2])
    local_log_phi[j] <- opt_global$par[5]
  }

  prob_matrix <- cbind(
    prob_A = local_prob_A,
    prob_B = local_prob_B
  )

  # rownames(prob_matrix) <- rownames(counts)

  y$zero_prob_matrix   <- 1- prob_matrix
  y$tagwise.disp <- exp(-local_log_phi)

  return(y)
}

#' Tag-wise Dispersion Estimation for Hurdle Negative Binomial Model
#'
#' Estimate gene-specific (tag-wise) dispersion parameters for a hurdle negative
#' binomial model using prior information derived from bin-level estimates.
#'
#' @param y A list object created by \code{\link{prepareDGE}} with size factors
#'   estimated, containing counts and sample information.
#' @param small.sample Logical. If \code{TRUE}, borrows information from genes
#'   with similar expression patterns. If \code{FALSE}, estimates dispersion
#'   using each gene independently. Setting \code{FALSE} is recommended for
#'   adaptive single-cell differential expression analysis.
#' @param prior Logical. Whether to use prior information for dispersion
#'   shrinkage when \code{small.sample = TRUE}. This argument is ignored when
#'   \code{small.sample = FALSE}.
#'
#' @return
#' The input \code{y} object augmented with:
#' \describe{
#'   \item{tagwise.disp}{Numeric vector of estimated gene-wise dispersions.}
#'   \item{zero_prob_matrix}{Numeric matrix of fixed zero probabilities for each
#'         gene and group.}
#' }
#'
#' @details
#' When \code{small.sample = TRUE}, the function calls
#' \code{tagwiseEst.smallSample}, which performs the following steps:
#' \enumerate{
#'   \item Retrieves bin-level prior estimates of zero probabilities and log-dispersion
#'         for each gene via \code{\link{priorEst}}.
#'   \item Fixes the zero probabilities and optimizes only the mean parameters
#'         and dispersion for each gene individually.
#'   \item Uses the internal function \code{nll_hurdle_fixed_P} to compute the negative
#'         log-likelihood with fixed zero probabilities.
#' }
#'
#' When \code{small.sample = FALSE}, the function calls
#' \code{tagwiseEst.largeSample}, which fits a hurdle negative binomial model
#' separately to each gene without borrowing information across genes.
#'
#' The resulting \code{tagwise.disp} will be used for downstream differential
#' expression analysis.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' head(y$tagwise.disp)
#' head(y$zero_prob_matrix)
#'
#' @export
tagwiseEst <- function(y, small.sample = TRUE, prior = TRUE) {
  if (small.sample) {
    return(tagwiseEst.smallSample(y, prior = prior))
  }

  return(tagwiseEst.largeSample(y))
}
