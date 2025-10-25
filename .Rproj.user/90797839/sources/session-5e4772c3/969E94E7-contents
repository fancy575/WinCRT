#' @useDynLib WinCRT, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
#' Win Ratio for Cluster-Randomized Trials (Hierarchical or Single Endpoint)
#'
#' Computes win/loss/tie tallies and win-ratio statistics in CRTs using a
#' hierarchical endpoint (tiers) or a single endpoint.
#' performs the pairwise counting; this R wrapper handles validation and summaries.
#'
#' @param data A data.frame/tibble in \emph{long} format (possibly multiple rows per subject).
#' @param id Name of subject ID column. Default: "id".
#' @param trt Name of treatment/arm indicator column (two arms). Default: "z".
#' @param cluster Name of cluster ID column. Default: "cluster".
#' @param outcome Name of the numeric outcome column. Default: "outcome".
#'   For survival-style data, this is a time; for a single endpoint, it is the scalar outcome.
#' @param tier Name of the tier code column. Default: "tier". Use \code{0} for censoring rows
#'   (survival-style), and \code{1..K} for event tiers with larger meaning higher priority.
#'   For a single endpoint, set \code{tier=1} for all rows and provide no censoring rows.
#' @param lower_better Logical; if \code{TRUE}, smaller outcome values mean better.
#'   The scale is flipped internally so that the kernel’s “earlier is worse” rule yields
#'   the correct win/loss direction. Default: \code{FALSE}.
#'
#' @return An object of class \code{"WinCRT"} with:
#' \itemize{
#'   \item \code{estimates}: \code{W_D}, \code{logW_R}, \code{logW_O}, \code{p_ties}
#'   \item \code{variances}: \code{VarW_D}, \code{VarlogW_R}, \code{VarlogW_O}
#'   \item \code{rho}: rank-based within-cluster correlation
#'   \item \code{counts}: totals across arms: \code{wins}, \code{losses}, \code{ties}, \code{T0}, \code{ties_cnt}
#'   \item \code{per_id}: per-subject table: \code{id, cluster, z, w, l, t, tie_diff, rank}
#'   \item \code{per_cluster}: per-cluster table with score \code{S} and \code{z}
#' }
#'
#' @examples
#' \dontrun{
#'   # Survival-style with tiers (0=censor, 1=hosp, 2=death)
#'   res <- wincrt(df, id="id", trt="z", cluster="cluster",
#'                 outcome="time", tier="event")
#'
#'   # Single endpoint (no censoring), smaller is better
#'   df$tier <- 1L
#'   res <- wincrt(df, outcome="score", tier="tier", lower_better=TRUE)
#' }
#'
#' @import dplyr Rcpp
#' @importFrom rlang .data !! sym
#' @export
wincrt <- function(data,
                   id,
                   trt,
                   cluster,
                   outcome,
                   tier,
                   lower_better = FALSE) {



  # ---------- validations ----------
  if (!is.data.frame(data)) data.frame(data)
  nm <- names(data)
  req <- c(id, trt, cluster, outcome, tier)
  miss <- setdiff(req, nm)
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

  df <- as.data.frame(data)

  # Canonical names for C++: id, z, y, tier
  df <- dplyr::rename(
    df,
    id      = !!rlang::sym(id),
    z       = !!rlang::sym(trt),
    cluster = !!rlang::sym(cluster),
    y       = !!rlang::sym(outcome),
    tier    = !!rlang::sym(tier)
  )

  if (any(is.na(df$id))) stop("`id` contains NA.")
  if (any(is.na(df$z)))  stop("`z` contains NA.")
  if (dplyr::n_distinct(df$z) != 2) stop("`z` must have exactly two distinct values (two arms).")
  if (!is.numeric(df$y))  stop("`outcome` must be numeric.")
  if (!all(df$tier >= 0, na.rm = TRUE))
    stop("`tier` must be >= 0 (0=censor, 1..K=tiers with larger=more important).")

  # Flip scale if smaller outcome is better
  if (isTRUE(lower_better)) {
    df$y <- -df$y
  }

  # ---------- per-subject unique table (no Nc) ----------
  per_id <- df |>
    dplyr::distinct(.data$id, .data$cluster, .data$z)

  n <- nrow(per_id)
  if (n < 2) stop("Need at least 2 distinct subjects after preprocessing.")

  # Keep the id order stable for mapping kernel outputs to per_id rows
  # Ensure df (long) rows are not required to be sorted the same way:
  # wr_hier is expected to produce a 4*n vector in the per_id order.
  # If your wr_hier relies on first-appearance order, make sure df is ordered by per_id id order.
  # We enforce that by matching:
  id_index <- match(df$id, per_id$id)
  if (any(is.na(id_index))) stop("Internal mapping error: some long rows do not map to distinct subjects.")

  # ---------- call C++ kernel ----------
  # wr_hier requires columns: id, z, y, tier
  ans <- wr_hier(df[, c("id", "z", "y", "tier")])
  if (length(ans) != 4L * n) {
    stop("Internal error: wr_hier returned length ", length(ans), " but expected ", 4L * n, ".")
  }

  w        <- ans[seq_len(n)]
  l        <- ans[(n + 1):(2 * n)]
  t_all    <- ans[(2 * n + 1):(3 * n)]
  tie_diff <- ans[(3 * n + 1):(4 * n)]

  per_id$w <- as.numeric(w)
  per_id$l <- as.numeric(l)
  per_id$t <- as.numeric(t_all)
  per_id$tie_diff <- as.numeric(tie_diff)
  per_id$rank <- per_id$w + 1 + per_id$t / 2

  # ---------- arm mapping ----------
  z_vals <- sort(unique(per_id$z))
  arm_C  <- z_vals[1]
  arm_T  <- z_vals[2]
  ti <- which(per_id$z == arm_T)
  ci <- which(per_id$z == arm_C)

  N_T <- length(ti); N_C <- length(ci)
  if (N_T == 0 || N_C == 0) stop("Both arms must be present after preprocessing.")

  # Cross-arm totals
  ties_cnt <- sum(per_id$tie_diff[ti])
  T0       <- as.numeric(N_T) * as.numeric(N_C)

  W_D_sum  <- sum((per_id$w - per_id$l)[ti])
  Wval     <- (T0 - ties_cnt + W_D_sum) / 2
  Lval     <- (T0 - ties_cnt - W_D_sum) / 2
  if (Wval < 0) Wval <- 0
  if (Lval < 0) Lval <- 0

  p_ties <- ties_cnt / T0
  W_D    <- W_D_sum / T0

  # logs with small guard if needed
  if (Lval == 0 || Wval == 0) {
    eps <- .Machine$double.eps
    logW_R <- log((Wval + eps) / (Lval + eps))
    logW_O <- log(((Wval + 0.5 * ties_cnt) + eps) / ((Lval + 0.5 * ties_cnt) + eps))
  } else {
    logW_R <- log(Wval / Lval)
    logW_O <- log((Wval + 0.5 * ties_cnt) / (Lval + 0.5 * ties_cnt))
  }

  # ---------- per-cluster scores (no Nc) ----------
  S_tbl <- per_id |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      S = sum(.data$w - .data$l),
      z = dplyr::first(.data$z),
      .groups = "drop"
    )

  M1 <- sum(S_tbl$z == arm_T)
  M0 <- sum(S_tbl$z == arm_C)
  M  <- M1 + M0
  if (M1 == 0 || M0 == 0) stop("Both arms must be present at the cluster level.")
  if (M1 < 2 || M0 < 2) warning("At least one arm has fewer than 2 clusters; variance estimates may be unstable/undefined.")

  # Unbiased per-arm cluster variance (your ALT formula; does NOT use Nc)
  # VarW_D = ((M1*M0)/(M*N_T*N_C))^2 * ( SS_T/(M1*(M1-1)) + SS_C/(M0*(M0-1)) )
  S_T <- S_tbl$S[S_tbl$z == arm_T]; S_T_bar <- mean(S_T); SS_T <- sum((S_T - S_T_bar)^2)
  S_C <- S_tbl$S[S_tbl$z == arm_C]; S_C_bar <- mean(S_C); SS_C <- sum((S_C - S_C_bar)^2)

  VarW_D <- if (M1 > 1 && M0 > 1) {
    ((M1 * M0) / (M * as.numeric(N_T) * as.numeric(N_C)))^2 *
      ( SS_T / (M1 * (M1 - 1)) + SS_C / (M0 * (M0 - 1)) )
  } else NA_real_

  VarlogW_R <- if (!is.na(VarW_D)) {
    (2 / (1 - p_ties) / (1 - (W_D / (1 - p_ties))^2))^2 * VarW_D
  } else NA_real_

  VarlogW_O <- if (!is.na(VarW_D)) {
    (2 / (1 - W_D^2))^2 * VarW_D
  } else NA_real_

  # ---------- rho from ranks ----------
  n_subj <- nrow(per_id)
  per_id$weight <- 1 / n_subj
  F_bar <- mean(per_id$rank)
  denom <- sum(per_id$weight * (per_id$rank - F_bar)^2)

  rank_sum <- per_id |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      product = {
        r <- rank - F_bar
        m <- length(r)
        if (m < 2) 0 else 2 * sum(utils::combn(r, 2, prod)) / (m * (m - 1)) * sum(weight)
      },
      .groups = "drop"
    )
  rho <- if (denom > 0) sum(rank_sum$product) / denom else NA_real_

  # ---------- assemble ----------
  res <- list(
    estimates = list(
      W_D    = W_D,
      logW_R = logW_R,
      logW_O = logW_O,
      p_ties = p_ties
    ),
    variances = list(
      VarW_D     = VarW_D,
      VarlogW_R  = VarlogW_R,
      VarlogW_O  = VarlogW_O
    ),
    rho = rho,
    counts = list(
      wins     = as.numeric(Wval),
      losses   = as.numeric(Lval),
      ties     = as.numeric(ties_cnt),
      T0       = as.numeric(T0),
      ties_cnt = as.numeric(ties_cnt)
    ),
    per_id = per_id[, c("id","cluster","z","w","l","t","tie_diff","rank")],
    per_cluster = S_tbl
  )
  class(res) <- "WinCRT"
  res
}

#' Summarize a WinCRT fit (prints results)
#'
#' Presents inferential results for a \code{WinCRT} object, including the estimate,
#' standard error, test statistic, p-value, and confidence interval for a chosen
#' estimand under either a t- or z-reference distribution. The output is printed
#' in a compact table; the function returns \code{NULL} invisibly.
#'
#' @description
#' This method computes Wald-style inference for win-based estimands from
#' \code{\link{wincrt}}. You can request one estimand or all of them. Confidence
#' intervals are shown in a single column as \code{"(low, high)"}. All numeric
#' values are formatted to three significant digits. When \code{test = "t"}, the
#' reference uses \eqn{M-2} degrees of freedom, where \eqn{M} is the total number
#' of clusters across arms; \code{test = "z"} uses the standard normal.
#'
#' @param object A fitted object of class \code{WinCRT}, returned by \code{\link{wincrt}}.
#' @param test Reference distribution: \code{"t"} (default; df = total clusters minus 2)
#'   or \code{"z"} (standard normal).
#' @param estimand Target parameter to display (default \code{"logWR"}):
#'   \code{"logWR"} (log win ratio), \code{"logWO"} (log win odds),
#'   \code{"WD"} (net win proportion), or \code{"all"} to display all three.
#' @param alpha Significance level for confidence intervals, a number in \code{(0, 1)}
#'   (default \code{0.05} for 95\% CIs).
#' @param alternative Alternative hypothesis for the Wald test:
#'   \code{"two.sided"} (default), \code{"greater"} (parameter > 0),
#'   or \code{"less"} (parameter < 0).
#' @param ... Unused; included for S3 compatibility.
#'
#' @details
#' Null hypotheses are \eqn{\mathrm{logWR}=0}, \eqn{\mathrm{logWO}=0}, and \eqn{\mathrm{WD}=0}.
#' Reported counts (clusters and subjects per arm, total wins/losses/ties, number of
#' cross-arm pairs, tie proportion, and rank ICC) are taken from the fitted object.
#'
#' @return Invisibly returns \code{NULL}. Results are printed to the console.
#'
#' @seealso \code{\link{wincrt}} for fitting.
#'
#' @examples
#' \dontrun{
#' fit <- wincrt(df, id = "id", trt = "z", cluster = "cluster",
#'               outcome = "outcome", tier = "tier")
#'
#' ## Default estimand is logWR:
#' summary(fit)
#'
#' ## All estimands with z-reference and one-sided alternative:
#' summary(fit, test = "z", estimand = "all", alpha = 0.025, alternative = "greater")
#' }
#'
#' @export
#' @method summary WinCRT

summary.WinCRT <- function(object,
                           test = c("t", "z"),
                           estimand = c("logWR", "logWO", "WD", "all"),
                           alpha = 0.05,
                           alternative = c("two.sided", "greater", "less"),
                           ...) {
  stopifnot(inherits(object, "WinCRT"))
  test        <- match.arg(test)
  estimand    <- match.arg(estimand)   # default "logWR"
  alternative <- match.arg(alternative)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0,1).")

  ests <- object$estimates
  vars <- object$variances
  per_id <- object$per_id
  per_cluster <- object$per_cluster

  # cluster counts (df used internally for p-values but not printed)
  z_vals <- sort(unique(per_cluster$z))
  arm_C  <- z_vals[1]; arm_T <- z_vals[2]
  M1 <- sum(per_cluster$z == arm_T)
  M0 <- sum(per_cluster$z == arm_C)
  M  <- M1 + M0
  df <- max(M - 2, 1)

  # helpers
  qfun <- switch(test,
                 "t" = function(p) stats::qt(p, df = df),
                 "z" = function(p) stats::qnorm(p))
  pfun <- switch(test,
                 "t" = function(x) stats::pt(x, df = df, lower.tail = FALSE),
                 "z" = function(x) stats::pnorm(x, lower.tail = FALSE))

  build_row <- function(name, est, var) {
    if (is.na(var)) {
      return(data.frame(Estimand = name, Estimate = est, SE = NA_real_,
                        `p-value` = NA_real_, CI = NA_character_,
                        check.names = FALSE))
    }
    se <- sqrt(var)
    stat <- est / se
    alpha2 <- if (alternative == "two.sided") alpha/2 else alpha
    crit <- qfun(1 - alpha2)

    if (alternative == "two.sided") {
      ci_low  <- est - crit * se
      ci_high <- est + crit * se
      pval <- 2 * pfun(abs(stat))
    } else if (alternative == "greater") {
      ci_low  <- est - crit * se
      ci_high <- Inf
      pval <- pfun(stat)
    } else {
      ci_low  <- -Inf
      ci_high <- est + crit * se
      pval <- pfun(-stat)
    }

    data.frame(Estimand = name, Estimate = est, SE = se,
               `p-value` = pval,
               CI = paste0("(", ci_low, ", ", ci_high, ")"),
               check.names = FALSE)
  }

  add_row <- function(which, acc) {
    if (which == "logWR") acc[[length(acc)+1]] <- build_row("logWR", ests$logW_R, vars$VarlogW_R)
    if (which == "logWO") acc[[length(acc)+1]] <- build_row("logWO", ests$logW_O, vars$VarlogW_O)
    if (which == "WD")    acc[[length(acc)+1]] <- build_row("WD",    ests$W_D,    vars$VarW_D)
    acc
  }
  rows <- list()
  if (estimand == "all") {
    rows <- add_row("logWR", rows); rows <- add_row("logWO", rows); rows <- add_row("WD", rows)
  } else {
    rows <- add_row(estimand, rows)
  }
  tab <- do.call(rbind, rows)

  # header (no df shown)
  counts <- object$counts
  n_per_arm <- table(per_id$z)
  N_T <- as.integer(n_per_arm[names(n_per_arm)==as.character(arm_T)])
  N_C <- as.integer(n_per_arm[names(n_per_arm)==as.character(arm_C)])

  cat("Win Ratio Summary (", toupper(test), "-test",
      ", alpha=", alpha, ", alternative=", alternative, ")\n", sep = "")
  cat("Clusters: M1=", M1, ", M0=", M0, ", M=", M, "\n", sep = "")
  cat("Subjects: n1=", N_T, ", n0=", N_C, "\n", sep = "")
  cat("Totals (between arms): Wins=", counts$wins, ", Losses=", counts$losses,
      ", Ties=", counts$ties_cnt,"\n", sep = "")
  cat("p_tie: ", formatC(ests$p_ties, digits = 3, format = "fg"),
      "; rho (rank ICC): ", formatC(object$rho, digits = 3, format = "fg"), "\n\n", sep = "")

  # 3-digit formatting for all numeric fields; CI stays in one column "(low, high)"
  fmt3 <- function(x) {
    if (is.infinite(x)) return("Inf")
    if (is.na(x)) return(NA_character_)
    formatC(x, digits = 3, format = "fg")
  }
  fmt_p <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 1e-3) formatC(p, digits = 3, format = "e") else formatC(p, digits = 3, format = "fg")
  }
  # re-build CI with formatted bounds
  tab$CI <- {
    get_bounds <- function(s) {
      if (is.na(s)) return(c(NA_real_, NA_real_))
      nums <- as.numeric(strsplit(gsub("[() ]", "", s), ",")[[1]])
      if (length(nums) != 2) c(NA_real_, NA_real_) else nums
    }
    bounds <- t(vapply(tab$CI, get_bounds, numeric(2)))
    paste0("(", fmt3(bounds[,1]), ", ", fmt3(bounds[,2]), ")")
  }
  tab$Estimate  <- vapply(tab$Estimate,  fmt3, character(1))
  tab$SE        <- vapply(tab$SE,        fmt3, character(1))
  tab$`p-value` <- vapply(tab$`p-value`, fmt_p, character(1))

  # final columns (no df, no Statistic)
  tab <- tab[, c("Estimand","Estimate","SE","p-value","CI")]
  print(tab, row.names = FALSE)
  invisible(NULL)
}
