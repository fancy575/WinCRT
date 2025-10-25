#' Ordinal single-endpoint CRT dataset: orData
#'
#' An example dataset for win–loss–tie analyses in cluster-randomized trials
#' with a **single ordinal endpoint** (no censoring, no hierarchical tiers).
#' Each row corresponds to one subject.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{cluster}{Cluster identifier (integer/factor).}
#'   \item{z}{Treatment arm indicator (two levels, e.g., 0 = control, 1 = treatment).}
#'   \item{Nc}{Cluster size (number of subjects in the subject’s cluster).}
#'   \item{y}{\strong{Ordinal outcome}. Store as an ordered numeric code reflecting category order
#'            (e.g., 1 < 2 < 3).}
#'   \item{id}{Subject identifier (integer).}
#' }
#'
#' @details
#' This dataset represents an \emph{ordinal} single-endpoint setting.
#'
#' @source Synthetic example for the \pkg{WinCRT} package.
#'
#' @usage data(orData)
#'
#' @examples
#' data(orData)
#' head(orData)
"orData"
