# The data frame form, shaped after metafor::conv.2x2.
#
# conv.2x2 converts a summary into one table and appends four columns of cell
# counts. enum.2x2 enumerates every table the summary admits and appends the
# outcome, the number of compatible tables, and the span of each cell across
# them; the four point columns are filled only where the printed digits
# determine the table, and are NA otherwise. Swapping one call for the other
# therefore shows what the enumeration buys: the same four columns, plus what
# they leave open.

.getx <- function(nm, mf, data, enclos) {
  if (is.null(mf[[nm]])) return(NULL)
  eval(mf[[nm]], data, enclos)
}

.as_character_col <- function(x, nm) {
  if (is.null(x)) return(NULL)
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) {
    stop_invalid_input(sprintf(
      paste0("%s must be the strings the sources printed, not %s; ",
             "a float cannot distinguish '0.10' from '0.1'"), nm, typeof(x)))
  }
  x
}

#' Enumerate the 2x2 tables a column of published statistics allows
#'
#' The data frame form of [recover2x2()], shaped after [metafor::conv.2x2()] so
#' that one call can be swapped for the other. `conv.2x2` reconstructs a single
#' table per row by optimization and appends four columns of cell counts;
#' `enum.2x2` enumerates every integer table whose statistics round back to the
#' printed digits and appends the outcome, the number of compatible tables, and
#' the span of each cell across them. The four point columns are filled only
#' where the digits determine the table and are `NA` otherwise, so the
#' difference between the two functions is visible in the same columns.
#'
#' Cell and marginal names are those of `conv.2x2`: `ai` positive on both
#' criteria, `bi` positive on the first only, `ci` positive on the second only,
#' `di` negative on both, with `n1i = ai + bi`, `n2i = ai + ci` and
#' `ni = ai + bi + ci + di`. `conv.2x2` takes an odds ratio, a phi coefficient
#' or a chi-square statistic; this takes unweighted Cohen's kappa, which is what
#' agreement studies report, as the string the source printed.
#'
#' Arguments are evaluated inside `data` when it is given, as in `conv.2x2`.
#' No row signals: a row whose figures cannot describe any table comes back with
#' status `"impossible"` and a reason.
#'
#' @param kappa Unweighted Cohen's kappa for each row, as the strings the
#'   sources printed. Character, never numeric.
#' @param ni Vector of total sample sizes.
#' @param n1i,n2i Vectors of marginal counts for the outcome of interest on the
#'   first and second variable.
#' @param p1i,p2i Optional vectors of marginals as the strings the sources
#'   printed, for rows where the count was not reported. Give `n1i` or `p1i` for
#'   a row, not both.
#' @param agreement Optional vector of observed agreement, as printed. Where a
#'   report prints it alongside kappa it enters as a further constraint.
#' @param marginals.as.percent,agreement.as.percent Logical, recycled over rows;
#'   whether those printed figures are percentages.
#' @param data Optional data frame containing the variables given above.
#' @param include Optional logical or numeric vector giving the subset of rows
#'   to enumerate. Rows left out come back with `status` `NA`.
#' @param var.names Character vector of four names for the reconstructed cells.
#'   The range columns are these names with `.min` and `.max` appended.
#' @param append Logical; whether to return `data` together with the new
#'   columns.
#' @param replace Character string or logical controlling how existing columns
#'   named in `var.names` are treated, as in [metafor::conv.2x2()]. With
#'   `"ifna"` (or `FALSE`) only missing values are filled; with `"all"` (or
#'   `TRUE`) every uniquely determined cell is written.
#' @param ... Unused, for compatibility.
#' @return A data frame. Where `data` was given and `append = TRUE` it is the
#'   input frame with the new columns; otherwise it is the new columns alone.
#'   The list of [recover2x2()] objects, one per row, is attached as the
#'   `"recoveries"` attribute, so the full set of tables for any row remains
#'   reachable.
#' @seealso [recover2x2()] for one comparison at a time, and
#'   [metafor::conv.2x2()] for the point reconstruction.
#' @examples
#' dat <- data.frame(
#'   study = c("delirium pooled", "cohort B", "cohort C"),
#'   ni    = c(768, 20306, 370),
#'   n1i   = c(158, 866, 165),
#'   n2i   = c(466, 1603, 160),
#'   kappa = c("0.29", "0.22", "0.48"))
#'
#' out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat)
#' out[, c("study", "status", "ntables", "ai", "ai.min", "ai.max")]
#'
#' # The full set for any row is still reachable.
#' attr(out, "recoveries")[[2]]
#'
#' # Scalar form, for one comparison.
#' enum.2x2(kappa = "0.29", ni = 768, n1i = 158, n2i = 466)
#' @export
enum.2x2 <- function(kappa, ni, n1i, n2i, p1i, p2i, agreement,
                     marginals.as.percent = FALSE, agreement.as.percent = FALSE,
                     data, include,
                     var.names = c("ai", "bi", "ci", "di"),
                     append = TRUE, replace = "ifna", ...) {
  mf <- match.call()
  has.data <- !missing(data) && !is.null(data)
  enclos <- parent.frame()
  env <- if (has.data) data else list()

  kappa <- .as_character_col(.getx("kappa", mf, env, enclos), "kappa")
  ni <- .getx("ni", mf, env, enclos)
  n1i <- .getx("n1i", mf, env, enclos)
  n2i <- .getx("n2i", mf, env, enclos)
  p1i <- .as_character_col(.getx("p1i", mf, env, enclos), "p1i")
  p2i <- .as_character_col(.getx("p2i", mf, env, enclos), "p2i")
  agreement <- .as_character_col(.getx("agreement", mf, env, enclos), "agreement")
  include <- .getx("include", mf, env, enclos)

  if (is.null(ni)) stop_invalid_input("ni must be given")
  if (length(var.names) != 4L || !is.character(var.names) || anyDuplicated(var.names)) {
    stop_invalid_input("var.names must be four distinct names")
  }

  k <- max(length(ni), length(kappa), length(n1i), length(n2i),
           length(p1i), length(p2i), length(agreement))
  grow <- function(x) if (is.null(x)) NULL else rep_len(x, k)
  ni <- grow(ni); kappa <- grow(kappa)
  n1i <- grow(n1i); n2i <- grow(n2i)
  p1i <- grow(p1i); p2i <- grow(p2i); agreement <- grow(agreement)
  mpct <- rep_len(marginals.as.percent, k)
  apct <- rep_len(agreement.as.percent, k)

  keep <- rep(TRUE, k)
  if (!is.null(include)) {
    if (is.logical(include)) {
      keep <- rep_len(include, k)
      keep[is.na(keep)] <- FALSE
    } else {
      keep <- rep(FALSE, k)
      keep[include[!is.na(include)]] <- TRUE
    }
  }

  rows <- lapply(seq_len(k), function(i) {
    if (!keep[i]) return(NULL)
    list(ni = ni[i], n1i = n1i[i], n2i = n2i[i], p1i = p1i[i], p2i = p2i[i],
         kappa = kappa[i], agreement = agreement[i],
         marginals.as.percent = mpct[i], agreement.as.percent = apct[i])
  })
  recoveries <- lapply(rows, function(r) if (is.null(r)) NULL else recover_many(list(r))[[1L]])

  na4 <- function() rep(NA_real_, k)
  out <- data.frame(status = rep(NA_character_, k), ntables = rep(NA_integer_, k),
                    stringsAsFactors = FALSE)
  for (cell in var.names) {
    out[[cell]] <- na4()
    out[[paste0(cell, ".min")]] <- na4()
    out[[paste0(cell, ".max")]] <- na4()
  }
  out$reason <- rep(NA_character_, k)

  for (i in seq_len(k)) {
    r <- recoveries[[i]]
    if (is.null(r)) next
    out$status[i] <- r$status
    out$ntables[i] <- nrow(r$tables)
    out$reason[i] <- r$reason
    if (nrow(r$tables) == 0L) next
    for (j in seq_along(var.names)) {
      col <- c("ai", "bi", "ci", "di")[j]
      out[[paste0(var.names[j], ".min")]][i] <- min(r$tables[[col]])
      out[[paste0(var.names[j], ".max")]][i] <- max(r$tables[[col]])
      if (nrow(r$tables) == 1L) out[[var.names[j]]][i] <- r$tables[[col]][1L]
    }
  }

  if (has.data && isTRUE(append)) {
    res <- data
    all_at <- isTRUE(replace) || identical(replace, "all")
    for (nm in names(out)) {
      if (nm %in% var.names && nm %in% names(res) && !all_at) {
        fill <- is.na(res[[nm]]) & !is.na(out[[nm]])
        res[[nm]][fill] <- out[[nm]][fill]
      } else if (nm %in% var.names && nm %in% names(res)) {
        fill <- !is.na(out[[nm]])
        res[[nm]][fill] <- out[[nm]][fill]
      } else {
        res[[nm]] <- out[[nm]]
      }
    }
    out <- res
  }
  attr(out, "recoveries") <- recoveries
  out
}
