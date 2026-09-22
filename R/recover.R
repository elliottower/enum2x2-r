# Every integer 2x2 table consistent with what a source printed.

UNIQUE <- "unique"
SET <- "set"
INFEASIBLE <- "infeasible"
INSUFFICIENT <- "insufficient"
# A batch cannot signal, so a figure that cannot describe any table needs a
# status of its own: counting it as 'insufficient' would inflate the number of
# sources that under-reported.
IMPOSSIBLE <- "impossible"

.fmt <- function(x) {
  if (is.null(x)) return("NULL")
  out <- if (is.character(x)) {
    sprintf("'%s'", x)
  } else if (is.logical(x)) {
    ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE"))
  } else {
    format(x, scientific = FALSE, trim = TRUE)
  }
  paste(out, collapse = ", ")
}

.side_word <- c(a = "first", b = "second")
.side_count <- c(a = "n1i", b = "n2i")
.side_printed <- c(a = "p1i", b = "p2i")

# Validate one declared figure and return its parse. Mirrors the order of the
# checks in the Python package, because the message a caller sees is part of
# what a port has to reproduce.
.validate_decimal <- function(name, value) {
  p <- .dec_parse(value)
  if (!isTRUE(p$ok)) {
    stop_invalid_input(sprintf("%s = %s is not a decimal number", name, .fmt(value)))
  }
  if (!isTRUE(p$finite)) {
    stop_invalid_input(sprintf("%s = %s is not a finite number", name, .fmt(value)))
  }
  p
}

.marginal_counts <- function(n, count, printed, as.percent, side) {
  if (!is.null(count)) {
    if (is.logical(count) || !.is_whole_number(count)) {
      stop_invalid_input(sprintf("%s must be an integer count, got %s",
                                 .side_count[[side]], .fmt(count)))
    }
    if (count < 0 || count > n) {
      stop_invalid_input(sprintf("%s = %s is outside [0, %s]",
                                 .side_count[[side]], .fmt(count), .fmt(n)))
    }
    return(count)
  }
  b <- .printed_bounds_double(printed, as.percent)
  if (b[2L] < 0 || b[1L] > 1) {
    # No proportion rounds to this, so there is no marginal at all. That is a
    # defect in the call, not a property of the source, and reporting it as an
    # empty candidate set would attach a false reason to it.
    stop_invalid_input(sprintf(
      "%s = %s is outside [0, 1] as a proportion%s", .side_printed[[side]],
      .fmt(printed), if (isTRUE(as.percent)) " (given as a percentage)" else ""))
  }
  counts_rounding_to(printed, n, as.percent)
}

.empty_tables <- function() {
  data.frame(ai = numeric(0), bi = numeric(0), ci = numeric(0), di = numeric(0))
}

# The enumeration itself. For each candidate pair of marginal counts, the cells
# are a one-parameter family in ai, and membership is a pair of integer
# cross-multiplications evaluated over the whole family at once.
.enumerate <- function(n, a_counts, b_counts, kiv, giv) {
  pieces <- list()
  for (na in a_counts) {
    for (nb in b_counts) {
      chance <- na * nb + (n - na) * (n - nb)
      den <- n * n - chance
      # A candidate pair can make expected agreement exactly 1, leaving kappa
      # undefined for that pair alone. It contributes no table; it is not a
      # reason to abandon the row.
      if (den == 0) next
      lo <- max(0, na + nb - n)
      hi <- min(na, nb)
      if (hi < lo) next
      ai <- seq.int(lo, hi)
      agree <- n - na - nb + 2 * ai
      num <- agree * n - chance
      keep <- .prod_le(kiv$lo_num, den, kiv$den, num) &
        .prod_le(kiv$den, num, kiv$hi_num, den)
      if (!is.null(giv)) {
        keep <- keep & .prod_le(giv$lo_num, n, giv$den, agree) &
          .prod_le(giv$den, agree, giv$hi_num, n)
      }
      if (any(keep)) {
        k <- ai[keep]
        pieces[[length(pieces) + 1L]] <-
          data.frame(ai = k, bi = na - k, ci = nb - k, di = n - na - nb + k)
      }
    }
  }
  if (!length(pieces)) .empty_tables() else do.call(rbind, pieces)
}

#' Enumerate every 2x2 table a published agreement statistic still allows
#'
#' Given a sample size, both positive marginals and Cohen's kappa as the strings
#' a source printed, `recover2x2` returns every non-negative integer 2x2 table
#' whose statistics round back to those digits. Rounding is treated as a
#' constraint rather than as noise: a printed figure names a closed interval,
#' and membership in it is decided by integer cross-multiplication, so there is
#' no tolerance parameter.
#'
#' Cell names follow [metafor::conv.2x2()]. `ai` is the count positive on both
#' criteria, `bi` positive on the first and negative on the second, `ci`
#' negative on the first and positive on the second, `di` negative on both;
#' `n1i = ai + bi` and `n2i = ai + ci` are the two positive marginals and
#' `ni = ai + bi + ci + di` the sample size.
#'
#' Marginals are given either as exact counts (`n1i`, `n2i`) or as the strings a
#' source printed (`p1i`, `p2i`). A printed string carries its own precision, so
#' `"0.10"` and `"0.1"` are different inputs and figures must not be passed as
#' numbers; a numeric `kappa` is refused with that reason.
#'
#' The status is one of four. `"unique"` and `"set"` mean the figures admit one
#' table or several. `"infeasible"` means they are individually possible but
#' jointly admit none, and `reason` names the figure that excluded. Both are
#' properties of the source and are reported rather than signaled.
#' `"insufficient"` means a required figure was never published, which is a
#' different finding from `"infeasible"` and is kept apart from it. A defect in
#' the call --- a marginal above `ni`, a kappa outside \[-1, 1\], two marginals
#' for one criterion --- signals `enum2x2_invalid_input` instead.
#'
#' @param ni Sample size, a positive whole number.
#' @param n1i,n2i Positive marginal counts for the first and second criterion.
#' @param p1i,p2i Positive marginals as the strings the source printed, for a
#'   criterion whose count was not reported. Give `n1i` or `p1i`, not both.
#' @param kappa Unweighted Cohen's kappa, as the string the source printed.
#' @param agreement Observed agreement, as the string the source printed. Where
#'   a report prints it alongside kappa it enters as a further constraint and
#'   narrows the set.
#' @param marginals.as.percent Logical; `p1i` and `p2i` are percentages.
#' @param agreement.as.percent Logical; `agreement` is a percentage.
#' @return An object of class `"enum2x2"`: a list with `status`, `tables` (a
#'   data frame of `ai`, `bi`, `ci`, `di`, one row per compatible table), `ni`,
#'   `reason`, and `published` (the figures the call declared). Use
#'   [cell_ranges()] for the span of each cell, [n_tables()] for how many
#'   survived, and [unique_table()] where exactly one did.
#' @seealso [enum.2x2()] for the data frame form, and [metafor::conv.2x2()] for
#'   the point reconstruction this is the set-valued counterpart to.
#' @examples
#' # A pooled series of 768 patients, two readings of the DSM-5 delirium
#' # criteria, published at kappa = 0.29. The digits determine the table.
#' r <- recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29")
#' r
#' unique_table(r)
#'
#' # Coarser reporting on a larger sample leaves a set.
#' s <- recover2x2(20306, n1i = 866, n2i = 1603, kappa = "0.22")
#' n_tables(s)
#' cell_ranges(s)
#'
#' # Figures that do not cohere are reported, not signaled.
#' recover2x2(370, n1i = 165, n2i = 160, kappa = "0.48",
#'            agreement = "73", agreement.as.percent = TRUE)
#' @export
recover2x2 <- function(ni, n1i = NULL, n2i = NULL, p1i = NULL, p2i = NULL,
                       kappa = NULL, agreement = NULL,
                       marginals.as.percent = FALSE,
                       agreement.as.percent = FALSE) {
  if (is.logical(ni) || !.is_whole_number(ni) || ni <= 0) {
    stop_invalid_input(sprintf("ni must be a positive integer, got %s", .fmt(ni)))
  }
  if (ni * ni > .MAX_EXACT) {
    stop_invalid_input(sprintf(
      "ni = %s is too large to square exactly; this implementation requires ni^2 <= 2^52",
      .fmt(ni)))
  }
  declared <- list(p1i = p1i, p2i = p2i, kappa = kappa, agreement = agreement)
  for (nm in names(declared)) {
    v <- declared[[nm]]
    if (!is.null(v) && !(is.character(v) && length(v) == 1L && !is.na(v))) {
      stop_invalid_input(sprintf(
        paste0("%s must be the string the source printed, not %s; ",
               "a float cannot distinguish '0.10' from '0.1'"), nm, typeof(v)))
    }
  }
  if (is.null(kappa)) {
    return(.recovery(INSUFFICIENT, ni = ni, reason =
      "no closing statistic was reported; kappa is needed alongside the marginals"))
  }
  for (side in c("a", "b")) {
    count <- if (side == "a") n1i else n2i
    printed <- if (side == "a") p1i else p2i
    if (!is.null(count) && !is.null(printed)) {
      stop_invalid_input(sprintf("give %s or %s for the %s criterion, not both",
                                 .side_count[[side]], .side_printed[[side]],
                                 .side_word[[side]]))
    }
    if (is.null(count) && is.null(printed)) {
      return(.recovery(INSUFFICIENT, ni = ni, reason = sprintf(
        "neither %s nor %s was reported for the %s criterion",
        .side_count[[side]], .side_printed[[side]], .side_word[[side]])))
    }
  }
  parsed_kappa <- .validate_decimal("kappa", kappa)
  if (!is.null(agreement)) .validate_decimal("agreement", agreement)
  if (!is.null(p1i)) .validate_decimal("p1i", p1i)
  if (!is.null(p2i)) .validate_decimal("p2i", p2i)
  if (.dec_cmp_abs_one(parsed_kappa) > 0L) {
    stop_invalid_input(sprintf("kappa = %s is outside [-1, 1], the values it can take", kappa))
  }

  a_counts <- .marginal_counts(ni, n1i, p1i, marginals.as.percent, "a")
  b_counts <- .marginal_counts(ni, n2i, p2i, marginals.as.percent, "b")

  kiv <- exact_interval(kappa)
  giv <- if (is.null(agreement)) NULL else exact_interval(agreement, agreement.as.percent)
  tables <- .enumerate(ni, a_counts, b_counts, kiv, giv)

  published <- list(n1i = n1i, n2i = n2i, p1i = p1i, p2i = p2i,
                    kappa = kappa, agreement = agreement)
  published <- published[!vapply(published, is.null, logical(1))]
  # Without the flags an archived result cannot be replayed from what it records.
  if (!is.null(p1i) || !is.null(p2i)) {
    published$marginals.as.percent <- isTRUE(marginals.as.percent)
  }
  if (!is.null(agreement)) {
    published$agreement.as.percent <- isTRUE(agreement.as.percent)
  }

  if (nrow(tables) > 0L) {
    return(.recovery(if (nrow(tables) == 1L) UNIQUE else SET, tables = tables,
                     ni = ni, published = published))
  }

  # Nothing survived. Say which figure excluded, rather than blaming the kappa.
  if (!is.null(agreement)) {
    without <- recover2x2(ni, n1i = n1i, n2i = n2i, p1i = p1i, p2i = p2i,
                          kappa = kappa,
                          marginals.as.percent = marginals.as.percent)
    if (nrow(without$tables) > 0L) {
      return(.recovery(INFEASIBLE, ni = ni, published = published, reason = sprintf(
        "the marginals and kappa admit %d table(s); adding the published agreement admits none",
        nrow(without$tables))))
    }
  }
  span_lo <- Inf
  span_hi <- -Inf
  for (na in a_counts) {
    for (nb in b_counts) {
      lo <- tryCatch(kappa_min(na / ni, nb / ni),
                     enum2x2_undefined_statistic = function(e) NULL)
      if (is.null(lo)) next
      span_lo <- min(span_lo, lo)
      span_hi <- max(span_hi, kappa_max(na / ni, nb / ni))
    }
  }
  if (is.finite(span_lo)) {
    k <- rounding_interval(kappa)
    if (k[2L] < span_lo || k[1L] > span_hi) {
      return(.recovery(INFEASIBLE, ni = ni, published = published, reason = sprintf(
        "the published kappa lies outside [%.3f, %.3f], the range these marginals permit",
        span_lo, span_hi)))
    }
  }
  .recovery(INFEASIBLE, ni = ni, published = published,
            reason = "the marginals permit this kappa, but no integer table attains it")
}

.recovery <- function(status, tables = .empty_tables(), ni = NA_real_,
                      reason = NA_character_, published = list()) {
  rownames(tables) <- NULL
  structure(list(status = status, tables = tables, ni = ni,
                 reason = reason, published = published),
            class = "enum2x2")
}

#' Recover a series of published comparisons
#'
#' `recover_many` applies [recover2x2()] to a list of comparisons and never
#' signals: a row whose figures cannot describe any table comes back with status
#' `"impossible"` and the reason, so one malformed row does not stop the rest,
#' and a caller can count what was skipped. Counting those rows as
#' `"insufficient"` would inflate how many sources under-reported.
#'
#' @param rows A list of named lists or a data frame; names not among the
#'   arguments of [recover2x2()] are ignored, so a spreadsheet of published
#'   comparisons can carry a label column.
#' @return A list of `"enum2x2"` objects, one per row, in order.
#' @examples
#' out <- recover_many(list(
#'   list(study = "delirium", ni = 768, n1i = 158, n2i = 466, kappa = "0.29"),
#'   list(study = "impossible", ni = 100, n1i = 101, n2i = 50, kappa = "0.5"),
#'   list(study = "under-reported", ni = 768, n1i = 158, kappa = "0.29")))
#' vapply(out, function(r) r$status, character(1))
#' @export
recover_many <- function(rows) {
  if (is.data.frame(rows)) {
    rows <- lapply(seq_len(nrow(rows)), function(i) as.list(rows[i, , drop = FALSE]))
  }
  accepted <- setdiff(names(formals(recover2x2)), "")
  lapply(rows, function(row) {
    row <- as.list(row)
    row <- row[names(row) %in% accepted]
    row <- row[!vapply(row, function(v) length(v) == 0L ||
                         (length(v) == 1L && is.na(v)), logical(1))]
    tryCatch(do.call(recover2x2, row),
             enum2x2_invalid_input = function(e)
               .recovery(IMPOSSIBLE, ni = if (is.null(row$ni)) NA_real_ else row$ni,
                         reason = conditionMessage(e)),
             error = function(e)
               .recovery(INSUFFICIENT, ni = if (is.null(row$ni)) NA_real_ else row$ni,
                         reason = conditionMessage(e)))
  })
}

# --------------------------------------------------------------- the accessors

#' What a recovery determined
#'
#' `n_tables` is how many tables survived, `cell_ranges` the smallest and
#' largest value each cell takes across them, and `unique_table` the single
#' table where the figures determine exactly one. `unique_table` signals rather
#' than guessing where they do not.
#'
#' @param x An object of class `"enum2x2"` from [recover2x2()].
#' @param ... Unused, for compatibility with the generics.
#' @param row.names,optional Passed to the `as.data.frame` generic.
#' @return `n_tables` a single integer; `cell_ranges` a data frame with columns
#'   `cell`, `min` and `max`; `unique_table` a one-row data frame;
#'   `as.data.frame` the full set of compatible tables.
#' @examples
#' r <- recover2x2(20306, n1i = 866, n2i = 1603, kappa = "0.22")
#' n_tables(r)
#' cell_ranges(r)
#' head(as.data.frame(r))
#'
#' unique_table(recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29"))
#' @export
n_tables <- function(x) {
  stopifnot(inherits(x, "enum2x2"))
  nrow(x$tables)
}

#' @rdname n_tables
#' @export
cell_ranges <- function(x) {
  stopifnot(inherits(x, "enum2x2"))
  if (nrow(x$tables) == 0L) {
    stop_invalid_input("no tables survived, so no cell has a range")
  }
  data.frame(cell = c("ai", "bi", "ci", "di"),
             min = vapply(x$tables, min, numeric(1)),
             max = vapply(x$tables, max, numeric(1)),
             row.names = NULL)
}

#' @rdname n_tables
#' @export
unique_table <- function(x) {
  stopifnot(inherits(x, "enum2x2"))
  if (!identical(x$status, UNIQUE)) {
    stop_invalid_input(sprintf(
      paste0("the report is '%s', so it does not determine one table; ",
             "use as.data.frame() for all %d of them"), x$status, nrow(x$tables)))
  }
  x$tables
}

#' @rdname n_tables
#' @export
as.data.frame.enum2x2 <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$tables, row.names = row.names, optional = optional, ...)
}

#' @rdname n_tables
#' @export
print.enum2x2 <- function(x, ...) {
  headline <- switch(x$status,
                     unique = "unique",
                     set = "set-identified",
                     infeasible = "infeasible",
                     insufficient = "insufficient inputs",
                     impossible = "impossible figures",
                     x$status)
  cat(sprintf("enum2x2 recovery -- %s\n", headline))
  if (!is.na(x$ni)) cat(sprintf("  ni = %s\n", .fmt(x$ni)))
  n <- nrow(x$tables)
  if (n == 1L) {
    t <- x$tables
    cat(sprintf("  1 compatible table: ai = %s, bi = %s, ci = %s, di = %s\n",
                .fmt(t$ai), .fmt(t$bi), .fmt(t$ci), .fmt(t$di)))
  } else if (n > 1L) {
    r <- cell_ranges(x)
    num <- function(v) format(v, scientific = FALSE, trim = TRUE)
    cat(sprintf("  %d compatible tables\n", n))
    cat(paste0("  ", paste(sprintf("%s %s-%s", r$cell, num(r$min), num(r$max)),
                           collapse = ", "), "\n"))
  } else {
    cat("  no compatible table\n")
  }
  if (!is.na(x$reason)) cat(sprintf("  %s\n", x$reason))
  invisible(x)
}
