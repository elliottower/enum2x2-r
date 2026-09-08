# Printed figures, read as intervals.
#
# A source prints "0.29", not the value 0.29. The digits fix an interval, and
# the interval is what the enumeration filters on, so the string is parsed into
# an exact coefficient and exponent and never into a double.

# Parse a decimal string the way Python's decimal.Decimal does: into a sign, a
# coefficient of digits, and an exponent, so that value = sign * digits * 10^exp
# and the printed precision is recoverable from the exponent alone. "0.10" and
# "0.1" therefore differ here, which is the point.
.dec_parse <- function(s) {
  s <- trimws(s)
  s <- gsub("(?<=[0-9])_(?=[0-9])", "", s, perl = TRUE)
  if (grepl("^[+-]?(inf|infinity)$", s, ignore.case = TRUE) ||
      grepl("^[+-]?s?nan[0-9]*$", s, ignore.case = TRUE)) {
    return(list(ok = TRUE, finite = FALSE))
  }
  m <- regexec("^([+-]?)([0-9]*)(\\.([0-9]*))?([eE]([+-]?[0-9]+))?$", s)
  g <- regmatches(s, m)[[1L]]
  if (length(g) == 0L) return(list(ok = FALSE, finite = NA))
  int_part <- g[3L]
  frac_part <- g[5L]
  if (nchar(int_part) + nchar(frac_part) == 0L) return(list(ok = FALSE, finite = NA))
  exp_part <- if (nchar(g[7L])) as.numeric(g[7L]) else 0
  digits <- sub("^0+", "", paste0(int_part, frac_part))
  list(ok = TRUE, finite = TRUE,
       sign = if (identical(g[2L], "-")) -1 else 1,
       digits = digits,
       exp = exp_part - nchar(frac_part))
}

# Compare the magnitude of a parsed decimal with 1, exactly and without forming
# the value, so that a coefficient or exponent far beyond what a double holds
# still gives the right answer.
.dec_cmp_abs_one <- function(p) {
  if (!nchar(p$digits)) return(-1L)
  len <- nchar(p$digits)
  ord <- len + p$exp
  if (ord <= 0) return(-1L)
  if (ord >= 2) return(1L)
  if (identical(p$digits, paste0("1", strrep("0", len - 1L)))) 0L else 1L
}

.dec_as_double <- function(p) {
  if (!nchar(p$digits)) return(0)
  p$sign * as.numeric(paste0(p$digits, "e", format(p$exp, scientific = FALSE)))
}

#' The interval a printed figure stands for, in exact integers
#'
#' A published figure is a decimal string, so it names a closed interval of
#' values rather than a value. `exact_interval` returns that interval as two
#' integer numerators over a common integer denominator, which is what
#' membership is decided on: a table's statistic is also a ratio of integers, so
#' the comparison is exact and there is no tolerance to choose.
#'
#' The interval is closed at both ends deliberately. A value falling exactly on
#' a boundary rounds up under one convention and down under another, and
#' published sources do not state which they used. Admitting both ends can only
#' widen the candidate set; excluding one end can drop the true table.
#'
#' All three integers are carried in doubles and must not exceed `2^52`, which
#' allows any figure with at most 15 significant digits printed to at most 13
#' decimal places (11 when `as.percent` is `TRUE`).
#'
#' @param printed The literal string the source printed, such as `"0.29"`.
#'   `"0.10"` and `"0.1"` are different inputs.
#' @param as.percent Logical; if `TRUE`, `printed` is a percentage and the
#'   interval is scaled to a proportion.
#' @return A list with elements `lo_num`, `hi_num` and `den`; the interval is
#'   `[lo_num / den, hi_num / den]`.
#' @seealso [rounding_interval()] for the same interval in double precision.
#' @examples
#' exact_interval("0.29")
#' exact_interval("0.290")
#' exact_interval("4.26", as.percent = TRUE)
#' @export
exact_interval <- function(printed, as.percent = FALSE) {
  p <- .dec_parse(printed)
  if (!isTRUE(p$ok)) {
    stop_invalid_input(sprintf("%s is not a decimal number", .quote(printed)))
  }
  if (!isTRUE(p$finite)) {
    stop_invalid_input(sprintf("%s is not a finite number", .quote(printed)))
  }
  if (nchar(p$digits) > 15L) {
    stop_invalid_input(sprintf(
      "%s carries more than 15 significant digits, beyond what this implementation holds exactly",
      .quote(printed)))
  }
  m <- p$sign * (if (nchar(p$digits)) as.numeric(p$digits) else 0)
  e <- p$exp
  if (e >= 0) {
    pw <- 10^e
    lo_num <- (2 * m - 1) * pw
    hi_num <- (2 * m + 1) * pw
    den <- 2
  } else {
    lo_num <- 2 * m - 1
    hi_num <- 2 * m + 1
    den <- 2 * 10^(-e)
  }
  if (isTRUE(as.percent)) den <- den * 100
  worst <- max(abs(lo_num), abs(hi_num), den)
  if (!is.finite(worst) || worst > .MAX_EXACT) {
    stop_invalid_input(sprintf(
      "%s needs more precision than this implementation holds exactly (the interval bounds must fit in 2^52)",
      .quote(printed)))
  }
  list(lo_num = lo_num, hi_num = hi_num, den = den)
}

#' The interval a printed figure stands for, in double precision
#'
#' The closed interval of exact values that round to a printed string, as two
#' doubles. `printed` is the literal text of the source, so `"0.10"` and `"0.1"`
#' give different intervals, which is why figures are declared as strings.
#'
#' This is the readable form of [exact_interval()]. Membership is decided on the
#' exact form; use this one to see what an interval is.
#'
#' @inheritParams exact_interval
#' @return A numeric vector of length two, the lower and upper bound.
#' @examples
#' rounding_interval("0.1")
#' rounding_interval("0.10")
#' rounding_interval("4.26", as.percent = TRUE)
#' @export
rounding_interval <- function(printed, as.percent = FALSE) {
  iv <- exact_interval(printed, as.percent)
  c(iv$lo_num / iv$den, iv$hi_num / iv$den)
}

#' Does a value round to a printed string?
#'
#' `rounds_to` compares a double against the interval of [rounding_interval()],
#' with a slack of `1e-12` at each end to absorb representation error in the
#' value. `exactly_rounds_to` takes the value as a ratio of two integers and
#' decides membership by integer cross-multiplication, with no slack and no
#' tolerance parameter; this is what the enumeration uses.
#'
#' @param value A numeric value.
#' @param num,den Integer-valued numerators and denominators of the value, with
#'   `den` positive.
#' @inheritParams exact_interval
#' @return A logical vector.
#' @examples
#' rounds_to(0.105, "0.10")
#' rounds_to(0.105, "0.11")
#' rounds_to(0.12, "0.10")
#'
#' # Cohen's kappa for the table (2, 0, 2, 20) is exactly 5/8, which sits on the
#' # upper endpoint of "0.62" and the lower endpoint of "0.63".
#' k <- exact_kappa(2, 0, 2, 20)
#' exactly_rounds_to(k$num, k$den, "0.62")
#' exactly_rounds_to(k$num, k$den, "0.63")
#' @export
rounds_to <- function(value, printed, as.percent = FALSE) {
  b <- rounding_interval(printed, as.percent)
  b[1L] - 1e-12 <= value & value <= b[2L] + 1e-12
}

#' @rdname rounds_to
#' @export
exactly_rounds_to <- function(num, den, printed, as.percent = FALSE) {
  if (any(den <= 0)) stop_invalid_input("den must be positive")
  iv <- exact_interval(printed, as.percent)
  .prod_le(iv$lo_num, den, iv$den, num) & .prod_le(iv$den, num, iv$hi_num, den)
}

#' Every count on a sample whose proportion rounds to a printed string
#'
#' A marginal printed as a proportion or a percentage names an interval, and on
#' a sample of `n` that interval admits a contiguous run of integer counts. The
#' bounds are found by searching on exact integer comparisons rather than by
#' dividing, so nothing is lost to truncation at either end.
#'
#' @param printed The literal string the source printed.
#' @param n Sample size, a positive whole number.
#' @inheritParams exact_interval
#' @return An increasing numeric vector of counts, possibly of length zero.
#' @examples
#' counts_rounding_to("0.14", 200)
#'
#' # At n = 240 the marginal "0.512" admits exactly one count, and it is the one
#' # a truncating upper bound would lose.
#' counts_rounding_to("0.512", 240)
#' counts_rounding_to("4.26", 20306, as.percent = TRUE)
#' @export
counts_rounding_to <- function(printed, n, as.percent = FALSE) {
  if (!.is_whole_number(n) || n <= 0) {
    stop_invalid_input(sprintf("n must be a positive whole number, got %s", format(n)))
  }
  iv <- exact_interval(printed, as.percent)
  first <- .least_c_above(iv$lo_num, n, iv$den, 0, n + 1)
  last <- .greatest_c_below(iv$hi_num, n, iv$den, -1, n)
  if (last >= first) seq.int(first, last) else numeric(0)
}

# The same bounds as rounding_interval, reached through a decimal string rather
# than through exact_interval, so that a printed exponent past the range a
# double covers gives an infinite bound instead of an error. Whether a printed
# marginal can be a proportion at all is decided on these, so "1e400" has to
# reach the answer "outside [0, 1]" rather than a complaint about precision.
# Both routes convert the same exact rational and are correctly rounded, so for
# every representable figure they return identical doubles.
.printed_bounds_double <- function(printed, as.percent = FALSE) {
  p <- .dec_parse(printed)
  if (nchar(p$digits) > 15L) return(rounding_interval(printed, as.percent))
  m <- if (nchar(p$digits)) p$sign * as.numeric(p$digits) else 0
  e <- format(p$exp - (if (isTRUE(as.percent)) 2 else 0),
              scientific = FALSE, trim = TRUE)
  c(as.numeric(paste0(sprintf("%.1f", m - 0.5), "e", e)),
    as.numeric(paste0(sprintf("%.1f", m + 0.5), "e", e)))
}

.quote <- function(x) {
  if (is.character(x)) sprintf("'%s'", x) else format(x)
}
