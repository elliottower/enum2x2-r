# Cohen's kappa and the bounds the marginals impose on it.
#
# Cell names follow metafor::conv.2x2 throughout: ai is the count positive on
# both criteria, bi positive on the first and negative on the second, ci
# negative on the first and positive on the second, di negative on both. The
# marginal totals are n1i = ai + bi and n2i = ai + ci, and ni = ai+bi+ci+di.

# Euclid, with the remainder corrected for a quotient that rounded, so the
# result is exact for any integer-valued doubles inside 2^53.
.gcd <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  while (b > 0) {
    r <- a - floor(a / b) * b
    if (r < 0) r <- r + b
    if (r >= b) r <- r - b
    a <- b
    b <- r
  }
  a
}

#' Expected agreement between two binary criteria
#'
#' The agreement two criteria would reach by chance alone at the given positive
#' rates: `p1 * p2 + (1 - p1) * (1 - p2)`. This is the `p_e` of Cohen's kappa.
#'
#' @param p1,p2 Positive rates of the first and second criterion, in \[0, 1\].
#' @return A numeric vector.
#' @examples
#' expected_agreement(0.3, 0.3)
#' expected_agreement(0.0426, 0.0789)
#' @export
expected_agreement <- function(p1, p2) {
  p1 * p2 + (1 - p1) * (1 - p2)
}

#' Cohen's kappa, raw agreement and directional asymmetry of a 2x2 table
#'
#' `kappa2x2` returns unweighted Cohen's kappa in double precision.
#' `exact_kappa` returns the same quantity as a ratio of two integers in lowest
#' terms, which is what membership in a rounding interval is decided on.
#' `agreement2x2` is the observed proportion of agreement, `(ai + di) / ni`.
#' `asymmetry2x2` is `abs(bi - ci) / (bi + ci)`: zero when the disagreement is
#' evenly split between the two criteria and one when it runs entirely one way.
#'
#' Kappa is undefined where expected agreement is exactly one, which happens
#' when both criteria are positive for everyone or negative for everyone; those
#' calls signal `enum2x2_undefined_statistic` rather than returning `NaN`.
#'
#' @param ai Count positive on both criteria.
#' @param bi Count positive on the first criterion and negative on the second.
#' @param ci Count negative on the first criterion and positive on the second.
#' @param di Count negative on both criteria.
#' @return `kappa2x2`, `agreement2x2` and `asymmetry2x2` are vectorized over
#'   their arguments and return a numeric vector; `asymmetry2x2` gives
#'   `NA_real_` where the criteria never disagree. `exact_kappa` takes one table
#'   and returns a list with integer-valued `num` and `den` in lowest terms.
#' @examples
#' # The pooled delirium series: two readings of the DSM-5 criteria on 768
#' # patients, one reading's positives wholly inside the other's.
#' kappa2x2(158, 0, 308, 302)
#' agreement2x2(158, 0, 308, 302)
#' asymmetry2x2(0, 308)
#' exact_kappa(158, 0, 308, 302)
#' @export
kappa2x2 <- function(ai, bi, ci, di) {
  n <- ai + bi + ci + di
  if (any(n <= 0)) stop_invalid_input("the table is empty")
  p_o <- (ai + di) / n
  p_e <- expected_agreement((ai + bi) / n, (ai + ci) / n)
  if (any(abs(1 - p_e) < 1e-12)) {
    stop_undefined_statistic("expected agreement is 1, so kappa is undefined")
  }
  (p_o - p_e) / (1 - p_e)
}

#' @rdname kappa2x2
#' @export
exact_kappa <- function(ai, bi, ci, di) {
  n <- ai + bi + ci + di
  if (n <= 0) stop_invalid_input("the table is empty")
  if (n * n > .MAX_EXACT) {
    stop_invalid_input("the sample is too large to square exactly (ni^2 must fit in 2^52)")
  }
  chance <- (ai + bi) * (ai + ci) + (ci + di) * (bi + di)
  den <- n * n - chance
  if (den == 0) {
    stop_undefined_statistic("expected agreement is 1, so kappa is undefined")
  }
  num <- (ai + di) * n - chance
  g <- .gcd(num, den)
  if (g == 0) g <- 1
  list(num = num / g, den = den / g)
}

#' @rdname kappa2x2
#' @export
agreement2x2 <- function(ai, bi, ci, di) {
  n <- ai + bi + ci + di
  if (any(n <= 0)) stop_invalid_input("the table is empty")
  (ai + di) / n
}

#' @rdname kappa2x2
#' @export
asymmetry2x2 <- function(bi, ci) {
  d <- bi + ci
  ifelse(d == 0, NA_real_, abs(bi - ci) / d)
}

#' The range of kappa two marginals permit
#'
#' Positive rates bound Cohen's kappa from both sides. `kappa_max` uses the
#' largest observed agreement the marginals allow, `1 - abs(p1 - p2)`, and
#' `kappa_min` the smallest. A published kappa outside `[kappa_min, kappa_max]`
#' cannot have come from a table with these marginals, which is a more specific
#' finding than reporting that no table survived.
#'
#' Both are vectorized, and both signal `enum2x2_undefined_statistic` where
#' expected agreement is one.
#'
#' @param p1,p2 Positive rates of the first and second criterion, in \[0, 1\].
#' @return A numeric vector.
#' @examples
#' # Equal marginals permit kappa = 1; unequal marginals do not.
#' kappa_max(0.3, 0.3)
#' kappa_max(158 / 768, 466 / 768)
#' kappa_min(158 / 768, 466 / 768)
#' @export
kappa_max <- function(p1, p2) {
  p_e <- expected_agreement(p1, p2)
  p_o_max <- pmin(p1, p2) + pmin(1 - p1, 1 - p2)
  if (any(abs(1 - p_e) < 1e-12)) {
    stop_undefined_statistic("expected agreement is 1, so kappa_max is undefined")
  }
  (p_o_max - p_e) / (1 - p_e)
}

#' @rdname kappa_max
#' @export
kappa_min <- function(p1, p2) {
  p_e <- expected_agreement(p1, p2)
  p_o_min <- pmax(0, 1 - p1 - p2) + pmax(0, p1 + p2 - 1)
  if (any(abs(1 - p_e) < 1e-12)) {
    stop_undefined_statistic("expected agreement is 1, so kappa_min is undefined")
  }
  (p_o_min - p_e) / (1 - p_e)
}
