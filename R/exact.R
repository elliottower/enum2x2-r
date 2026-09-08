# Exact integer arithmetic, without a dependency.
#
# Membership in a rounding interval is a comparison between two rationals: a
# decimal string the source printed, and a ratio of integers a table produces.
# Deciding it exactly means comparing p/q against r/s by cross-multiplication,
# and the products can leave the range a double holds exactly. R has no rational
# type and no integer wider than 32 bits, so this file supplies what the
# comparison needs.
#
# Every quantity here is an integer carried in a double, which is exact up to
# 2^53. Only the cross-products can exceed that, so only they are special-cased:
# .prod_le compares a*b with c*d in double arithmetic where the result is
# provably exact, and in base-10^7 limb arithmetic where it is not. The gmp
# package would do the same job through libgmp, at the cost of a hard
# dependency on a compiled package; the limb arithmetic below is thirty lines
# and runs only on the inputs that need it.

# Largest magnitude at which a double product is taken as exact. If the true
# product of two integers exceeds 2^53 then the computed product exceeds 2^52,
# so a computed product at or below 2^52 came from a true product below 2^53,
# which every double represents exactly.
.MAX_EXACT <- 2^52

.BI_BASE <- 1e7

# ---------------------------------------------------------------- big integers

# A non-negative integer as a vector of base-10^7 limbs, least significant
# first. Limb products stay below 10^14 and limb sums below 2^53, so every
# intermediate is exact.

.bi_from <- function(x) {
  if (x < 0) stop("negative value passed to .bi_from")
  if (x == 0) return(0)
  limbs <- numeric(0L)
  while (x > 0) {
    r <- x %% .BI_BASE
    limbs <- c(limbs, r)
    x <- (x - r) / .BI_BASE
  }
  limbs
}

.bi_trim <- function(a) {
  while (length(a) > 1L && a[length(a)] == 0) a <- a[-length(a)]
  a
}

.bi_mul <- function(a, b) {
  out <- numeric(length(a) + length(b))
  for (i in seq_along(a)) {
    if (a[i] == 0) next
    carry <- 0
    for (j in seq_along(b)) {
      t <- out[i + j - 1L] + a[i] * b[j] + carry
      r <- t %% .BI_BASE
      carry <- (t - r) / .BI_BASE
      out[i + j - 1L] <- r
    }
    k <- i + length(b)
    while (carry > 0) {
      t <- out[k] + carry
      r <- t %% .BI_BASE
      carry <- (t - r) / .BI_BASE
      out[k] <- r
      k <- k + 1L
    }
  }
  .bi_trim(out)
}

.bi_cmp <- function(a, b) {
  a <- .bi_trim(a)
  b <- .bi_trim(b)
  if (length(a) != length(b)) return(if (length(a) > length(b)) 1L else -1L)
  for (i in rev(seq_along(a))) {
    if (a[i] != b[i]) return(if (a[i] > b[i]) 1L else -1L)
  }
  0L
}

# ------------------------------------------------------------- cross-products

# Sign of a*b - c*d, exactly, for one set of integer-valued doubles.
.prod_cmp_exact <- function(a, b, c, d) {
  sab <- sign(a) * sign(b)
  scd <- sign(c) * sign(d)
  if (sab != scd) return(if (sab > scd) 1L else -1L)
  if (sab == 0) return(0L)
  m <- .bi_cmp(.bi_mul(.bi_from(abs(a)), .bi_from(abs(b))),
               .bi_mul(.bi_from(abs(c)), .bi_from(abs(d))))
  if (sab > 0) m else -m
}

#' Exact comparison of two integer products
#'
#' Tests `a * b <= c * d` for integer-valued doubles, elementwise, without
#' relying on the products being representable. Where both products are small
#' enough that double arithmetic is exact the comparison is made directly;
#' otherwise it falls back to big-integer arithmetic.
#'
#' @param a,b,c,d integer-valued numeric vectors, recycled to a common length.
#' @return A logical vector.
#' @noRd
.prod_le <- function(a, b, c, d) {
  lens <- c(length(a), length(b), length(c), length(d))
  if (min(lens) == 0L) return(logical(0L))
  n <- max(lens)
  a <- rep_len(a, n); b <- rep_len(b, n)
  c <- rep_len(c, n); d <- rep_len(d, n)
  ab <- a * b
  cd <- c * d
  ok <- abs(ab) <= .MAX_EXACT & abs(cd) <= .MAX_EXACT
  res <- logical(n)
  res[ok] <- ab[ok] <= cd[ok]
  slow <- which(!ok)
  for (i in slow) res[i] <- .prod_cmp_exact(a[i], b[i], c[i], d[i]) <= 0L
  res
}

# ------------------------------------------------------------------- division

# Smallest integer c in [lo, hi] with a * m <= c * b, or hi + 1 if there is
# none. Written as a search on the exact product comparison rather than as a
# division, so no quotient is ever formed in floating point.
.least_c_above <- function(a, m, b, lo, hi) {
  if (.prod_le(a, m, hi, b)) {
    while (lo < hi) {
      mid <- lo + floor((hi - lo) / 2)
      if (.prod_le(a, m, mid, b)) hi <- mid else lo <- mid + 1
    }
    lo
  } else {
    hi + 1
  }
}

# Largest integer c in [lo, hi] with c * b <= a * m, or lo - 1 if there is none.
.greatest_c_below <- function(a, m, b, lo, hi) {
  if (.prod_le(lo, b, a, m)) {
    while (lo < hi) {
      mid <- hi - floor((hi - lo) / 2)
      if (.prod_le(mid, b, a, m)) lo <- mid else hi <- mid - 1
    }
    lo
  } else {
    lo - 1
  }
}

# ---------------------------------------------------------------- input types

.is_whole_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == floor(x) && abs(x) <= .MAX_EXACT
}
