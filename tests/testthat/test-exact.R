# The arithmetic the membership test rests on.
#
# R has no rational type, so every comparison in this package is a
# cross-multiplication of integers carried in doubles, with a big-integer
# fallback where a product would leave the range a double holds exactly. If that
# layer is wrong, everything above it is wrong quietly, so it is checked on its
# own rather than only through the enumeration.

prod_le <- enum2x2:::.prod_le
prod_cmp_exact <- enum2x2:::.prod_cmp_exact
bi_from <- enum2x2:::.bi_from
bi_mul <- enum2x2:::.bi_mul
bi_cmp <- enum2x2:::.bi_cmp
least_c_above <- enum2x2:::.least_c_above
greatest_c_below <- enum2x2:::.greatest_c_below
MAX_EXACT <- enum2x2:::.MAX_EXACT

test_that("big-integer multiplication agrees with double multiplication where doubles are exact", {
  # Below 2^53 a double product is the true product, so it is an oracle for the
  # limb arithmetic. Many draws, no seed.
  for (i in 1:400) {
    a <- floor(runif(1, 0, 2^26))
    b <- floor(runif(1, 0, 2^26))
    expect_equal(bi_cmp(bi_mul(bi_from(a), bi_from(b)), bi_from(a * b)), 0L)
  }
  expect_equal(bi_cmp(bi_from(0), bi_from(0)), 0L)
  expect_equal(bi_cmp(bi_mul(bi_from(0), bi_from(12345)), bi_from(0)), 0L)
  expect_equal(bi_cmp(bi_from(1), bi_from(0)), 1L)
  expect_equal(bi_cmp(bi_from(9999999), bi_from(10000000)), -1L)
})

test_that("big-integer multiplication is commutative and distributes", {
  for (i in 1:200) {
    a <- floor(runif(1, 0, 2^50))
    b <- floor(runif(1, 0, 2^50))
    expect_equal(bi_cmp(bi_mul(bi_from(a), bi_from(b)),
                        bi_mul(bi_from(b), bi_from(a))), 0L)
    # (2a)*b == a*(2b), a relation the limb code cannot satisfy by accident.
    expect_equal(bi_cmp(bi_mul(bi_from(2 * a), bi_from(b)),
                        bi_mul(bi_from(a), bi_from(2 * b))), 0L)
  }
})

test_that("the exact product comparison orders products doubles cannot hold", {
  for (i in 1:400) {
    p <- floor(runif(1, 2^40, 2^50))
    q <- floor(runif(1, 2^40, 2^50))
    expect_gt(p * q, MAX_EXACT)
    expect_equal(prod_cmp_exact(p, q, p, q), 0L)
    expect_equal(prod_cmp_exact(p, q, p, q + 1), -1L)
    expect_equal(prod_cmp_exact(p, q + 1, p, q), 1L)
    expect_true(prod_le(p, q, p, q))
    expect_true(prod_le(p, q, p, q + 1))
    expect_false(prod_le(p, q + 1, p, q))
  }
})

test_that("the exact product comparison gets the signs right", {
  expect_equal(prod_cmp_exact(-2^40, 2^40, 2^40, 2^40), -1L)
  expect_equal(prod_cmp_exact(-2^40, -2^40, 2^40, 2^40), 0L)
  expect_equal(prod_cmp_exact(-2^40, 2^40, -2^40, 2^41), 1L)
  expect_equal(prod_cmp_exact(0, 2^50, -1, 2^50), 1L)
  expect_equal(prod_cmp_exact(0, 2^50, 0, -2^50), 0L)
  expect_true(prod_le(-1, 2^50, 0, 1))
  expect_false(prod_le(1, 2^50, 0, 1))
})

test_that("the fast and the exact path agree wherever both are defined", {
  for (i in 1:300) {
    a <- floor(runif(1, -2^25, 2^25))
    b <- floor(runif(1, -2^25, 2^25))
    c <- floor(runif(1, -2^25, 2^25))
    d <- floor(runif(1, -2^25, 2^25))
    expect_identical(prod_le(a, b, c, d), a * b <= c * d)
    expect_identical(prod_le(a, b, c, d), prod_cmp_exact(a, b, c, d) <= 0L)
  }
})

test_that("the product comparison is vectorized and recycles", {
  a <- c(3, -4, 0, 2^40)
  expect_identical(prod_le(a, 2, 6, 1), c(TRUE, TRUE, TRUE, FALSE))
  expect_identical(prod_le(2, a, 1, 6), c(TRUE, TRUE, TRUE, FALSE))
  expect_identical(prod_le(numeric(0), 1, 1, 1), logical(0))
})

test_that("the bounds search finds the same integers division would, without dividing", {
  for (i in 1:300) {
    a <- floor(runif(1, -1000, 1000))
    m <- floor(runif(1, 1, 500))
    b <- floor(runif(1, 1, 500))
    # Search bounds wide enough to contain the quotient, so the answer is the
    # quotient rather than a clamp.
    lo <- -1e6
    hi <- 1e6
    expect_equal(least_c_above(a, m, b, lo, hi), ceiling(a * m / b))
    expect_equal(greatest_c_below(a, m, b, lo, hi), floor(a * m / b))
  }
})

test_that("the bounds search reports an empty range rather than clamping into one", {
  # No integer in [0, 5] is at or above 100, and none is at or below -3.
  expect_equal(least_c_above(100, 1, 1, 0, 5), 6)
  expect_equal(greatest_c_below(-3, 1, 1, 0, 5), -1)
  expect_equal(least_c_above(0, 1, 1, 0, 5), 0)
  expect_equal(greatest_c_below(5, 1, 1, 0, 5), 5)
})

test_that("the enumeration really reaches the big-integer path", {
  # A kappa printed to twelve decimals makes the cross-products exceed what a
  # double holds, so the fallback is exercised rather than being dead code. The
  # answer must be the same table the two-decimal printing gives.
  k12 <- sprintf("%.12f", kappa2x2(158, 0, 308, 302))
  iv <- exact_interval(k12)
  expect_gt(abs(iv$hi_num) * 768^2, MAX_EXACT)
  r <- recover2x2(768, n1i = 158, n2i = 466, kappa = k12)
  expect_equal(r$status, "unique")
  expect_equal(as_key(r$tables), "158/0/308/302")
  expect_equal(as_key(recover2x2(768, n1i = 158, n2i = 466,
                                 kappa = "0.29")$tables), "158/0/308/302")
})

test_that("the interval of a printed figure is exact in integers", {
  expect_equal(exact_interval("0.29"), list(lo_num = 57, hi_num = 59, den = 200))
  expect_equal(exact_interval("0.290"), list(lo_num = 579, hi_num = 581, den = 2000))
  expect_equal(exact_interval("29"), list(lo_num = 57, hi_num = 59, den = 2))
  expect_equal(exact_interval("2.9e1"), list(lo_num = 57, hi_num = 59, den = 2))
  expect_equal(exact_interval("-0.29"), list(lo_num = -59, hi_num = -57, den = 200))
  expect_equal(exact_interval("4.26", as.percent = TRUE),
               list(lo_num = 851, hi_num = 853, den = 20000))
  expect_equal(exact_interval("0"), list(lo_num = -1, hi_num = 1, den = 2))
  expect_equal(exact_interval(".5"), list(lo_num = 9, hi_num = 11, den = 20))
})

test_that("the two routes to an interval in double precision agree", {
  # rounding_interval divides exact integers; the marginal gate builds a decimal
  # string instead, so that an exponent past the double range gives an infinite
  # bound rather than an error. Where both are defined they must be identical to
  # the last bit.
  bounds <- enum2x2:::.printed_bounds_double
  for (s in c("0.1", "0.10", "0.29", "0.512", "4.26", "66.4", "-0.20", "0",
              "1", "0.0", "1.0", "0.000001", "12.345", "0.9999")) {
    for (pct in c(FALSE, TRUE)) {
      expect_identical(bounds(s, pct), rounding_interval(s, pct), info = s)
    }
  }
  expect_equal(bounds("1e400", FALSE), c(Inf, Inf))
})
