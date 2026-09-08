# Does the answer depend on which rounding convention the source used?
#
# The interval is closed at both ends because a value falling exactly on a
# boundary rounds up under one convention and down under another, and published
# sources do not state which they used. The consequence is testable: re-run the
# enumeration under an endpoint rule that excludes the lower bound, one that
# excludes the upper, and one that excludes both, and the answer is unchanged
# except for tables whose statistic sits exactly on a bound. Those tables are
# the reason for the choice, so they are counted rather than assumed away.

# The same membership test as the package's, with the endpoint rule as an
# argument. Written here rather than exposed, because the package deliberately
# offers no way to open the interval.
compatible_under <- function(n, n1i, n2i, kappa, rule) {
  iv <- exact_interval(kappa)
  lo_ok <- function(num, den) {
    if (rule %in% c("closed", "open-upper")) num * iv$den >= iv$lo_num * den
    else num * iv$den > iv$lo_num * den
  }
  hi_ok <- function(num, den) {
    if (rule %in% c("closed", "open-lower")) num * iv$den <= iv$hi_num * den
    else num * iv$den < iv$hi_num * den
  }
  out <- character(0)
  for (ai in max(0, n1i + n2i - n):min(n1i, n2i)) {
    cells <- c(ai, n1i - ai, n2i - ai, n - n1i - n2i + ai)
    chance <- n1i * n2i + (n - n1i) * (n - n2i)
    den <- n * n - chance
    if (den == 0) next
    num <- (cells[1] + cells[4]) * n - chance
    if (lo_ok(num, den) && hi_ok(num, den)) {
      out <- c(out, paste(cells, collapse = "/"))
    }
  }
  sort(out)
}

on_a_bound <- function(n, n1i, n2i, kappa) {
  iv <- exact_interval(kappa)
  chance <- n1i * n2i + (n - n1i) * (n - n2i)
  den <- n * n - chance
  if (den == 0) return(FALSE)
  any(vapply(max(0, n1i + n2i - n):min(n1i, n2i), function(ai) {
    num <- (ai + (n - n1i - n2i + ai)) * n - chance
    num * iv$den == iv$lo_num * den || num * iv$den == iv$hi_num * den
  }, logical(1)))
}

test_that("the four endpoint rules agree except where a statistic sits on a bound", {
  n <- 24
  differed <- 0L
  differed_off_bound <- character(0)
  cases <- 0L
  # A grid across the marginals, with 2 and 4 included because 2/0/2/20 is the
  # case whose kappa, 5/8, lands exactly on a printed bound.
  grid <- c(2, 4, 9, 16, 22)
  for (n1i in grid) {
    for (n2i in grid) {
      chance <- n1i * n2i + (n - n1i) * (n - n2i)
      if (n * n - chance == 0) next
      for (k in sprintf("%.2f", seq(-0.5, 0.95, by = 0.01))) {
        closed <- compatible_under(n, n1i, n2i, k, "closed")
        others <- list(compatible_under(n, n1i, n2i, k, "open-lower"),
                       compatible_under(n, n1i, n2i, k, "open-upper"),
                       compatible_under(n, n1i, n2i, k, "open"))
        expect_identical(closed, as_key(recover2x2(n, n1i = n1i, n2i = n2i,
                                                   kappa = k)$tables))
        cases <- cases + 1L
        same <- all(vapply(others, identical, logical(1), closed))
        if (!same) {
          differed <- differed + 1L
          if (!on_a_bound(n, n1i, n2i, k)) {
            differed_off_bound <- c(differed_off_bound,
                                    sprintf("%d/%d %s", n1i, n2i, k))
          }
        }
      }
    }
  }
  expect_gt(cases, 3000)
  # Where the rules disagree, a statistic is exactly on a bound. Nowhere else.
  expect_equal(head(differed_off_bound, 5), character(0))
  # And the disagreement is real: a half-open rule would change the answer.
  expect_gt(differed, 0)
})

test_that("a table on a bound survives under the closed rule and is lost otherwise", {
  # Kappa for 2/0/2/20 is exactly 5/8, the upper bound of "0.62" and the lower
  # bound of "0.63". Under a half-open rule one of the two printings loses it.
  expect_equal(compatible_under(24, 2, 4, "0.62", "closed"), "2/0/2/20")
  expect_equal(compatible_under(24, 2, 4, "0.63", "closed"), "2/0/2/20")
  expect_equal(compatible_under(24, 2, 4, "0.62", "open-upper"), character(0))
  expect_equal(compatible_under(24, 2, 4, "0.63", "open-lower"), character(0))
  expect_true(on_a_bound(24, 2, 4, "0.62"))
  expect_false(on_a_bound(24, 2, 4, "0.61"))
})

test_that("the closed rule can only widen the set, never narrow it", {
  n <- 20
  for (n1i in c(3, 7, 11)) {
    for (n2i in c(4, 9, 15)) {
      for (k in c("0.1", "0.30", "-0.20", "0.55")) {
        closed <- compatible_under(n, n1i, n2i, k, "closed")
        for (rule in c("open-lower", "open-upper", "open")) {
          expect_true(all(compatible_under(n, n1i, n2i, k, rule) %in% closed))
        }
      }
    }
  }
})
