# A second implementation, written from the definitions.
#
# Agreement between two implementations that share their arithmetic establishes
# only that the arithmetic is consistent, so nothing below is taken from the
# package: the interval comes from counting the digits of the printed string,
# kappa is composed as (p_o - p_e) / (1 - p_e) out of rational operations
# defined here, the printed figures are formatted from that kappa, and
# membership is a cross-multiplication of those rationals.
#
# The sweeps run at n = 16, 20 and 24. Rationals are reduced to lowest terms
# after each composition, so a kappa denominator never exceeds n^2 and an
# interval denominator never exceeds twice the printed scale. The last test
# asserts those bounds rather than assuming them: they are what makes plain
# double arithmetic exact here, and what lets this check avoid the package's
# big-integer layer as well.

# --- exact rationals, vectorized, num over den with den > 0

rat <- function(num, den) {
  k <- max(length(num), length(den))
  num <- rep_len(num, k)
  den <- rep_len(den, k)
  s <- ifelse(den < 0, -1, 1)
  list(num = num * s, den = den * s)
}
rat_add <- function(x, y) rat(x$num * y$den + y$num * x$den, x$den * y$den)
rat_sub <- function(x, y) rat(x$num * y$den - y$num * x$den, x$den * y$den)
rat_mul <- function(x, y) rat(x$num * y$num, x$den * y$den)
rat_div <- function(x, y) rat(x$num * y$den, x$den * y$num)
rat_value <- function(x) x$num / x$den

# Modulo on integer-valued doubles, corrected so a rounded quotient cannot put
# the remainder outside [0, b).
vmod <- function(a, b) {
  r <- a - floor(a / b) * b
  r <- ifelse(r < 0, r + b, r)
  ifelse(r >= b, r - b, r)
}
rat_reduce <- function(x) {
  a <- abs(x$num)
  b <- abs(x$den)
  while (any(b > 0)) {
    t <- b
    b <- ifelse(b > 0, vmod(a, b), 0)
    a <- ifelse(t > 0, t, a)
  }
  g <- ifelse(a == 0, 1, a)
  rat(x$num / g, x$den / g)
}

# --- what a printed decimal string stands for, from its digits alone

ind_interval <- function(printed) {
  neg <- substr(printed, 1L, 1L) == "-"
  body <- sub("^-", "", printed)
  parts <- strsplit(body, ".", fixed = TRUE)[[1L]]
  frac <- if (length(parts) > 1L) parts[2L] else ""
  scale <- 10^nchar(frac)
  value <- as.numeric(paste0(parts[1L], frac))
  if (neg) value <- -value
  half <- rat(1, 2 * scale)
  list(lo = rat_reduce(rat_sub(rat(value, scale), half)),
       hi = rat_reduce(rat_add(rat(value, scale), half)))
}

# --- kappa, composed from the definition rather than from a rearrangement

ind_expected_agreement <- function(a, b, c, d) {
  n <- a + b + c + d
  one <- rat(1, 1)
  row <- rat(a + b, n)
  col <- rat(a + c, n)
  rat_reduce(rat_add(rat_mul(row, col), rat_mul(rat_sub(one, row), rat_sub(one, col))))
}

ind_kappa <- function(a, b, c, d) {
  n <- a + b + c + d
  one <- rat(1, 1)
  p_o <- rat(a + d, n)
  p_e <- ind_expected_agreement(a, b, c, d)
  rat_reduce(rat_div(rat_sub(p_o, p_e), rat_sub(one, p_e)))
}

rat_within <- function(x, iv) {
  (iv$lo$num * x$den <= x$num * iv$lo$den) & (x$num * iv$hi$den <= iv$hi$num * x$den)
}

reimplementation_agrees <- function(n, dp_marg, dp_kappa) {
  tab <- all_tables_of(n)
  ok <- non_degenerate(tab, n)
  row <- rat(tab$ai + tab$bi, n)
  col <- rat(tab$ai + tab$ci, n)
  p_e <- ind_expected_agreement(tab$ai, tab$bi, tab$ci, tab$di)
  kap <- ind_kappa(tab$ai, tab$bi, tab$ci, tab$di)
  # Tables where kappa does not exist, and the degenerate margins the sweep
  # leaves out, take no part in either set.
  eligible <- ok & (p_e$num != p_e$den)

  checked <- 0L
  mismatches <- character(0)
  for (i in which(ok)) {
    printed_a <- sprintf("%.*f", dp_marg, rat_value(rat(row$num[i], row$den[i])))
    printed_b <- sprintf("%.*f", dp_marg, rat_value(rat(col$num[i], col$den[i])))
    printed_k <- sprintf("%.*f", dp_kappa, rat_value(rat(kap$num[i], kap$den[i])))
    ia <- ind_interval(printed_a)
    ib <- ind_interval(printed_b)
    ik <- ind_interval(printed_k)
    keep <- eligible & rat_within(row, ia) & rat_within(col, ib) & rat_within(kap, ik)
    keep[is.na(keep)] <- FALSE
    mine <- as_key(tab[keep, , drop = FALSE])
    theirs <- as_key(recover2x2(n, p1i = printed_a, p2i = printed_b,
                                kappa = printed_k)$tables)
    self <- paste(tab$ai[i], tab$bi[i], tab$ci[i], tab$di[i], sep = "/")
    if (!identical(mine, theirs)) {
      mismatches <- c(mismatches, sprintf("%s: only here {%s}, only there {%s}",
                                          self, paste(setdiff(mine, theirs), collapse = ","),
                                          paste(setdiff(theirs, mine), collapse = ",")))
    } else if (!self %in% mine) {
      mismatches <- c(mismatches, paste0("absent from its own set: ", self))
    }
    checked <- checked + 1L
  }
  list(checked = checked, mismatches = mismatches)
}

test_that("a reimplementation from the definitions returns the same set at n = 16", {
  res <- reimplementation_agrees(16, 2, 1)
  expect_equal(res$checked, 905)
  expect_equal(head(res$mismatches, 3), character(0))
})

test_that("a reimplementation from the definitions returns the same set at n = 20", {
  res <- reimplementation_agrees(20, 2, 2)
  expect_gt(res$checked, 1000)
  expect_equal(head(res$mismatches, 3), character(0))
})

test_that("a reimplementation from the definitions returns the same set at n = 24", {
  res <- reimplementation_agrees(24, 3, 2)
  expect_gt(res$checked, 2000)
  expect_equal(head(res$mismatches, 3), character(0))
})

test_that("the reimplementation reproduces kappa where it is already known", {
  # 2/0/2/20 has kappa exactly 5/8; 158/0/308/302 is the delirium table.
  expect_equal(rat_value(ind_kappa(2, 0, 2, 20)), 0.625)
  expect_equal(rat_value(ind_kappa(158, 0, 308, 302)), kappa2x2(158, 0, 308, 302))
  # Reversing both labels leaves kappa unchanged, and so does swapping the two
  # criteria; both are properties of the definition, not of this code.
  expect_equal(rat_value(ind_kappa(800, 10, 90, 100)),
               rat_value(ind_kappa(100, 90, 10, 800)))
  expect_equal(rat_value(ind_kappa(100, 10, 90, 800)),
               rat_value(ind_kappa(100, 90, 10, 800)))
  # Perfect agreement is 1 and independence is 0, exactly.
  expect_equal(rat_value(ind_kappa(30, 0, 0, 70)), 1)
  expect_equal(rat_value(ind_kappa(9, 21, 21, 49)), 0)
})

test_that("the reimplementation's intervals are the ones the digits imply", {
  iv <- ind_interval("0.29")
  expect_equal(rat_value(iv$lo), 0.285)
  expect_equal(rat_value(iv$hi), 0.295)
  expect_equal(rat_value(ind_interval("0.290")$lo), 0.2895)
  expect_equal(rat_value(ind_interval("-0.20")$lo), -0.205)
  expect_equal(rat_value(ind_interval("-0.20")$hi), -0.195)
  expect_equal(rat_value(ind_interval("1")$lo), 0.5)
})

test_that("reducing to lowest terms is exact", {
  # The reduction is what keeps every later product inside a double, so it is
  # checked on its own: gcd against a known value, and a reduced fraction that
  # still names the value it started from.
  x <- rat_reduce(rat(2 * 3 * 5 * 7 * 11, 3 * 5 * 13))
  expect_equal(x$num, 2 * 7 * 11)
  expect_equal(x$den, 13)
  expect_equal(rat_reduce(rat(0, 7)), list(num = 0, den = 1))
  expect_equal(vmod(7962624, 331776), 0)
  expect_equal(vmod(2^52 + 3, 5), (2^52 + 3) %% 5)
  big <- rat_reduce(rat(7962624, 331776))
  expect_equal(rat_value(big), 24)
})

test_that("every quantity the reimplementation forms is exact in a double", {
  # The independence of this check rests on its own arithmetic being exact, so
  # the magnitudes are asserted rather than assumed. In lowest terms a kappa
  # denominator divides n^2 - chance, so it cannot exceed n^2.
  for (n in c(16, 20, 24)) {
    tab <- all_tables_of(n)
    kap <- ind_kappa(tab$ai, tab$bi, tab$ci, tab$di)
    finite <- is.finite(kap$num) & is.finite(kap$den) & kap$den != 0
    worst_num <- max(abs(kap$num[finite]))
    worst_den <- max(abs(kap$den[finite]))
    expect_lte(worst_den, n^2)
    expect_lte(worst_num, n^2)
    for (printed in c("0.1", "0.99", "0.001", "-0.20", "0.512")) {
      iv <- ind_interval(printed)
      expect_lt(worst_num * max(abs(c(iv$lo$den, iv$hi$den))), 2^53)
      expect_lt(worst_den * max(abs(c(iv$lo$num, iv$hi$num))), 2^53)
      expect_lt(n * max(abs(c(iv$lo$den, iv$hi$den))), 2^53)
    }
  }
})
