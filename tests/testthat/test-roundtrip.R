# Round-trip: print a known table, then recover from what was printed.
#
# For every table in a sweep, compute its statistics, round them to a stated
# precision, feed the strings back in, and require the table itself to be among
# those returned. A recovery that lost the table it was generated from would be
# unsound whatever else it got right, and this is the property that has to hold
# for every table rather than for a chosen one.

round_trip <- function(tab, n, dp_marg, dp_kappa, agreement_dp = NULL) {
  missing_from_own_set <- character(0)
  never_empty <- TRUE
  for (i in seq_len(nrow(tab))) {
    cells <- as.numeric(tab[i, c("ai", "bi", "ci", "di")])
    row <- printed_row(as.list(cells), n, dp_marg, dp_kappa)
    args <- list(n, p1i = row$p1i, p2i = row$p2i, kappa = row$kappa)
    if (!is.null(agreement_dp)) {
      args$agreement <- sprintf("%.*f", agreement_dp, (cells[1] + cells[4]) / n)
    }
    r <- do.call(recover2x2, args)
    if (n_tables(r) == 0L) never_empty <- FALSE
    self <- paste(cells, collapse = "/")
    if (!self %in% as_key(r$tables)) {
      missing_from_own_set <- c(missing_from_own_set, self)
    }
  }
  list(missing = missing_from_own_set, never_empty = never_empty)
}

test_that("every table at n = 17 is in the set recovered from its own printed figures", {
  n <- 17
  tab <- all_tables_of(n)
  tab <- tab[non_degenerate(tab, n), , drop = FALSE]
  res <- round_trip(tab, n, 4, 2)
  expect_gt(nrow(tab), 1000)
  expect_equal(head(res$missing, 5), character(0))
  expect_true(res$never_empty)
})

test_that("a published agreement alongside kappa never excludes the true table", {
  n <- 19
  tab <- all_tables_of(n)
  tab <- tab[non_degenerate(tab, n), , drop = FALSE]
  res <- round_trip(tab, n, 3, 2, agreement_dp = 3)
  expect_equal(head(res$missing, 5), character(0))
  expect_true(res$never_empty)
})

test_that("the round trip holds at printed precisions from coarse to fine", {
  n <- 14
  tab <- all_tables_of(n)
  tab <- tab[non_degenerate(tab, n), , drop = FALSE]
  for (dp_marg in 2:4) {
    for (dp_kappa in 1:3) {
      res <- round_trip(tab, n, dp_marg, dp_kappa)
      expect_equal(head(res$missing, 3), character(0))
    }
  }
})

test_that("the round trip holds at the sample sizes real reports use", {
  # Sizes and shapes drawn without a seed, so the property is asserted across
  # every draw rather than across one arrangement of them.
  draws <- 250
  missing_from_own_set <- character(0)
  widths <- integer(0)
  for (i in seq_len(draws)) {
    n <- sample(c(120L, 370L, 768L, 2000L, 7328L, 20306L), 1L)
    n1 <- sample.int(n - 1L, 1L)
    n2 <- sample.int(n - 1L, 1L)
    lo <- max(0L, n1 + n2 - n)
    hi <- min(n1, n2)
    ai <- lo + sample.int(hi - lo + 1L, 1L) - 1L
    cells <- c(ai, n1 - ai, n2 - ai, n - n1 - n2 + ai)
    p_e <- (n1 / n) * (n2 / n) + (1 - n1 / n) * (1 - n2 / n)
    if (abs(1 - p_e) < 1e-12) next
    dp <- sample(2:3, 1L)
    k <- sprintf("%.*f", dp, kappa2x2(cells[1], cells[2], cells[3], cells[4]))
    r <- recover2x2(n, n1i = n1, n2i = n2, kappa = k)
    self <- paste(cells, collapse = "/")
    if (!self %in% as_key(r$tables)) {
      missing_from_own_set <- c(missing_from_own_set, paste(self, k))
    }
    widths <- c(widths, n_tables(r))
  }
  expect_equal(head(missing_from_own_set, 5), character(0))
  expect_gt(length(widths), 200)
  expect_true(all(widths >= 1))
})

test_that("finer printed digits never widen the set and coarser never narrow it", {
  n <- 200
  cells <- c(40, 15, 10, 135)
  sizes <- vapply(2:6, function(dp) {
    row <- printed_row(as.list(cells), n, dp, dp)
    n_tables(recover2x2(n, p1i = row$p1i, p2i = row$p2i, kappa = row$kappa))
  }, numeric(1))
  expect_true(all(diff(sizes) <= 0))
  expect_gt(sizes[1], sizes[length(sizes)])
  expect_equal(sizes[length(sizes)], 1)
})

test_that("an extra published statistic can only narrow the set", {
  n <- 200
  cells <- c(40, 15, 10, 135)
  row <- printed_row(as.list(cells), n, 2, 1)
  without <- recover2x2(n, p1i = row$p1i, p2i = row$p2i, kappa = row$kappa)
  with_po <- recover2x2(n, p1i = row$p1i, p2i = row$p2i, kappa = row$kappa,
                        agreement = sprintf("%.4f", (cells[1] + cells[4]) / n))
  expect_lte(n_tables(with_po), n_tables(without))
  expect_gte(n_tables(with_po), 1)
  expect_true(all(as_key(with_po$tables) %in% as_key(without$tables)))
})

test_that("exact counts identify what printed percentages leave open", {
  cells <- c(125, 28, 49, 7126)
  n <- sum(cells)
  k <- sprintf("%.3f", kappa2x2(cells[1], cells[2], cells[3], cells[4]))
  loose <- recover2x2(n, p1i = "2.1", p2i = "2.4", kappa = k,
                      marginals.as.percent = TRUE)
  tight <- recover2x2(n, n1i = cells[1] + cells[2], n2i = cells[1] + cells[3],
                      kappa = k)
  expect_true(paste(cells, collapse = "/") %in% as_key(loose$tables))
  expect_equal(as_key(tight$tables), paste(cells, collapse = "/"))
  expect_lt(n_tables(tight), n_tables(loose))
})
