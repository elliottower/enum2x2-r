# Soundness and completeness over a finite state space.
#
# Brute force every integer table of size n, keep those whose printed summaries
# match the row's, and require that set to equal what the enumerator returns.
# Set equality, not membership: a check that only asks whether the true table
# survived cannot catch a rival being wrongly excluded, or a stranger wrongly
# kept.
#
# n = 18 and n = 22 are the sizes the Python package sweeps, at the same printed
# precisions: 1,330 and 2,300 tables, of which 1,258 and 2,212 have a kappa.

sweep_is_exactly_the_compatible_set <- function(n, dp_marg, dp_kappa) {
  tab <- all_tables_of(n)
  ok <- non_degenerate(tab, n)
  kap <- rep(NA_real_, nrow(tab))
  kap[ok] <- kappa_of(tab[ok, , drop = FALSE], n)
  na <- (tab$ai + tab$bi) / n
  nb <- (tab$ai + tab$ci) / n

  seen <- 0L
  mismatches <- character(0)
  for (i in which(ok)) {
    cells <- as.numeric(tab[i, c("ai", "bi", "ci", "di")])
    row <- printed_row(as.list(cells), n, dp_marg, dp_kappa)
    hit <- ok &
      rounds_to(na, row$p1i) &
      rounds_to(nb, row$p2i) &
      rounds_to(kap, row$kappa)
    hit[is.na(hit)] <- FALSE
    brute <- as_key(tab[hit, , drop = FALSE])
    got <- as_key(recover2x2(n, p1i = row$p1i, p2i = row$p2i,
                             kappa = row$kappa)$tables)
    if (!identical(got, brute)) {
      mismatches <- c(mismatches, sprintf(
        "%s: enumerated {%s}, compatible {%s}", paste(cells, collapse = "/"),
        paste(setdiff(got, brute), collapse = ","),
        paste(setdiff(brute, got), collapse = ",")))
    }
    seen <- seen + 1L
  }
  list(seen = seen, mismatches = mismatches, tables = nrow(tab))
}

test_that("at n = 18 the candidate set is exactly the compatible set", {
  res <- sweep_is_exactly_the_compatible_set(18, 2, 1)
  expect_equal(res$tables, 1330)
  expect_equal(res$seen, 1258)
  expect_equal(head(res$mismatches, 3), character(0))
})

test_that("at n = 22 the candidate set is exactly the compatible set", {
  res <- sweep_is_exactly_the_compatible_set(22, 3, 2)
  expect_equal(res$tables, 2300)
  expect_equal(res$seen, 2212)
  expect_equal(head(res$mismatches, 3), character(0))
})

test_that("the brute-force comparison would catch a dropped or an extra table", {
  # The check above is only worth what it would reject. Take a row whose printed
  # digits leave more than one table, then drop one and add a stranger, and
  # confirm each is caught. Marginals printed to one decimal are what leaves the
  # row open at this sample size.
  n <- 18
  tab <- all_tables_of(n)
  ok <- non_degenerate(tab, n)
  kap <- rep(NA_real_, nrow(tab))
  kap[ok] <- kappa_of(tab[ok, , drop = FALSE], n)
  row <- list(p1i = "0.2", p2i = "0.4", kappa = "0.2")

  hit <- ok & rounds_to((tab$ai + tab$bi) / n, row$p1i) &
    rounds_to((tab$ai + tab$ci) / n, row$p2i) & rounds_to(kap, row$kappa)
  hit[is.na(hit)] <- FALSE
  brute <- as_key(tab[hit, , drop = FALSE])
  got <- as_key(recover2x2(n, p1i = row$p1i, p2i = row$p2i,
                           kappa = row$kappa)$tables)

  expect_identical(got, brute)
  expect_gt(length(brute), 1)
  expect_false(identical(got[-1], brute))
  expect_false(identical(sort(c(got, "0/0/0/18")), brute))
})
