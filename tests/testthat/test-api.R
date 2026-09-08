# The public surface: what a caller touches, and how it fails.
#
# A figure that cannot describe any table is a defect in the call and is
# signalled. A source that reported too little, or reported figures that do not
# cohere, is described in the returned object. Keeping those apart is the whole
# reason the object has a status, so both directions are asserted here.

DELIRIUM <- list(768, n1i = 158, n2i = 466, kappa = "0.29")

test_that("a marginal outside the sample is refused", {
  expect_error(recover2x2(100, n1i = 101, n2i = 50, kappa = "0.5"),
               class = "enum2x2_invalid_input")
  expect_error(recover2x2(100, n1i = -1, n2i = 50, kappa = "0.5"), "outside")
  expect_error(recover2x2(100, n1i = 50, n2i = 101, kappa = "0.5"), "outside")
  expect_error(recover2x2(100, n1i = 50.5, n2i = 50, kappa = "0.5"),
               "integer count")
  expect_error(recover2x2(100, n1i = TRUE, n2i = 50, kappa = "0.5"),
               "integer count")
})

test_that("a kappa beyond the unit interval is refused", {
  expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = "1.4"),
               "outside \\[-1, 1\\]")
  expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = "-1.4"),
               "outside \\[-1, 1\\]")
  # Exactly at the endpoints is admissible; it is the marginals that then decide.
  expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = "1.0"), NA)
  expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = "-1.0"), NA)
})

test_that("a non-positive or non-integer sample is refused", {
  expect_error(recover2x2(0, n1i = 0, n2i = 0, kappa = "0.5"), "positive integer")
  expect_error(recover2x2(-5, n1i = 0, n2i = 0, kappa = "0.5"), "positive integer")
  expect_error(recover2x2(7.5, n1i = 0, n2i = 0, kappa = "0.5"), "positive integer")
  expect_error(recover2x2(TRUE, n1i = 1, n2i = 1, kappa = "0.5"), "positive integer")
  expect_error(recover2x2(c(10, 20), n1i = 1, n2i = 1, kappa = "0.5"),
               "positive integer")
})

test_that("a number where a printed string belongs is refused", {
  # 0.10 and 0.1 are the same double and imply different intervals, so the
  # distinction cannot be recovered after coercion.
  expect_error(recover2x2(768, n1i = 158, n2i = 466, kappa = 0.29),
               "the string the source printed")
  expect_error(recover2x2(768, p1i = 0.206, n2i = 466, kappa = "0.29"),
               "the string the source printed")
  expect_error(recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29",
                          agreement = 0.6), "the string the source printed")
})

test_that("a figure that is not a decimal number is refused", {
  for (bad in c("abc", "", "one half")) {
    expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = bad),
                 "not a decimal number")
  }
  for (bad in c("nan", "inf", "-inf", "Infinity", "NaN")) {
    expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = bad),
                 "not a finite number")
  }
  # A decimal parse alone would accept this one; it is the range that rejects it.
  expect_error(recover2x2(100, n1i = 50, n2i = 50, kappa = "1e400"),
               "outside \\[-1, 1\\]")
})

test_that("two marginals for one criterion is a defect in the call", {
  expect_error(recover2x2(768, n1i = 158, p1i = "20.6", n2i = 466, kappa = "0.29",
                          marginals.as.percent = TRUE), "not both")
  expect_error(recover2x2(768, n1i = 158, n2i = 466, p2i = "60.7", kappa = "0.29",
                          marginals.as.percent = TRUE), "not both")
})

test_that("a printed marginal that is not a proportion is refused", {
  expect_error(recover2x2(100, p1i = "150", p2i = "50", kappa = "0.5",
                          marginals.as.percent = TRUE),
               "outside \\[0, 1\\] as a proportion")
  expect_error(recover2x2(100, p1i = "-50", p2i = "50", kappa = "0.5",
                          marginals.as.percent = TRUE), "as a percentage")
  expect_error(recover2x2(100, p1i = "1e400", n2i = 50, kappa = "0.5"),
               "outside \\[0, 1\\] as a proportion")
})

test_that("a source that reported too little is described, not signalled", {
  r <- recover2x2(768, n1i = 158, n2i = 466)
  expect_equal(r$status, "insufficient")
  expect_match(r$reason, "kappa was not reported")

  s <- recover2x2(768, n1i = 158, kappa = "0.29")
  expect_equal(s$status, "insufficient")
  expect_match(s$reason, "n2i")

  t <- recover2x2(768, n2i = 466, kappa = "0.29")
  expect_match(t$reason, "n1i")

  # A missing kappa is reported before a missing marginal, because kappa is the
  # figure the enumeration inverts.
  expect_match(recover2x2(768)$reason, "kappa was not reported")
})

test_that("a determined report returns one table and names it", {
  r <- do.call(recover2x2, DELIRIUM)
  expect_equal(r$status, "unique")
  expect_equal(n_tables(r), 1)
  expect_equal(unique_table(r), data.frame(ai = 158, bi = 0, ci = 308, di = 302))
  expect_equal(r$ni, 768)
  expect_true(is.na(r$reason))
  expect_equal(r$published$kappa, "0.29")
  expect_equal(r$published$n1i, 158)
  expect_null(r$published$p1i)
})

test_that("an underdetermined report returns the set and its ranges", {
  r <- recover2x2(20306, n1i = 866, n2i = 1603, kappa = "0.22")
  expect_equal(r$status, "set")
  expect_equal(n_tables(r), 11)
  cr <- cell_ranges(r)
  expect_equal(cr$cell, c("ai", "bi", "ci", "di"))
  expect_equal(cr$min[cr$cell == "ai"], 320)
  expect_equal(cr$max[cr$cell == "ai"], 330)
  expect_true(all(rowSums(r$tables) == 20306))
  expect_true(all(r$tables$ai + r$tables$bi == 866))
  expect_true(all(r$tables$ai + r$tables$ci == 1603))
})

test_that("asking a set for the table says so rather than guessing", {
  r <- recover2x2(20306, n1i = 866, n2i = 1603, kappa = "0.22")
  expect_error(unique_table(r), "does not determine one table")
  expect_error(cell_ranges(recover2x2(768, n1i = 158, n2i = 466, kappa = "0.95")),
               "no cell has a range")
})

test_that("every returned table reproduces the figures it was recovered from", {
  calls <- list(DELIRIUM,
                list(20306, n1i = 866, n2i = 1603, kappa = "0.22"),
                list(7328, n1i = 174, n2i = 175, kappa = "0.668"),
                list(768, p1i = "66.4", p2i = "20.6", kappa = "0.22",
                     marginals.as.percent = TRUE))
  for (call in calls) {
    r <- do.call(recover2x2, call)
    expect_gt(n_tables(r), 0)
    for (i in seq_len(n_tables(r))) {
      cells <- as.numeric(r$tables[i, ])
      k <- exact_kappa(cells[1], cells[2], cells[3], cells[4])
      expect_true(exactly_rounds_to(k$num, k$den, r$published$kappa))
      if (!is.null(r$published$n1i)) {
        expect_equal(cells[1] + cells[2], r$published$n1i)
      } else {
        expect_true(rounds_to((cells[1] + cells[2]) / r$ni, r$published$p1i,
                              r$published$marginals.as.percent))
      }
    }
  }
})

test_that("a published agreement is applied as a constraint, not as a hint", {
  # The delirium table agrees on 460 of 768, which is 0.599 to three decimals.
  with_true <- recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29",
                          agreement = "0.599")
  expect_equal(n_tables(with_true), 1)
  with_wrong <- recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29",
                           agreement = "0.900")
  expect_equal(with_wrong$status, "infeasible")
  expect_match(with_wrong$reason, "adding the published agreement admits none")
  expect_equal(n_tables(recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29",
                                    agreement = "60",
                                    agreement.as.percent = TRUE)), 1)
})

test_that("percentage marginals agree with the counts they stand for", {
  by_count <- recover2x2(768, n1i = 510, n2i = 158, kappa = "0.22")
  by_pct <- recover2x2(768, p1i = "66.4", p2i = "20.6", kappa = "0.22",
                       marginals.as.percent = TRUE)
  expect_equal(as.data.frame(by_count), as.data.frame(by_pct))
  expect_equal(by_count$status, by_pct$status)
})

test_that("the object supports the accessors it documents", {
  r <- do.call(recover2x2, DELIRIUM)
  expect_s3_class(r, "enum2x2")
  expect_s3_class(as.data.frame(r), "data.frame")
  expect_equal(names(as.data.frame(r)), c("ai", "bi", "ci", "di"))
  expect_output(print(r), "unique")
  expect_output(print(r), "158")
  expect_output(print(recover2x2(20306, n1i = 866, n2i = 1603, kappa = "0.22")),
                "11 compatible tables")
  expect_output(print(recover2x2(768, n1i = 158, n2i = 466, kappa = "0.95")),
                "infeasible")
  expect_output(print(recover2x2(768, n1i = 158, kappa = "0.29")),
                "insufficient")
})

test_that("agreement and asymmetry describe the recovered table", {
  expect_equal(agreement2x2(158, 0, 308, 302), 460 / 768)
  expect_equal(asymmetry2x2(0, 308), 1)
  expect_equal(asymmetry2x2(50, 50), 0)
  expect_true(is.na(asymmetry2x2(0, 0)))
  expect_error(agreement2x2(0, 0, 0, 0), class = "enum2x2_invalid_input")
})

# ------------------------------------------------------------------- batches

test_that("one bad row does not stop the others", {
  out <- recover_many(list(
    list(ni = 768, n1i = 158, n2i = 466, kappa = "0.29"),
    list(ni = 100, n1i = 101, n2i = 50, kappa = "0.5"),
    list(ni = 7328, n1i = 174, n2i = 175, kappa = "0.668")))
  expect_equal(vapply(out, function(r) r$status, character(1)),
               c("unique", "impossible", "unique"))
  expect_match(out[[2]]$reason, "outside")
})

test_that("a batch survives a column the function does not take", {
  out <- recover_many(list(list(study = "Meagher 2014", ni = 768, n1i = 158,
                                n2i = 466, kappa = "0.29")))
  expect_equal(out[[1]]$status, "unique")
})

test_that("a cell that cannot be read is not counted as an under-reporting source", {
  # A figure this code cannot use and a figure the source never published are
  # different findings, and counting them together would misreport how much of
  # the literature under-reports.
  out <- recover_many(list(list(ni = 370, n1i = 500, n2i = 160, kappa = "0.48"),
                           list(ni = 768, n1i = 158, kappa = "0.29"),
                           list(ni = 100, n1i = 50, n2i = 50, kappa = "NA")))
  expect_equal(vapply(out, function(r) r$status, character(1)),
               c("impossible", "insufficient", "impossible"))
  expect_match(out[[3]]$reason, "not a decimal number")
})

test_that("a batch returns one result per row, in order, and never signals", {
  rows <- list(list(ni = 768, n1i = 158, n2i = 466, kappa = "0.29"),
               list(ni = 768, n1i = 158, kappa = "0.29"),
               list(ni = 768, n1i = 158, n2i = 466, kappa = "0.95"),
               list(ni = 20306, n1i = 866, n2i = 1603, kappa = "0.22"),
               list(n1i = 1, n2i = 1, kappa = "0.5"))
  out <- recover_many(rows)
  expect_length(out, 5)
  expect_equal(vapply(out, function(r) r$status, character(1)),
               c("unique", "insufficient", "infeasible", "set", "insufficient"))
  expect_true(all(vapply(out, function(r) r$status, character(1)) %in%
                    c("unique", "set", "infeasible", "insufficient", "impossible")))
})

test_that("a batch reads a data frame of published comparisons", {
  dat <- data.frame(ni = c(768, 20306), n1i = c(158, 866), n2i = c(466, 1603),
                    kappa = c("0.29", "0.22"), stringsAsFactors = FALSE)
  out <- recover_many(dat)
  expect_equal(vapply(out, function(r) r$status, character(1)),
               c("unique", "set"))
})
