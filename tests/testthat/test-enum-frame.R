# The data frame form, and where it differs from metafor::conv.2x2.
#
# The interface is deliberately the one a conv.2x2 user already knows, so what
# is asserted here is that the call shape carries over and that the answer does
# not: conv.2x2 always fills ai, bi, ci, di, while enum.2x2 fills them only
# where the printed digits determine them and reports the span otherwise.

corpus <- function() {
  data.frame(study = c("delirium pooled", "cohort B", "cohort C"),
             ni = c(768, 20306, 370),
             n1i = c(158, 866, 165),
             n2i = c(466, 1603, 160),
             kappa = c("0.29", "0.22", "0.48"),
             stringsAsFactors = FALSE)
}

test_that("arguments are evaluated inside the data frame, as conv.2x2 does", {
  dat <- corpus()
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat)
  expect_equal(nrow(out), 3)
  expect_true(all(names(dat) %in% names(out)))
  expect_equal(out$study, dat$study)
  expect_equal(out$status, c("unique", "set", "unique"))
  expect_equal(out$ntables, c(1L, 11L, 1L))
})

test_that("the point columns are filled only where the digits determine them", {
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = corpus())
  expect_equal(out$ai, c(158, NA, 115))
  expect_equal(out$bi, c(0, NA, 50))
  expect_equal(out$ci, c(308, NA, 45))
  expect_equal(out$di, c(302, NA, 160))
  # Rows 1 and 3 are determined; only row 2, on a sample of 20,306 with kappa
  # printed to two decimals, is left open.
  expect_equal(out$ntables, c(1L, 11L, 1L))
  expect_true(is.na(out$ai[2]))
  expect_false(any(is.na(out[1, c("ai", "bi", "ci", "di")])))
})

test_that("the range columns span the compatible set and contain the point", {
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = corpus())
  for (cell in c("ai", "bi", "ci", "di")) {
    lo <- out[[paste0(cell, ".min")]]
    hi <- out[[paste0(cell, ".max")]]
    expect_true(all(lo <= hi))
    determined <- !is.na(out[[cell]])
    expect_equal(lo[determined], out[[cell]][determined])
    expect_equal(hi[determined], out[[cell]][determined])
  }
  expect_equal(out$ai.min[2], 320)
  expect_equal(out$ai.max[2], 330)
  expect_gt(out$ai.max[2] - out$ai.min[2], 0)
})

test_that("the full set for any row is reachable from the returned frame", {
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = corpus())
  recs <- attr(out, "recoveries")
  expect_length(recs, 3)
  expect_s3_class(recs[[2]], "enum2x2")
  expect_equal(n_tables(recs[[2]]), out$ntables[2])
  expect_equal(min(recs[[2]]$tables$ai), out$ai.min[2])
})

test_that("var.names renames the point columns and the ranges with them", {
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = corpus(),
                  var.names = c("both", "first", "second", "neither"))
  expect_true(all(c("both", "both.min", "both.max", "neither.max") %in% names(out)))
  expect_false("ai" %in% names(out))
  expect_equal(out$both[1], 158)
  expect_equal(out$second[1], 308)
  expect_error(enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i,
                        data = corpus(), var.names = c("a", "b")),
               "four distinct names")
})

test_that("append and replace behave as conv.2x2 documents", {
  dat <- corpus()
  bare <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat,
                   append = FALSE)
  expect_false("study" %in% names(bare))
  expect_true("status" %in% names(bare))

  # An existing cell column is filled only where it is missing, unless told
  # otherwise.
  dat$ai <- c(999, NA, NA)
  keep <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat)
  expect_equal(keep$ai[1], 999)
  over <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat,
                   replace = "all")
  expect_equal(over$ai[1], 158)
  expect_equal(enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat,
                        replace = TRUE)$ai[1], 158)
})

test_that("include restricts which rows are enumerated", {
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = corpus(),
                  include = c(TRUE, FALSE, TRUE))
  expect_equal(out$status, c("unique", NA, "unique"))
  expect_true(is.na(out$ntables[2]))
  by_index <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i,
                       data = corpus(), include = c(1, 3))
  expect_equal(by_index$status, c("unique", NA, "unique"))
})

test_that("the scalar form takes one comparison without a data frame", {
  out <- enum.2x2(kappa = "0.29", ni = 768, n1i = 158, n2i = 466)
  expect_equal(nrow(out), 1)
  expect_equal(out$status, "unique")
  expect_equal(out$ai, 158)
  expect_equal(out$di.min, 302)
  # Variables in the calling frame are found the way conv.2x2 finds them.
  k <- "0.29"
  n <- 768
  expect_equal(enum.2x2(kappa = k, ni = n, n1i = 158, n2i = 466)$ai, 158)
})

test_that("printed marginals and an agreement column are accepted", {
  dat <- data.frame(ni = c(768, 20306), p1i = c("20.6", "4.26"),
                    p2i = c("60.7", "7.89"), kappa = c("0.29", "0.22"),
                    stringsAsFactors = FALSE)
  out <- enum.2x2(kappa = kappa, ni = ni, p1i = p1i, p2i = p2i, data = dat,
                  marginals.as.percent = TRUE)
  expect_equal(out$status, c("unique", "set"))
  expect_equal(out$ai[1], 158)

  with_po <- enum.2x2(kappa = "0.29", ni = 768, n1i = 158, n2i = 466,
                      agreement = "0.599")
  expect_equal(with_po$ntables, 1L)
})

test_that("a row that cannot describe a table does not stop the frame", {
  dat <- data.frame(ni = c(768, 100, 768), n1i = c(158, 101, 158),
                    n2i = c(466, 50, NA), kappa = c("0.29", "0.5", "0.29"),
                    stringsAsFactors = FALSE)
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat)
  expect_equal(out$status, c("unique", "impossible", "insufficient"))
  expect_match(out$reason[2], "outside")
  expect_equal(out$ntables, c(1L, 0L, 0L))
})

test_that("a numeric kappa column is refused, in the frame form too", {
  dat <- data.frame(ni = 768, n1i = 158, n2i = 466, kappa = 0.29)
  expect_error(enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat),
               "the strings the sources printed")
})

test_that("the frame form and the scalar form agree row by row", {
  dat <- corpus()
  out <- enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat)
  for (i in seq_len(nrow(dat))) {
    r <- recover2x2(dat$ni[i], n1i = dat$n1i[i], n2i = dat$n2i[i],
                    kappa = dat$kappa[i])
    expect_equal(out$status[i], r$status)
    expect_equal(out$ntables[i], n_tables(r))
    expect_equal(out$ai.min[i], min(r$tables$ai))
    expect_equal(out$di.max[i], max(r$tables$di))
  }
})

test_that("conv.2x2 fills the cells this leaves open, and gets them wrong", {
  # The contrast the package exists to make. metafor reconstructs one table from
  # the phi coefficient of the same study; enum.2x2 reports the span the printed
  # kappa leaves. The point estimate lands inside the span but is not the table.
  skip_if_not_installed("metafor")
  cells <- c(ai = 326, bi = 540, ci = 1277, di = 18163)
  n <- sum(cells)
  phi <- (cells[["ai"]] * cells[["di"]] - cells[["bi"]] * cells[["ci"]]) /
    sqrt((cells[["ai"]] + cells[["bi"]]) * (cells[["ci"]] + cells[["di"]]) *
           (cells[["ai"]] + cells[["ci"]]) * (cells[["bi"]] + cells[["di"]]))
  point <- metafor::conv.2x2(ri = round(phi, 2), ni = n,
                             n1i = cells[["ai"]] + cells[["bi"]],
                             n2i = cells[["ai"]] + cells[["ci"]])
  set <- enum.2x2(kappa = sprintf("%.2f", kappa2x2(cells[1], cells[2],
                                                   cells[3], cells[4])),
                  ni = n, n1i = cells[["ai"]] + cells[["bi"]],
                  n2i = cells[["ai"]] + cells[["ci"]])
  expect_equal(nrow(point), 1)
  expect_false(is.na(point$ai))
  # The enumerated span contains the true table; the point estimate is a single
  # number with no indication of how far off it is.
  expect_lte(set$ai.min, cells[["ai"]])
  expect_gte(set$ai.max, cells[["ai"]])
  expect_gt(set$ai.max - set$ai.min, 0)
  expect_true(is.na(set$ai))
})
