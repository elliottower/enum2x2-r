# Boundaries, degenerate margins, and the endpoint convention.

test_that("printed digits fix an interval, not a value", {
  expect_equal(rounding_interval("0.1"), c(0.05, 0.15))
  expect_equal(rounding_interval("0.10"), c(0.095, 0.105))
  expect_equal(rounding_interval("0.100"), c(0.0995, 0.1005))
  expect_equal(rounding_interval("4.26", as.percent = TRUE), c(0.04255, 0.04265))
  expect_equal(rounding_interval("-0.20"), c(-0.205, -0.195))
  # A trailing zero is information, so the two are not the same input.
  expect_false(identical(rounding_interval("0.1"), rounding_interval("0.10")))
})

test_that("both endpoints are admitted, because the convention is not stated", {
  # 0.105 rounds to 0.10 under half-even and to 0.11 under half-up, and a source
  # does not say which it used.
  expect_true(rounds_to(0.105, "0.10"))
  expect_true(rounds_to(0.105, "0.11"))
  expect_false(rounds_to(0.12, "0.10"))
})

test_that("a statistic exactly on an endpoint survives either printing", {
  # The table 2/0/2/20 at n = 24 has kappa exactly 5/8, which a source prints as
  # "0.62" under half-even and "0.63" under half-up. A half-open interval loses
  # the table under one of the two.
  k <- exact_kappa(2, 0, 2, 20)
  expect_equal(k, list(num = 5, den = 8))
  expect_equal(kappa2x2(2, 0, 2, 20), 0.625)
  expect_true(exactly_rounds_to(k$num, k$den, "0.62"))
  expect_true(exactly_rounds_to(k$num, k$den, "0.63"))
  expect_equal(unique_table(recover2x2(24, n1i = 2, n2i = 4, kappa = "0.62")),
               data.frame(ai = 2, bi = 0, ci = 2, di = 20))
  expect_equal(unique_table(recover2x2(24, n1i = 2, n2i = 4, kappa = "0.63")),
               data.frame(ai = 2, bi = 0, ci = 2, di = 20))
})

test_that("membership is exact where a tolerance would be wrong", {
  # 59/200 is the upper endpoint of "0.29". A value one part in 10^13 above it is
  # outside the interval and is excluded, although the slack in rounds_to admits
  # it. The enumeration uses the exact test, which is why there is no tolerance
  # to choose.
  num <- 590000000000020
  den <- 2000000000000000
  expect_gt(num / den, 0.295)
  expect_true(rounds_to(num / den, "0.29"))
  expect_false(exactly_rounds_to(num, den, "0.29"))
  expect_true(exactly_rounds_to(59, 200, "0.29"))
  expect_true(exactly_rounds_to(57, 200, "0.29"))
  expect_false(exactly_rounds_to(5699, 20000, "0.29"))
  expect_true(exactly_rounds_to(5700, 20000, "0.29"))
})

test_that("an interval declared to twelve decimals is still decided exactly", {
  iv <- exact_interval("0.000000000000")
  expect_equal(iv, list(lo_num = -1, hi_num = 1, den = 2e12))
  expect_true(exactly_rounds_to(1, 2e12, "0.000000000000"))
  expect_false(exactly_rounds_to(1000001, 2e18, "0.000000000000"))
  expect_false(exactly_rounds_to(1, 1e12, "0.000000000000"))
})

test_that("the top of a marginal interval is not truncated away", {
  # At n = 240 the marginal "0.512" admits exactly one count, 123, and
  # int(0.5125 * 240) is 122, so a truncating upper bound loses it and the row
  # wrongly reports that no table exists.
  expect_equal(counts_rounding_to("0.512", 240), 123)
  expect_equal(counts_rounding_to("0.14", 200), c(27, 28, 29))
  expect_true(63 %in% counts_rounding_to("0.3", 180))
  r <- recover2x2(240, p1i = "0.512", n2i = 100, kappa = "0.18")
  expect_true(r$status %in% c("unique", "set"))
  expect_true(all(r$tables$ai + r$tables$bi == 123))
})

test_that("a marginal of zero or of the whole sample is handled, not special-cased", {
  # A printed marginal of "0.0" admits a candidate pair whose expected agreement
  # is exactly one, leaving kappa undefined for that pair alone. The row still
  # has compatible tables and returns them.
  r <- recover2x2(24, p1i = "0.0", p2i = "0.0", kappa = "0.00")
  expect_setequal(as_key(r$tables), c("0/0/1/23", "0/1/0/23"))

  # The same at the top of the range: every case is negative on neither.
  s <- recover2x2(24, p1i = "1.0", p2i = "1.0", kappa = "0.00")
  expect_setequal(as_key(s$tables), c("23/1/0/0", "23/0/1/0"))

  # An exact count of zero is a legitimate marginal, not a defect in the call:
  # it is reported as admitting no table rather than signaled.
  z <- recover2x2(24, n1i = 0, n2i = 0, kappa = "0.00")
  expect_equal(z$status, "infeasible")
  expect_equal(n_tables(z), 0)
  expect_equal(n_tables(recover2x2(24, n1i = 24, n2i = 24, kappa = "0.00")), 0)
  # One margin at zero and the other not still leaves kappa defined.
  expect_equal(n_tables(recover2x2(24, n1i = 0, n2i = 4, kappa = "0.00")), 1)
})

test_that("kappa is undefined where expected agreement is one", {
  expect_error(kappa2x2(100, 0, 0, 0), class = "enum2x2_undefined_statistic")
  expect_error(exact_kappa(100, 0, 0, 0), class = "enum2x2_undefined_statistic")
  expect_error(kappa2x2(0, 0, 0, 100), class = "enum2x2_undefined_statistic")
  expect_error(kappa_max(0, 0), class = "enum2x2_undefined_statistic")
  expect_error(kappa_min(1, 1), class = "enum2x2_undefined_statistic")
  expect_error(kappa2x2(0, 0, 0, 0), class = "enum2x2_invalid_input")
})

test_that("kappa reaches its marginal maximum and no further", {
  # Equal marginals permit kappa = 1 exactly; unequal marginals cap it below 1,
  # and a published value above the cap admits no table.
  expect_equal(kappa_max(0.3, 0.3), 1)
  expect_equal(kappa_max(0.5, 0.5), 1)
  expect_lt(kappa_max(0.0426, 0.0789), 1)
  for (p in list(c(0.0426, 0.0789), c(0.011, 0.017), c(0.446, 0.432), c(0.7, 0.2))) {
    p_e <- expected_agreement(p[1], p[2])
    expect_equal(kappa_max(p[1], p[2]), ((1 - abs(p[1] - p[2])) - p_e) / (1 - p_e))
  }
  reachable <- recover2x2(1000, p1i = "0.200", p2i = "0.400",
                          kappa = sprintf("%.2f", kappa_max(0.2, 0.4)))
  expect_gt(n_tables(reachable), 0)
  expect_equal(n_tables(recover2x2(1000, p1i = "0.200", p2i = "0.400",
                                   kappa = "1.00")), 0)
})

test_that("kappa reaches its marginal minimum, which the marginals also fix", {
  for (p in list(c(0.3, 0.3), c(0.0426, 0.0789), c(0.7, 0.2), c(0.9, 0.9))) {
    expect_lte(kappa_min(p[1], p[2]), 0)
    expect_gte(kappa_max(p[1], p[2]), 0)
    expect_lt(kappa_min(p[1], p[2]), kappa_max(p[1], p[2]))
  }
  # Built to minimise agreement: as little overlap as the marginals allow.
  n <- 1000
  n1 <- 300
  n2 <- 400
  ai <- max(0, n1 + n2 - n)
  expect_equal(kappa2x2(ai, n1 - ai, n2 - ai, n - n1 - n2 + ai),
               kappa_min(0.3, 0.4), tolerance = 1e-9)
})

test_that("the two marginals may be printed to different precisions", {
  fine <- recover2x2(768, p1i = "0.2057", p2i = "0.6068", kappa = "0.29")
  coarse <- recover2x2(768, p1i = "0.21", p2i = "0.61", kappa = "0.29")
  expect_lte(n_tables(fine), n_tables(coarse))
  mixed <- recover2x2(768, p1i = "0.2057", p2i = "0.61", kappa = "0.29")
  expect_lte(n_tables(fine), n_tables(mixed))
  expect_lte(n_tables(mixed), n_tables(coarse))
  # A count on one side and a printed proportion on the other.
  expect_gt(n_tables(recover2x2(768, n1i = 158, p2i = "0.6068", kappa = "0.29")), 0)
})

test_that("figures that admit no table are reported with the figure that excluded", {
  r <- recover2x2(370, n1i = 165, n2i = 160, kappa = "0.48",
                  agreement = "73", agreement.as.percent = TRUE)
  expect_equal(r$status, "infeasible")
  expect_match(r$reason, "adding the published agreement admits none")
  expect_equal(n_tables(r), 0)

  s <- recover2x2(768, n1i = 158, n2i = 466, kappa = "0.95")
  expect_equal(s$status, "infeasible")
  expect_match(s$reason, "range these marginals permit")

  # Marginals of 2% and 90% cannot support kappa 0.99.
  expect_equal(n_tables(recover2x2(1000, p1i = "0.0200", p2i = "0.9000",
                                    kappa = "0.99")), 0)
})
