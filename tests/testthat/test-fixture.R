# The backbone: this implementation against the Python one, case by case.
#
# Every case the Python test suite states by hand, every table in the exhaustive
# sweeps it runs at N <= 24, and a grid at the sample sizes real reports use.
# The assertion is on the whole outcome -- the status, the reason, the count,
# and the cells of every surviving table in the order they were enumerated -- so
# a divergence anywhere shows up here rather than in a summary statistic.

test_that("the fixture is present and covers every outcome", {
  cases <- fixture_cases()
  expect_gt(nrow(cases), 10000)
  expect_setequal(unique(cases$status),
                  c("unique", "set", "infeasible", "insufficient", "error"))
  for (s in c("unique", "set", "infeasible", "insufficient", "error")) {
    expect_gt(sum(cases$status == s), 0)
  }
  # Every table the fixture records belongs to a case it records.
  expect_true(all(fixture_tables()$id %in% cases$id))
  expect_equal(sum(cases$ntables), nrow(fixture_tables()))
})

test_that("every fixture case returns the Python status, count, reason and tables", {
  cases <- fixture_cases()
  tables <- fixture_tables()
  by_id <- split(tables[, c("ai", "bi", "ci", "di")], tables$id)

  bad_status <- character(0)
  bad_count <- character(0)
  bad_tables <- character(0)
  bad_reason <- character(0)

  for (i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    got <- run_fixture_case(row)
    if (!identical(got$status, row$status)) {
      bad_status <- c(bad_status, sprintf("%s: %s != %s", row$group,
                                          got$status, row$status))
      next
    }
    if (nrow(got$tables) != row$ntables) {
      bad_count <- c(bad_count, sprintf("%s: %d != %g", row$group,
                                        nrow(got$tables), row$ntables))
      next
    }
    want <- by_id[[as.character(row$id)]]
    if (is.null(want)) {
      want <- data.frame(ai = numeric(0), bi = numeric(0),
                         ci = numeric(0), di = numeric(0))
    }
    rownames(want) <- NULL
    if (!isTRUE(all.equal(as.matrix(got$tables), as.matrix(want),
                          check.attributes = FALSE))) {
      bad_tables <- c(bad_tables, row$group)
    }
    if (!identical(got$reason, python_reason_as_r(row$reason))) {
      bad_reason <- c(bad_reason, sprintf("%s: '%s' != '%s'", row$group,
                                          got$reason, row$reason))
    }
  }

  expect_equal(head(bad_status, 5), character(0))
  expect_equal(head(bad_count, 5), character(0))
  expect_equal(head(bad_tables, 5), character(0))
  expect_equal(head(bad_reason, 5), character(0))
})
