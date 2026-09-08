# The Python package as the oracle.
#
# tests/testthat/fixtures/ holds the outcome the Python implementation of
# enum2x2 produces for every case in its own test suite, plus the exhaustive
# sweeps those tests run and a grid at the sample sizes real reports use. The
# files are written by scripts/export_r_fixture.py in that repository and are
# read here; nothing in R regenerates them.

fixture_cases <- function() {
  d <- read.csv(gzfile(test_path("fixtures", "cases.csv.gz")),
                colClasses = "character", stringsAsFactors = FALSE)
  d$id <- as.numeric(d$id)
  d$ntables <- as.numeric(d$ntables)
  d
}

fixture_tables <- function() {
  d <- read.csv(gzfile(test_path("fixtures", "tables.csv.gz")),
                stringsAsFactors = FALSE)
  d[order(d$id, d$k), , drop = FALSE]
}

# The two implementations name the marginals differently: the Python package
# follows the statistic (n_a, p_b), this one follows metafor::conv.2x2 (n1i,
# p2i). Nothing else in a reason string differs, so a reason is compared after
# renaming, which keeps the assertion on the whole sentence.
python_reason_as_r <- function(x) {
  x <- gsub("n_a", "n1i", x, fixed = TRUE)
  x <- gsub("n_b", "n2i", x, fixed = TRUE)
  x <- gsub("p_a", "p1i", x, fixed = TRUE)
  x <- gsub("p_b", "p2i", x, fixed = TRUE)
  sub("^N must be", "ni must be", x)
}

# The fixture writes a figure the call never declared as this sentinel, so that a
# figure declared as an empty string stays distinguishable from an absent one.
ABSENT <- "<none>"

# Run one fixture row through the R implementation, returning the same shape the
# fixture records: a status, a reason, and the surviving tables in order.
run_fixture_case <- function(row) {
  num <- function(x) if (identical(x, ABSENT)) NULL else as.numeric(x)
  chr <- function(x) if (identical(x, ABSENT)) NULL else x
  args <- list(ni = as.numeric(row$n),
               n1i = num(row$n_a), n2i = num(row$n_b),
               p1i = chr(row$p_a), p2i = chr(row$p_b),
               kappa = chr(row$kappa), agreement = chr(row$agreement),
               marginals.as.percent = identical(row$marginals_as_percent, "TRUE"),
               agreement.as.percent = identical(row$agreement_as_percent, "TRUE"))
  args <- args[!vapply(args, is.null, logical(1))]
  tryCatch({
    r <- do.call(recover2x2, args)
    list(status = r$status,
         reason = if (is.na(r$reason)) "" else r$reason,
         tables = r$tables)
  }, enum2x2_error = function(e) {
    list(status = "error", reason = conditionMessage(e),
         tables = data.frame(ai = numeric(0), bi = numeric(0),
                             ci = numeric(0), di = numeric(0)))
  })
}
