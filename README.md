# enum2x2

**Every 2x2 table a published summary statistic still allows.**

An R port of the [enum2x2](https://github.com/elliottower/enum2x2) Python package,
with an interface shaped after `metafor::conv.2x2`. It ports the kappa path: a sample
size, both marginals and a printed kappa, with observed agreement as a further check.
The Python package also accepts McNemar, phi and other closing statistics, and figures
read as truncated.

A study comparing two binary criteria on one population computes its numbers from a
2x2 table, then publishes the numbers and not the table. The quantity a reader wants
— how many cases the two criteria classify differently, and *in which direction* — is
gone. It can often be recovered, because the same reports print the sample size, both
marginal totals and a statistic such as kappa, and together those leave a few integer
tables, often just one.

## The problem, concretely

A pooled series of 768 patients, two readings of the DSM-5 delirium criteria,
published at **60% agreement** and **kappa = 0.29**. That reads like two definitions
that disagree a fair amount.

Here is the table underneath it:

|  | relaxed + | relaxed − |
|---|---:|---:|
| **strict +** | 158 | 0 |
| **strict −** | 308 | 302 |

The strict class is nested inside the relaxed one. No patient is strict-positive and
relaxed-negative; the 308
they differ on all fall the same way, relaxed-positive and strict-negative. One reading
is simply three times wider than the other. Observed agreement is the highest these two
marginals permit, so kappa = 0.29 is the highest kappa they allow.

None of that is in "60% agreement and kappa = 0.29", which is what the paper printed.

## Install

```r
# install.packages("remotes")
remotes::install_github("elliottower/enum2x2-r")
```

Base R only, no dependencies. R >= 3.5.0.

## Use

```r
library(enum2x2)

r <- recover2x2(768, n1i = 158, n2i = 466, kappa = "0.29")
r
#> enum2x2 recovery -- unique
#>   ni = 768
#>   1 compatible table: ai = 158, bi = 0, ci = 308, di = 302

unique_table(r)
#>    ai bi  ci  di
#> 1 158  0 308 302

asymmetry2x2(0, 308)      # 1: the disagreement runs entirely one way
#> [1] 1

s <- recover2x2(20306, n1i = 866, n2i = 1603, kappa = "0.22")
s
#> enum2x2 recovery -- set-identified
#>   ni = 20306
#>   11 compatible tables
#>   ai 320-330, bi 536-546, ci 1273-1283, di 18157-18167

recover2x2(370, p1i = "44.6", p2i = "43.2", marginals.as.percent = TRUE,
           kappa = "0.48", agreement = "73", agreement.as.percent = TRUE)$reason
#> [1] "the marginals and kappa admit 1 table(s); adding the published agreement admits none"
```

**Kappa is passed as the string the source printed**, not as a number: `"0.10"` and
`"0.1"` imply different intervals and a double cannot tell them apart. A numeric kappa
is refused with that reason.

Four statuses, and they are kept apart: `unique`, `set`, `infeasible` (the published
figures admit no common table), `insufficient` (a required figure was never published).
An impossible *call* — a marginal above `ni`, a kappa outside [-1, 1], two marginals
given for one criterion — signals `enum2x2_invalid_input`. An impossible *source* is
reported, never signaled. `recover_many` and `enum.2x2` cannot signal, so they carry a
fifth status, `impossible`, for a row whose figures could not describe any table.

Cell names are `metafor`'s: `ai` positive on both criteria, `bi` positive on the first
only, `ci` positive on the second only, `di` negative on both, with `n1i = ai + bi`,
`n2i = ai + ci` and `ni = ai + bi + ci + di`.

## Why rounding is the whole problem

With exact inputs the algebra is one line. For positive rates `p1` and `p2` on `ni`
observations,

    p_e  = p1*p2 + (1 - p1)*(1 - p2)
    p_o  = kappa*(1 - p_e) + p_e
    ai/ni = (p_o + p1 + p2 - 1) / 2

But nobody prints exact inputs. A marginal printed `4.26%` and a kappa printed `0.22`
are *intervals*, so the inputs identify a **set** of integer tables. This walks every
non-negative integer table summing to `ni` whose statistics round back to the printed
strings, and reports the set. Where a source also prints its raw agreement, that enters
as a further constraint and narrows the set.

Membership is decided in exact arithmetic. A printed figure is a decimal string and a
table's statistic is a ratio of integers, so both are rational and the comparison is a
cross-multiplication of integers, with no tolerance parameter to choose. R has no
rational type, so the integers are carried in doubles, which are exact below 2^53, and
the cross-products that would exceed that go through base-10^7 limb arithmetic in
`R/exact.R`. The `gmp` package would do the same job at the cost of a hard dependency on
a compiled package.

## Against `metafor::conv.2x2`

`conv.2x2` **conv**erts a summary into a single table. This **enum**erates every table
the summary admits. The names are meant to sit next to each other, and so are the
calls:

```r
# metafor: one table per row, from a phi coefficient
metafor::conv.2x2(ri = phi, ni = ni, n1i = n1i, n2i = n2i, data = dat)

# enum2x2: every table per row, from a printed kappa
enum.2x2(kappa = kappa, ni = ni, n1i = n1i, n2i = n2i, data = dat)
```

Both take a data frame plus column arguments, both honor `data`, `include`,
`var.names`, `append` and `replace`, and both append columns named `ai`, `bi`, `ci`,
`di`. The difference is what goes in them. `conv.2x2` always fills all four.
`enum.2x2` fills them only where the printed digits determine the table, leaves them
`NA` otherwise, and adds `status`, `ntables`, and `ai.min` / `ai.max` (and the same for
the other three cells) so the reader sees the span instead of a point:

```
            study status ntables  ai ai.min ai.max
1 delirium pooled unique       1 158    158    158
2          cohort    set      11  NA    320    330
3    third report unique       1 115    115    115
```

`conv.2x2` takes an odds ratio, a phi coefficient or a chi-square — not Cohen's kappa,
which is what agreement studies report — and on rounding its documentation is candid:

> In practice, when such accuracy measures are reported, the values are typically
> rounded to some extent. This introduces inaccuracies into the reconstruction. The
> present function uses optimization methods to reconstruct the table counts so that the
> discrepancy between the reported measures and the reconstructed ones are minimized.
> **This is not guaranteed to reconstruct the actual table exactly**, but should usually
> yield a close match.

The Python package measures what "close" costs, on 1,200 known tables at N = 2,000,
9,170 and 20,306 with one statistic printed to two decimals:

| | `conv.2x2` point estimate | `enum2x2` |
|---|---|---|
| exactly right | **21.8%** | — |
| contains the true table | — | **100%** |
| largest cell error | median 2, 90th pct 8, **max 20** | — |
| uncertainty reported | none | the range itself |

`conv.2x2` is not broken; it does what it says. But it hands you one table and no
indication of how far off it is. Where a finding rests on how many patients were
reclassified in each direction, an unstated error of up to 20 patients can change it.

## What it will not do

- It does not tell you which criterion is **correct**. Neither does the table.
- It does not separate disagreement *between* criteria from inconsistent *application*
  of either. That needs the same criterion applied twice to the same cases, which no
  aggregate report contains.
- An empty candidate set says the published figures admit no common table **under the
  definitions and rounding rules used here**. It does not say an error was made — a
  denominator, an analysis set, or a version of a variable may simply have gone
  unstated.
- Three or more categories are out of scope; the identity above is the binary case.

## Validation

The Python implementation is the specification, and this port is checked against it
rather than against itself.

1. **Against the Python package, case by case.** `tests/testthat/fixtures/` holds the
   outcome the Python implementation produces for 11,296 cases: every determinate case
   in its own test suite, every integer table in the exhaustive sweeps it runs at
   N ≤ 24, and a grid at N = 200, 768 and 2,000 where coarse printing leaves a set. The
   assertion is on the whole outcome — status, reason, count, and the cells of every
   surviving table in the order they were enumerated. The fixture is written by
   `scripts/export_r_fixture.py` in the Python repository and nothing in R regenerates
   it.
2. **Exhaustively, over the whole space.** The enumerator must agree, as a *set*, with
   a brute-force sweep over every non-degenerate integer table at N = 18 and N = 22 —
   1,258 and 2,212 of them. Set equality, not membership, so over-inclusion is caught
   too.
3. **Against a second implementation.** Written from the definitions in exact rational
   arithmetic in the test file, sharing no code with the package: the interval comes
   from counting the digits of the printed string, and kappa is composed as
   `(p_o - p_e) / (1 - p_e)` out of rational operations defined there. Over every table
   at N = 16, 20 and 24.
4. **Round-trip, over sweeps rather than examples.** Every table at N = 14, 17 and 19,
   at printed precisions from coarse to fine and with and without a published
   agreement, must be in the set recovered from its own printed figures; likewise for
   draws at N up to 20,306.
5. **Rounding-convention independence.** The enumeration is repeated under endpoint
   rules that include and exclude each bound. The answer is unchanged except for tables
   whose statistic sits exactly on a bound, which is why the interval is closed at both
   ends: the table (2, 0, 2, 20) at N = 24 has kappa exactly 5/8 and a half-open rule
   loses it under one of the two printings a source might use.

## The study this was built for

The delirium table above, and other published comparisons, are analyzed in a manuscript
on what a published agreement statistic conceals. That work lives in its own repository,
with the corpus, the preregistration frozen before the enumeration was run, the results,
and the sentence each published figure was read from. This package carries only the
method.

## Citation

See `CITATION.cff`.

## License

MIT.
