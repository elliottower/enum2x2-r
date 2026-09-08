# enum2x2 0.1.0

First release, a port of the Python package of the same name.

- `recover2x2()` enumerates every integer 2x2 table whose statistics round back to
  the digits a source printed, given a sample size, both positive marginals and
  Cohen's kappa. Membership is decided by exact integer cross-multiplication; there
  is no tolerance parameter and no package dependency.
- `enum.2x2()` is the data frame form, shaped after `metafor::conv.2x2()`: the same
  `data`, `include`, `var.names`, `append` and `replace` arguments, and the same four
  cell columns, filled only where the printed digits determine them and accompanied
  by the span of each cell and the number of compatible tables.
- Outcomes are reported as unique, set-identified, infeasible, or insufficient
  inputs, so a source whose published figures admit no common table is distinguished
  from one whose rounding leaves several.
- `kappa2x2()`, `exact_kappa()`, `kappa_max()`, `kappa_min()`, `expected_agreement()`,
  `agreement2x2()`, `asymmetry2x2()`, `exact_interval()`, `rounding_interval()`,
  `rounds_to()`, `exactly_rounds_to()` and `counts_rounding_to()` are exposed for use
  on their own.
