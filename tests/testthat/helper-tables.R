# Every integer 2x2 table of a given size, and what a paper reporting one would
# print. Used by the brute-force, independent-implementation and round-trip
# checks, none of which may take the enumerator's word for anything.

all_tables_of <- function(n) {
  g <- expand.grid(ai = 0:n, bi = 0:n, ci = 0:n)
  g$di <- n - g$ai - g$bi - g$ci
  g <- g[g$di >= 0, , drop = FALSE]
  rownames(g) <- NULL
  g
}

# The four margins strictly inside (0, n). Kappa is undefined outside that, and
# the sweeps skip those tables rather than special-casing them.
non_degenerate <- function(tab, n) {
  na <- tab$ai + tab$bi
  nb <- tab$ai + tab$ci
  na > 0 & na < n & nb > 0 & nb < n
}

# Kappa for a whole frame of tables at once, from the definition, in doubles.
# The margins are already known to be strictly inside (0, n), so expected
# agreement is below one and the quotient is defined throughout.
kappa_of <- function(tab, n) {
  p_o <- (tab$ai + tab$di) / n
  p_e <- ((tab$ai + tab$bi) / n) * ((tab$ai + tab$ci) / n) +
    (1 - (tab$ai + tab$bi) / n) * (1 - (tab$ai + tab$ci) / n)
  (p_o - p_e) / (1 - p_e)
}

# What a paper reporting this table would print.
printed_row <- function(cells, n, dp_marg, dp_kappa) {
  list(p1i = sprintf("%.*f", dp_marg, (cells[[1]] + cells[[2]]) / n),
       p2i = sprintf("%.*f", dp_marg, (cells[[1]] + cells[[3]]) / n),
       kappa = sprintf("%.*f", dp_kappa,
                       kappa2x2(cells[[1]], cells[[2]], cells[[3]], cells[[4]])))
}

as_key <- function(df) {
  if (nrow(df) == 0L) return(character(0))
  sort(paste(df$ai, df$bi, df$ci, df$di, sep = "/"))
}
