# Fix for two-sided Wilcoxon p-values in R 4.6.1 (normal approximation).
#
# stats:::.wilcox_test_two_pval_asymp / .wilcox_test_one_pval_asymp compute the
# two-sided p-value as 2 * min(p, 1 - p) with p = lower tail. When the first group
# is larger, 1 - p rounds to 0, so every p-value below ~1e-16 is reported as 0
# (e.g. YPD.1 pleiotropy, CNVs vs SNPs: 0 instead of 2.65e-29; the reverse argument
# order gave the right value). The patch computes the upper tail directly. It
# applies to wilcox.test, pairwise.wilcox.test and ggpubr, which all call these
# functions. Exact p-values (small samples without ties) are not affected.
#
# Use: source() this file after loading packages, before any Wilcoxon test.

local({
  fixed_switch <- quote(switch(alternative, less = F(z), greater = F(z, lower.tail = FALSE),
                               two.sided = 2 * min(F(z), F(z, lower.tail = FALSE))))
  ns <- asNamespace("stats")
  for (fn in c(".wilcox_test_two_pval_asymp", ".wilcox_test_one_pval_asymp")) {
    if (!exists(fn, envir = ns, inherits = FALSE)) next   # older R: no bug
    f <- get(fn, envir = ns)
    b <- body(f)
    b[[length(b)]] <- fixed_switch
    body(f) <- b
    unlockBinding(fn, ns); assign(fn, f, envir = ns); lockBinding(fn, ns)
  }
  # self-check: the two argument orders must give the same, non-zero p-value
  set.seed(1); x <- stats::rnorm(800, 1); y <- stats::rnorm(2000)
  p_xy <- stats::wilcox.test(x, y)$p.value; p_yx <- stats::wilcox.test(y, x)$p.value
  stopifnot(p_xy > 0, isTRUE(all.equal(p_xy, p_yx)))
})
