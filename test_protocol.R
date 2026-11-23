# Run from the folder containing the protocol scripts.
source("cohort_helpers.R")

# Equal weighting of sections within blocks and blocks within lungs.
unbalanced <- data.frame(
  section = c("s1", "s1", "s2", "s3", "s4"),
  block = c("b1", "b1", "b1", "b2", "b3"),
  lung = c("L1", "L1", "L1", "L1", "L2"),
  condition = c("control", "control", "control", "control", "IPF"),
  value = c(0, 2, 3, 10, 20)
)
means <- summarise_by_lung(unbalanced, "value")
stopifnot(identical(means$value, c(6, 20)))

# Complete separation of four versus four lungs has a two-sided P of 2/70.
donors <- data.frame(condition = rep(c("control", "IPF"), each = 4), value = 1:8)
stats <- permutation_test(donors, "value")
stopifnot(abs(stats$p_value - 2/70) < 1e-12)

# Plot whiskers extend to the range, as stated in the figure legends.
p <- plot_by_condition(donors, "value")
b <- ggplot2::ggplot_build(p)$data[[1]]
stopifnot(identical(b$ymin, c(1, 5)), identical(b$ymax, c(4, 8)))

# Matrix order is retained in the figure heatmaps.
m <- matrix(1:6, 2, dimnames = list(c("B", "A"), c("3", "1", "2")))
h <- matrispace::plot_heatmap(m)
stopifnot(identical(levels(h$data$row), rownames(m)),
          identical(levels(h$data$column), colnames(m)))
cat("Protocol checks passed.\n")
