# Reproduction checks

The cohort analysis was rerun with MatriSpace 0.3.0 using the deposited cohort
object and the 1,297 cell-matrix pairs. Comparisons with the figure source values
are recorded in `figure_checks.csv`.

Scores in Figure 4B, Moran's I in Figure 4C and the selected correlations in
Figure 5B agree within 1.4e-15. Pair fold changes in Figure 5C and alveolar
expression values in Figure 5F agree within 2.1e-8. The same 14 pairs were
selected, and the selection-adjusted P value for the two leading pairs is 1/70.

Moran's I uses binary neighbour weights. This corrects the row-standardized
weights in the earlier notebook to match Figure 4C. Boxplot whiskers span the
observed range, and heatmaps retain the supplied row and column order.

The source package built and installed successfully. Package tests passed;
eight tests requiring the optional nnSVG dependency were skipped. The protocol
checks passed against the installed package. Run `source("test_protocol.R")`
from this folder to check donor weighting, permutation tests and plot settings.

R CMD check could not complete in the validation environment because its
installed-package discovery returned no packages, including base packages.
The installation and test results above were obtained separately.
