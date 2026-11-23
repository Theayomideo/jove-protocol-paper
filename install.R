# Install the package version used in this protocol.
options(repos = c(CRAN = "https://cloud.r-project.org"))
required <- c("BiocManager", "remotes", "Seurat", "rmarkdown", "knitr", "scales", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
if (!requireNamespace("UCell", quietly = TRUE)) BiocManager::install("UCell", ask = FALSE)
remotes::install_local("matrispace_0.3.0.tar.gz", dependencies = NA, upgrade = "never", force = TRUE)
