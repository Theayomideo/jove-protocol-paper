# Spatial Matrisome Profiling of Fibrotic Human Lungs with MatriSpace

This repository contains the scripts used for the MatriSpace protocol on human
idiopathic pulmonary fibrosis (IPF). The protocol demonstrates analysis of an
individual IPF lung section using the web application and a cohort comparison
of four IPF and four healthy control lungs using the R package.

## Data

The processed data are deposited in Zenodo at
[10.5281/zenodo.24451144](https://doi.org/10.5281/zenodo.24451144).

- `IPF_single_section.rds`: the IPF section used in the web workflow, containing
  4,488 spots, the tissue image and cell type annotations.
- `IPF_cohort_MatriSpace.rds`: the cohort object, containing 98,681 spots across
  25 sections from eight lung donors.

The original data were published by Franzén et al. in
[Nature Genetics](https://doi.org/10.1038/s41588-024-01819-2) and are available
through BioStudies under accession
[S-BSST1410](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410).

## Single-section workflow

Open the [MatriSpace web application](https://izzilab-matrispace.share.connect.posit.cloud/)
and upload `IPF_single_section.rds`. Select `dominant_cell_group` as the
annotation and follow protocol section 1. The source sample identifier,
`IPF_4.TD042.B3.1`, is retained in the data description.

## Cohort workflow

The analysis uses R 4.6.1, Seurat 5.5.1 and MatriSpace 0.3.0. The source package
included here is the version used for this protocol. It is kept with the paper
so that the published MatriSpace repository remains unchanged.

1. Download the repository files and place `IPF_cohort_MatriSpace.rds` beside
   `cohort_analysis.Rmd`.
2. Open R in this folder and run `source("install.R")`.
3. Run `rmarkdown::render("cohort_analysis.Rmd")`.

The notebook follows protocol section 2 and produces the cohort plots, test
results and R session information. `cohort_helpers.R` contains the donor-level
summaries, permutation tests and figure settings. `cell_matrix_pairs.csv`
contains the 1,297 pairs used for the co-expression analysis.

Values are averaged across spots within sections, sections within tissue blocks,
and blocks within lungs. Each lung is one biological replicate. Spatial graphs
use a 120 micrometre radius within each section. Moran's I uses binary neighbour
weights, matching Figure 4C. Pair co-expression uses row-standardized weights.
The two top-ranked pairs are tested against the maximum fold change obtained
after repeating pair selection for all 70 assignments of the eight lungs.

## Rebuilding the cohort object

`build_cohort_object.R` reads the original BioStudies files from this layout:

- `data/spaceranger/<section>/raw_feature_bc_matrix.h5`
- `data/spaceranger/<section>/spatial/tissue_positions_list.csv`
- `data/spaceranger/<section>/spatial/scalefactors_json.json`
- `data/cell2location_habermann2020/<section>_spot_cell_abundances_5pc.csv`
- `data/hs_visium_metadata.tsv`
- `data/hs_visium_merged_histo_annotations.csv`

Run `source("build_cohort_object.R")` from the folder containing `data/`.

## MatriSpace

[Package repository](https://github.com/Theayomideo/matrispace) ·
[Web application](https://izzilab-matrispace.share.connect.posit.cloud/) ·
[MatriSpace paper](https://doi.org/10.64898/2026.04.26.720198)
