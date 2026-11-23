# Build IPF_cohort_MatriSpace.rds from the raw BioStudies S-BSST1410 files.

library(Seurat)
library(Matrix)
library(jsonlite)

spaceranger_dir <- "data/spaceranger"
cell2location_dir <- "data/cell2location_habermann2020"
metadata_file <- "data/hs_visium_metadata.tsv"
histology_file <- "data/hs_visium_merged_histo_annotations.csv"

cell_types <- c("AT1", "AT2", "Transitional AT2", "KRT5-/KRT17+", "Myofibroblasts",
                "HAS1 High Fibroblasts", "Fibroblasts", "Endothelial Cells", "Macrophages",
                "B Cells", "Plasma Cells", "T Cells")
haemoglobin <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

md <- read.delim(metadata_file, check.names = FALSE)
histology <- read.csv(histology_file)

data_list <- list()
meta_list <- list()
coord_list <- list()

for (i in seq_len(nrow(md))) {
  sid <- md$sample_id[i]
  path <- file.path(spaceranger_dir, sid)

  counts <- Read10X_h5(file.path(path, "raw_feature_bc_matrix.h5"))
  pos <- read.csv(file.path(path, "spatial/tissue_positions_list.csv"), header = FALSE,
                  col.names = c("barcode", "in_tissue", "row", "col", "pixel_row", "pixel_col"))
  pos <- pos[pos$in_tissue == 1, ]
  counts <- counts[, pos$barcode]

  n_count <- colSums(counts)
  n_gene <- colSums(counts > 0)
  pct_mt <- 100 * colSums(counts[grepl("^MT-", rownames(counts)), ]) / n_count
  pct_hb <- 100 * colSums(counts[rownames(counts) %in% haemoglobin, ]) / n_count
  keep <- n_count >= 500 & n_gene >= 200 & pct_mt < 25 & pct_hb < 20
  counts <- counts[, keep]
  pos <- pos[keep, ]

  spots <- paste(sid, pos$barcode, sep = "_")
  colnames(counts) <- spots
  data_list[[sid]] <- log1p(counts %*% Diagonal(x = 1e4 / n_count[keep]))
  dimnames(data_list[[sid]]) <- dimnames(counts)

  sf <- fromJSON(file.path(path, "spatial/scalefactors_json.json"))
  um_per_px <- 55 / sf$spot_diameter_fullres
  coord_list[[sid]] <- data.frame(x = pos$pixel_col * um_per_px, y = -pos$pixel_row * um_per_px,
                                  cell = spots)

  c2l <- read.csv(file.path(cell2location_dir, paste0(sid, "_spot_cell_abundances_5pc.csv")),
                  check.names = FALSE)
  ab <- as.matrix(c2l[match(paste0(sid, "_", pos$barcode), c2l$spot_id), -1])
  colnames(ab) <- sub("q05cell_abundance_w_sf_", "", colnames(ab), fixed = TRUE)
  prop <- ab[, cell_types] / rowSums(ab)
  colnames(prop) <- paste0("prop_", gsub("[^A-Za-z0-9]+", "_", cell_types))

  meta_list[[sid]] <- data.frame(
    section = sid,
    block = paste(md$subject_alias[i], md$tissue_alias[i], sep = "_"),
    lung = md$subject_alias[i],
    condition = md$condition[i],
    histology = histology$annotation[match(paste(sid, pos$barcode),
                                           paste(histology$sample_slide_id, histology$barcode))],
    nCount_Spatial = n_count[keep],
    percent_mt = pct_mt[keep],
    prop,
    row.names = spots
  )
}

genes <- Reduce(intersect, lapply(data_list, rownames))
expr <- do.call(cbind, lapply(data_list, function(x) x[genes, ]))
obj <- CreateSeuratObject(CreateAssay5Object(data = expr), assay = "Spatial")
obj <- AddMetaData(obj, do.call(rbind, unname(meta_list)))

for (k in seq_along(coord_list)) {
  fov <- CreateFOV(coord_list[[k]], type = "centroids", assay = "Spatial", key = paste0("s", k, "_"))
  obj[[gsub("-", "_", names(coord_list)[k])]] <- fov
}

saveRDS(obj, "IPF_cohort_MatriSpace.rds")
