# Helper functions for the MatriSpace cohort analysis (protocol section 2).

library(ggplot2)

`%||%` <- function(x, y) if (is.null(x)) y else x

design_cols <- c("section", "block", "lung", "condition")
condition_colours <- c(control = "#9DB4D3", IPF = "#22336B")

# ---------------------------------------------------------------- summaries and tests

# Average values per section, then per tissue block, then per lung.
summarise_by_lung <- function(x, features) {
  if (inherits(x, "Seurat")) {
    vals <- Seurat::FetchData(x, vars = features)
    x <- cbind(x[[]][, design_cols], vals)
  }
  sec <- aggregate(x[, features, drop = FALSE], x[, design_cols], mean)
  blk <- aggregate(sec[, features, drop = FALSE], sec[, c("block", "lung", "condition")], mean)
  out <- aggregate(blk[, features, drop = FALSE], blk[, c("lung", "condition")], mean)
  out <- out[order(out$condition != "control", out$lung), ]
  rownames(out) <- out$lung
  out
}

# Exact two-sided permutation test over all assignments of lungs to groups.
permutation_test <- function(df, features, group = "condition") {
  g <- df[[group]]
  ctrl <- g == "control"
  splits <- combn(nrow(df), sum(!ctrl))
  do.call(rbind, lapply(features, function(f) {
    v <- df[[f]]
    diff <- mean(v[!ctrl]) - mean(v[ctrl])
    null <- apply(splits, 2, function(j) mean(v[j]) - mean(v[-j]))
    data.frame(feature = f, mean_control = mean(v[ctrl]), mean_ipf = mean(v[!ctrl]),
               difference = diff, p_value = mean(abs(null) >= abs(diff) - 1e-12))
  }))
}

# Residualise each feature on log depth (quadratic) and mitochondrial fraction within sections.
adjust_for_depth <- function(obj, features) {
  md <- obj[[]]
  adj <- matrix(NA_real_, nrow(md), length(features),
                dimnames = list(rownames(md), paste0(features, "_adj")))
  for (s in unique(md$section)) {
    ix <- md$section == s
    X <- model.matrix(~ poly(log1p(nCount_Spatial), 2, raw = TRUE) + percent_mt, md[ix, ])
    adj[ix, ] <- qr.resid(qr(X), as.matrix(md[ix, features]))
  }
  Seurat::AddMetaData(obj, as.data.frame(adj))
}

# Moran's I per section on the stored MatriSpace graph, with a permutation P value.
section_morans_i <- function(obj, features, graph, n_perm = 199, seed = 1) {
  g <- matrispace::get_matrispace_graph(obj, graph)
  W <- matrispace::transform_weights(g$connectivities, method = "binary")
  md <- obj[[]][g$spot_ids, ]
  set.seed(seed)
  rows <- lapply(unique(md$section), function(s) {
    ix <- which(md$section == s)
    w <- W[ix, ix]
    z <- scale(as.matrix(md[ix, features]), scale = FALSE)
    moran <- function(z) length(ix) / sum(w) * colSums(z * as.matrix(w %*% z)) / colSums(z^2)
    obs <- moran(z)
    null <- replicate(n_perm, moran(z[sample(nrow(z)), , drop = FALSE]))
    p <- (1 + rowSums(matrix(null, length(features)) >= obs)) / (1 + n_perm)
    names(p) <- paste0(features, "_p")
    data.frame(md[ix[1], design_cols], t(obs), t(p), check.names = FALSE, row.names = NULL)
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------- figure-specific tables

collagen_families <- list(
  Fibrillar = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2"),
  `Beaded (VI)` = c("COL6A1", "COL6A2", "COL6A3"),
  FACIT = c("COL12A1", "COL14A1", "COL16A1", "COL21A1"),
  Network = c("COL8A1", "COL10A1"),
  Multiplexin = c("COL15A1", "COL18A1"),
  `Basement membrane (IV)` = c("COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5")
)

# log2 fold change of mean collagen expression (counts per 10,000) per lung vs the control mean.
collagen_lfc <- function(obj, pseudocount = 0.01) {
  genes <- unlist(collagen_families, use.names = FALSE)
  cp10k <- expm1(as.matrix(Seurat::FetchData(obj, vars = genes)))
  lung <- summarise_by_lung(cbind(obj[[]][, design_cols], cp10k), genes)
  ctrl <- colMeans(lung[lung$condition == "control", genes])
  lfc <- log2(sweep(as.matrix(lung[, genes]) + pseudocount, 2, ctrl + pseudocount, "/"))
  structure(lfc, legend = "log2 fold change", italic = "x")
}

# Spearman correlations between features within each section, summarised per lung and condition.
colocalisation <- function(obj, features) {
  md <- obj[[]]
  pairs <- combn(features, 2)
  pair_names <- paste(pairs[1, ], pairs[2, ], sep = "~")
  per_section <- lapply(unique(md$section), function(s) {
    r <- cor(md[md$section == s, features], method = "spearman")
    data.frame(md[md$section == s, design_cols][1, ],
               t(setNames(r[t(pairs)], pair_names)), check.names = FALSE, row.names = NULL)
  })
  lung <- summarise_by_lung(do.call(rbind, per_section), pair_names)
  reversed <- lung[, pair_names]
  names(reversed) <- paste(pairs[2, ], pairs[1, ], sep = "~")
  cells <- features[startsWith(features, "prop_")]
  to_matrix <- function(v) {
    m <- matrix(NA_real_, length(features), length(features), dimnames = list(features, features))
    m[t(pairs[2:1, ])] <- v
    m <- m[cells, cells]
    m[upper.tri(m, diag = TRUE)] <- NA
    structure(m[-1, -length(cells)], legend = "Spearman's \u03c1", limits = c(-0.5, 0.5))
  }
  list(control = to_matrix(colMeans(lung[lung$condition == "control", pair_names])),
       ipf = to_matrix(colMeans(lung[lung$condition != "control", pair_names])),
       sample = cbind(lung, reversed))
}

# Rank pairs whose lung-level co-expression separates IPF from control.
rank_pairs <- function(obj, min_coexpression = 0.05, max_per_ligand = 2, max_per_receptor = 3,
                       n = 14, pseudocount = 1e-3) {
  enr <- matrispace::find_lr_enrichment(obj, group_by = "section")
  wide <- reshape(enr[, c("pair", "group", "mean_activity")], idvar = "group",
                  timevar = "pair", direction = "wide")
  names(wide) <- sub("^mean_activity\\.", "", names(wide))
  pairs <- setdiff(names(wide), "group")
  wide <- merge(unique(obj[[]][, design_cols]), wide, by.x = "section", by.y = "group")
  lung <- summarise_by_lung(wide, pairs)
  ctrl <- lung[lung$condition == "control", pairs]
  ipf <- lung[lung$condition != "control", pairs]
  stats <- data.frame(pair = pairs, ctrl = colMeans(ctrl), ipf = colMeans(ipf))
  stats$lfc <- log2((stats$ipf + pseudocount) / (stats$ctrl + pseudocount))
  stats$separate <- apply(ipf, 2, min) > apply(ctrl, 2, max) | apply(ipf, 2, max) < apply(ctrl, 2, min)
  stats <- stats[stats$separate & pmax(stats$ctrl, stats$ipf) > min_coexpression & stats$lfc > 0, ]
  stats <- stats[order(-stats$lfc), ]
  lig <- sub("__.*", "", stats$pair)
  rec <- sub(".*__", "", stats$pair)
  chosen <- character(0)
  for (i in seq_len(nrow(stats))) {
    if (sum(lig[match(chosen, stats$pair)] == lig[i]) >= max_per_ligand) next
    if (sum(rec[match(chosen, stats$pair)] == rec[i]) >= max_per_receptor) next
    chosen <- c(chosen, stats$pair[i])
    if (length(chosen) == n) break
  }
  lfc <- log2(sweep(as.matrix(lung[, chosen]) + pseudocount, 2,
                    colMeans(ctrl[, chosen]) + pseudocount, "/"))
  list(pairs = chosen, lfc = structure(t(lfc), legend = "log2 fold change", italic = "y"),
       sample = lung[, c("lung", "condition", chosen)], stats = stats,
       all = lung, min_coexpression = min_coexpression, pseudocount = pseudocount)
}

# Selection-aware test: rerun the pair selection for every relabelling of the lungs and compare
# each selected pair's fold change with the largest fold change of any selected pair.
selection_test <- function(top, pairs) {
  X <- as.matrix(top$all[, setdiff(names(top$all), c("lung", "condition"))])
  X <- X[, colSums(is.na(X)) == 0]
  select <- function(ipf) {
    a <- X[!ipf, , drop = FALSE]; b <- X[ipf, , drop = FALSE]
    sep <- apply(b, 2, min) > apply(a, 2, max)
    expr <- pmax(colMeans(a), colMeans(b)) > top$min_coexpression
    lfc <- log2((colMeans(b) + top$pseudocount) / (colMeans(a) + top$pseudocount))
    list(n = sum(sep & expr), max = max(c(0, lfc[sep & expr])), lfc = lfc)
  }
  ipf <- top$all$condition != "control"
  obs <- select(ipf)
  null <- apply(combn(nrow(X), sum(ipf)), 2, function(j) {
    s <- select(seq_len(nrow(X)) %in% j); c(n = s$n, max = s$max)
  })
  structure(data.frame(feature = pairs, lfc = obs$lfc[pairs],
                       p_value = vapply(obs$lfc[pairs], function(v) mean(null["max", ] >= v), 0)),
            n_separating = obs$n, null_median = median(null["n", ]),
            null_max = sort(null["n", ], decreasing = TRUE)[2])
}

# Copy the co-expression score of one pair into the metadata as LIGAND__RECEPTOR.
add_pair_score <- function(obj, ligand, receptor) {
  act <- matrispace::get_lr_activity(obj)
  id <- matrispace:::.lr_pair_id(ligand, receptor)
  obj[[paste(ligand, receptor, sep = "__")]] <- as.numeric(act[id, colnames(obj)])
  obj
}

# ---------------------------------------------------------------- figure styling

category_labels <- c(collagens = "Collagens", proteoglycans = "Proteoglycans",
                     ecm_glycoproteins = "ECM glycoproteins", ecm_regulators = "ECM regulators",
                     secreted_factors = "Secreted factors", ecm_affiliated_proteins = "ECM-affiliated")
cell_labels <- c(AT1 = "AT1", AT2 = "AT2", Transitional_AT2 = "Transitional AT2",
                 KRT5_KRT17_ = "KRT5−/KRT17+", Myofibroblasts = "Myofibroblasts",
                 HAS1_High_Fibroblasts = "HAS1-high fibroblasts", Fibroblasts = "Fibroblasts",
                 Endothelial_Cells = "Endothelial", Macrophages = "Macrophages",
                 B_Cells = "B cells", Plasma_Cells = "Plasma cells", T_Cells = "T cells")

nice_label <- function(x) {
  x <- sub("_adj$", "", x)
  out <- x
  hit <- x %in% names(category_labels); out[hit] <- category_labels[x[hit]]
  cell <- startsWith(x, "prop_")
  out[cell] <- cell_labels[sub("^prop_", "", x[cell])]
  out <- sub("^HC_", "Control ", out)
  out <- sub("^IPF_", "IPF ", out)
  out <- sub("__", "–", out)
  pair <- grepl("~", x)
  out[pair] <- vapply(strsplit(x[pair], "~"), function(p) {
    p <- nice_label(p); p[p == "Collagens"] <- "Collagen score"; paste(p, collapse = " and ")
  }, "")
  out
}

y_label <- function(features) {
  if (all(grepl("_adj$", features))) return("Moran's I")
  if (all(grepl("~", features))) return("Spearman's \u03c1")
  if (all(grepl("__", features))) return("Co-expression")
  "Expression"
}

theme_figure <- function() {
  theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 10),
          axis.text = element_text(colour = "black"), legend.position = "none")
}

# Box plots of lung-level values by condition, with exact permutation P values.
plot_by_condition <- function(df, features, p = NULL) {
  long <- do.call(rbind, lapply(features, function(f)
    data.frame(feature = f, condition = df$condition, value = df[[f]])))
  long$feature <- factor(long$feature, levels = features, labels = nice_label(features))
  p <- p %||% permutation_test(df, features)$p_value
  top <- tapply(long$value, long$feature, max)
  span <- tapply(long$value, long$feature, function(v) diff(range(v)))
  ann <- data.frame(feature = factor(nice_label(features), levels = nice_label(features)),
                    y = top + 0.12 * span, label = sprintf("P = %.2g", p))
  italic <- all(grepl("__", features)) || all(grepl("^[A-Z0-9]+$", features))
  ggplot(long, aes(condition, value)) +
    geom_boxplot(aes(fill = condition), width = 0.5, coef = Inf,
                 outlier.shape = NA, colour = "black") +
    geom_point(shape = 21, fill = "white", size = 2) +
    geom_segment(data = ann, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE) +
    geom_text(data = ann, aes(x = 1.5, y = y, label = label), vjust = -0.6, size = 3.3,
              fontface = "italic", inherit.aes = FALSE) +
    facet_wrap(~ feature, scales = "free_y", ncol = min(3, length(features)),
               labeller = label_wrap_gen(width = 18)) +
    scale_fill_manual(values = condition_colours) +
    scale_x_discrete(labels = c(control = "Control", IPF = "IPF")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = NULL, y = y_label(features)) +
    theme_figure() +
    theme(strip.text = element_text(face = if (italic) "italic" else "plain"))
}

# MatriSpace heatmap with the figure palette (RdBu), labels and ordering.
plot_heatmap <- function(data, scale_by = "none", center = "auto", limits = NULL, ...) {
  lims <- limits %||% attr(data, "limits") %||% if (scale_by == "row") c(-2, 2) else NULL
  p <- matrispace::plot_heatmap(data, scale_by = scale_by, center = center, limits = lims,
                                na_colour = "white", ...)
  if (is.null(lims)) lims <- c(-1, 1) * max(abs(p$data$value), na.rm = TRUE)
  legend <- attr(data, "legend") %||% if (scale_by == "row") "z-score" else "value"
  italic <- attr(data, "italic") %||% ""
  suppressMessages(p +
    scale_fill_distiller(palette = "RdBu", limits = lims, oob = scales::squish,
                         na.value = "white", name = legend) +
    scale_x_discrete(labels = nice_label) + scale_y_discrete(labels = nice_label, limits = rev) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black",
                                     face = if (italic == "x") "italic" else "plain"),
          axis.text.y = element_text(colour = "black",
                                     face = if (italic == "y") "italic" else "plain"),
          legend.position = "bottom"))
}

# MatriSpace spatial map with the figure palette (reversed Spectral), title and 1 mm scale bar.
plot_spatial <- function(object, feature, limits = NULL, point_size = 2.2, ...) {
  p <- matrispace::plot_spatial(object, feature, limits = limits, point_size = point_size, ...)
  if (is.null(limits)) limits <- c(0, max(p$data$value, na.rm = TRUE))
  legend <- if (startsWith(feature, "prop_")) "Proportion" else
    if (grepl("__", feature)) "Co-expression" else "Signature score"
  x1 <- max(p$data$x); y0 <- min(p$data$y) - 150
  suppressMessages(p +
    scale_colour_distiller(palette = "Spectral", limits = limits, oob = scales::squish,
                           name = legend) +
    annotate("segment", x = x1 - 1000, xend = x1, y = y0, yend = y0, linewidth = 1) +
    annotate("text", x = x1 - 500, y = y0, label = "1 mm", vjust = -0.6, size = 3) +
    labs(title = paste(nice_label(feature), nice_label(unique(object$lung)), sep = ", ")) +
    theme(plot.title = element_text(size = 11, hjust = 0.5),
          legend.position = "bottom"))
}
