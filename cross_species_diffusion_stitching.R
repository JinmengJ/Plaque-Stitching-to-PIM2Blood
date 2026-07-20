# Cross-species macrophage state stitching on an LDLR mouse disease timeline
#
# This script implements the analysis used to project APOE mouse and human
# macrophage/monocyte query cells onto an LDLR mouse atherosclerosis time axis.
# 

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(destiny)
  library(dplyr)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

config <- list(
  input_dir = "/data2/jiajm/Atherosclerosis_Manuscript/figures/AS_analysis_0506",
  output_dir = "/data2/jiajm/Atherosclerosis_Manuscript/figures/AS_analysis_0506/ref",
  ldlr_rds = "Mm_integrated.rds",
  apoe_rds = "Mm.query1_new_maaping.rds",
  human_rds = c(
    coronary_1 = "Mm.query2_new_maaping.rds",
    carotid = "Mm.query3_new_maaping.rds",
    coronary_2 = "Mm.query4_new_maaping.rds"
  ),
  assay = "RNA",
  data_slot = "data",
  ldlr_stage_col = "TimePoints",
  ldlr_celltype_col = "cellgroup",
  query_celltype_col = "predicted.celltype",
  target_celltype = "cluster9",
  hvg_n = 3000L,
  geneset_n = 500L,
  dm_n_pcs = 50L,
  dm_n_dc_max = 20L,
  balanced_cells_per_stage = 500L,
  backend = "knn_balanced",
  knn_k = 8L,
  alpha_fdr = 0.05,
  seed = 1L
)

path_here <- function(input_dir, ...) {
  file.path(input_dir, ...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_seurat <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  readRDS(path)
}

get_expr <- function(object, assay = "RNA", slot = "data") {
  mat <- GetAssayData(object, assay = assay, slot = slot)
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  mat
}

as_dgc <- function(mat) {
  if (inherits(mat, "dgCMatrix")) mat else Matrix::Matrix(as.matrix(mat), sparse = TRUE)
}

fit_reference_z <- function(expr) {
  mu <- Matrix::rowMeans(expr)
  variance <- Matrix::rowSums(expr^2) / ncol(expr) - mu^2
  list(mu = mu, sd = sqrt(pmax(variance, 1e-8)))
}

apply_reference_z <- function(expr, fit) {
  expr <- as_dgc(expr)
  mu <- fit$mu[rownames(expr)]
  sd <- pmax(fit$sd[rownames(expr)], 1e-8)
  inv_sd <- 1 / sd

  for (j in seq_len(ncol(expr))) {
    if (expr@p[j + 1L] > expr@p[j]) {
      idx <- (expr@p[j] + 1L):expr@p[j + 1L]
      rows <- expr@i[idx] + 1L
      expr@x[idx] <- (expr@x[idx] - mu[rows]) * inv_sd[rows]
    }
  }
  expr
}

pairwise_time_variable_genes <- function(object, stage_col, features, p_cut = 0.20,
                                         logfc_cut = 0.20) {
  object <- subset(object, features = intersect(features, rownames(object)))
  stage <- factor(object@meta.data[[stage_col]], levels = sort(unique(object@meta.data[[stage_col]])))
  Idents(object) <- stage

  stage_levels <- levels(stage)
  gene_set <- character(0)

  for (i in seq_along(stage_levels)) {
    for (j in seq_along(stage_levels)) {
      if (i <= j) next

      markers <- tryCatch(
        FindMarkers(
          object,
          ident.1 = stage_levels[i],
          ident.2 = stage_levels[j],
          features = rownames(object),
          test.use = "wilcox",
          min.pct = 0,
          min.cells.group = 1,
          logfc.threshold = 0,
          return.thresh = 1,
          verbose = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(markers) || nrow(markers) == 0) next

      markers$gene <- rownames(markers)
      logfc_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
      p_col <- if ("p_val" %in% colnames(markers)) "p_val" else "p_val_adj"
      keep <- markers[[p_col]] < p_cut & abs(markers[[logfc_col]]) > logfc_cut
      gene_set <- union(gene_set, markers$gene[keep])
    }
  }

  unique(gene_set)
}

select_ldlr_timeline_genes <- function(ldlr, stage_col, assay, slot, hvg_n, target_n,
                                       output_path = NULL) {
  DefaultAssay(ldlr) <- assay
  ldlr <- FindVariableFeatures(ldlr, selection.method = "vst", nfeatures = hvg_n, verbose = FALSE)
  hvg_pool <- intersect(VariableFeatures(ldlr), rownames(ldlr))
  if (length(hvg_pool) < 300) stop("Too few highly variable genes in LDLR reference.")

  tv_genes <- pairwise_time_variable_genes(ldlr, stage_col, hvg_pool)
  genes_use <- intersect(tv_genes, hvg_pool)

  if (length(genes_use) < 300) {
    pool <- setdiff(hvg_pool, genes_use)
    expr <- get_expr(ldlr, assay = assay, slot = slot)
    variance <- apply(expr[pool, , drop = FALSE], 1, var)
    add_genes <- head(names(sort(variance, decreasing = TRUE)), target_n - length(genes_use))
    genes_use <- unique(c(genes_use, add_genes))
  }

  if (!is.null(output_path)) saveRDS(genes_use, output_path)

  list(
    ldlr = ldlr,
    hvg_pool = hvg_pool,
    time_variable_genes = tv_genes,
    genes_use = genes_use
  )
}

estimate_diffusion_sigma <- function(dm, x_ref, knn_k = 50L, sample_max = 2000L) {
  sigma <- tryCatch(methods::slot(dm, "sigma"), error = function(e) NA_real_)
  if (is.numeric(sigma) && is.finite(sigma) && sigma > 0) return(as.numeric(sigma))

  n <- nrow(x_ref)
  idx <- if (n > sample_max) sample.int(n, sample_max) else seq_len(n)
  distances <- as.matrix(dist(x_ref[idx, , drop = FALSE], method = "euclidean"))
  k <- min(knn_k + 1L, nrow(distances))
  knn_dist <- apply(distances, 1, function(v) sort(v, partial = k)[k])
  sigma <- stats::median(knn_dist, na.rm = TRUE)

  if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  as.numeric(sigma)
}

nystrom_project <- function(dm, x_ref, x_new, dc_use) {
  psi_ref <- as.matrix(eigenvectors(dm))[, seq_len(dc_use), drop = FALSE]
  lambda <- as.numeric(eigenvalues(dm))[seq_len(dc_use)]
  sigma <- estimate_diffusion_sigma(dm, x_ref)

  a2 <- rowSums(x_new * x_new)
  b2 <- rowSums(x_ref * x_ref)
  gram <- x_new %*% t(x_ref)
  d2 <- outer(a2, b2, "+") - 2 * gram
  kernel <- exp(-d2 / (sigma^2 + 1e-12))
  transition <- kernel / (rowSums(kernel) + 1e-12)

  as.matrix(sweep(transition %*% psi_ref, 2, lambda, "/"))
}

cosine_similarity <- function(a, b) {
  a_norm <- sqrt(rowSums(a * a)) + 1e-8
  b_norm <- sqrt(rowSums(b * b)) + 1e-8
  sim <- a %*% t(b)
  sim <- sim / a_norm
  sim <- sim / rep(b_norm, each = nrow(a))
  as.matrix(sim)
}

topk_index <- function(x, k) {
  k <- min(k, length(x))
  idx <- integer(k)
  tmp <- x
  for (i in seq_len(k)) {
    idx[i] <- which.max(tmp)
    tmp[idx[i]] <- -Inf
  }
  idx
}

effective_dimension <- function(x) {
  variance <- apply(x, 2, var)
  deff <- (sum(variance)^2) / (sum(variance^2) + 1e-12)
  as.integer(max(3, min(ncol(x), round(deff))))
}

qmax_cosine_null <- function(q, d, m) {
  a <- (d - 1) / 2
  y <- qbeta(q^(1 / m), a, a)
  2 * y - 1
}

pval_max_cosine <- function(r, d, m) {
  a <- (d - 1) / 2
  y <- (pmin(pmax(as.numeric(r), -1), 1) + 1) / 2
  f1 <- pbeta(y, a, a)
  pmax(0, 1 - f1^m)
}

build_ldlr_diffusion_reference <- function(ldlr, genes_use, config) {
  set.seed(config$seed)
  DefaultAssay(ldlr) <- config$assay

  genes_use <- intersect(genes_use, rownames(ldlr))
  if (length(genes_use) < 50) stop("Too few reference genes after intersection with LDLR object.")

  expr <- get_expr(ldlr, assay = config$assay, slot = config$data_slot)
  expr <- expr[genes_use, , drop = FALSE]
  z_fit <- fit_reference_z(expr)
  expr_z <- apply_reference_z(expr, z_fit)
  x_ref <- t(as.matrix(expr_z))

  dm <- destiny::DiffusionMap(x_ref, n_pcs = config$dm_n_pcs)
  dc_all <- as.matrix(eigenvectors(dm))
  dc_keep <- seq_len(min(config$dm_n_dc_max, ncol(dc_all)))
  dc_ref_all <- dc_all[, dc_keep, drop = FALSE]
  rownames(dc_ref_all) <- colnames(ldlr)

  ref_stage <- ldlr@meta.data[[config$ldlr_stage_col]]
  ref_stage <- factor(ref_stage, levels = sort(unique(ref_stage)))
  stage_levels <- levels(ref_stage)

  stage_centroids <- t(sapply(stage_levels, function(stage) {
    colMeans(dc_ref_all[ref_stage == stage, , drop = FALSE])
  }))
  rownames(stage_centroids) <- stage_levels

  balanced_idx <- unlist(lapply(stage_levels, function(stage) {
    cells <- which(ref_stage == stage)
    sample(cells, min(length(cells), config$balanced_cells_per_stage))
  }), use.names = FALSE)

  dc_ref_bal <- dc_ref_all[balanced_idx, , drop = FALSE]
  ref_stage_bal <- ref_stage[balanced_idx]

  list(
    ldlr = ldlr,
    genes_use = genes_use,
    z_fit = z_fit,
    dm = dm,
    x_ref = x_ref,
    dc_ref_all = dc_ref_all,
    ref_stage = ref_stage,
    stage_levels = stage_levels,
    stage_centroids = stage_centroids,
    dc_ref_bal = dc_ref_bal,
    ref_stage_bal = ref_stage_bal
  )
}

calibrate_mapping_parameters <- function(reference, config) {
  dc_ref_all <- reference$dc_ref_all
  dc_ref_bal_all <- reference$dc_ref_bal

  d_eff_all <- effective_dimension(dc_ref_all)
  dc_use <- max(6L, min(ncol(dc_ref_all), min(12L, d_eff_all)))

  dc_ref_bal <- dc_ref_bal_all[, seq_len(dc_use), drop = FALSE]
  d_eff <- effective_dimension(dc_ref_all[, seq_len(dc_use), drop = FALSE])
  m_ref <- nrow(dc_ref_bal)

  sim <- cosine_similarity(dc_ref_bal, dc_ref_bal)
  diag(sim) <- -Inf

  purity <- numeric(nrow(sim))
  gap <- numeric(nrow(sim))
  for (i in seq_len(nrow(sim))) {
    idx <- topk_index(sim[i, ], config$knn_k)
    labels <- reference$ref_stage_bal[idx]
    tab <- table(labels)
    purity[i] <- max(tab) / config$knn_k
    top_values <- sort(sim[i, idx], decreasing = TRUE)
    gap[i] <- top_values[1] - ifelse(length(top_values) >= 2, top_values[2], 0)
  }

  list(
    backend = config$backend,
    knn_k = config$knn_k,
    alpha_fdr = config$alpha_fdr,
    dc_use = dc_use,
    d_eff = d_eff,
    m_ref = m_ref,
    m_ref_sig = m_ref,
    thr_top1_global = qmax_cosine_null(0.90, d_eff, m_ref),
    purity_map_min = 0.50,
    purity_thr = min(max(as.numeric(quantile(purity, 0.70, na.rm = TRUE)), 0.60), 0.95),
    gap_thr = min(max(as.numeric(quantile(gap, 0.70, na.rm = TRUE)), 0.000), 0.050)
  )
}

project_query_to_reference <- function(query, reference, params, config) {
  DefaultAssay(query) <- config$assay
  expr <- get_expr(query, assay = config$assay, slot = config$data_slot)
  feats <- intersect(intersect(rownames(expr), reference$genes_use), names(reference$z_fit$mu))
  if (length(feats) < 200) stop("Too few query genes overlap the locked LDLR feature set.")

  expr_z <- apply_reference_z(
    expr[feats, , drop = FALSE],
    list(mu = reference$z_fit$mu[feats], sd = reference$z_fit$sd[feats])
  )

  dc_query <- nystrom_project(
    dm = reference$dm,
    x_ref = reference$x_ref[, feats, drop = FALSE],
    x_new = t(as.matrix(expr_z)),
    dc_use = params$dc_use
  )
  rownames(dc_query) <- colnames(query)
  dc_query
}

map_query_stage <- function(query, reference, params, config, query_name = "query") {
  dc_query <- project_query_to_reference(query, reference, params, config)
  dc_ref <- reference$dc_ref_bal[, seq_len(params$dc_use), drop = FALSE]

  sim <- cosine_similarity(dc_query[, seq_len(params$dc_use), drop = FALSE], dc_ref)

  pred <- character(nrow(sim))
  purity <- numeric(nrow(sim))
  top1 <- numeric(nrow(sim))
  gap12 <- numeric(nrow(sim))

  for (i in seq_len(nrow(sim))) {
    v <- sim[i, ]
    idx <- topk_index(v, params$knn_k)
    labels <- reference$ref_stage_bal[idx]
    tab <- table(labels)
    pred[i] <- names(which.max(tab))
    purity[i] <- max(tab) / params$knn_k
    ordered <- sort(v, decreasing = TRUE)
    top1[i] <- ordered[1]
    gap12[i] <- ordered[1] - ifelse(length(ordered) >= 2, ordered[2], 0)
  }

  p_raw <- pval_max_cosine(top1, d = params$d_eff, m = params$m_ref_sig)
  fdr <- p.adjust(p_raw, method = "BH")

  is_mapped <- top1 >= params$thr_top1_global & purity >= params$purity_map_min
  mapped_stage <- ifelse(is_mapped, pred, "Query-late")

  stage_category <- rep("Late-specific", length(is_mapped))
  mapped_idx <- which(is_mapped)
  if (length(mapped_idx) > 0) {
    is_sig <- fdr[mapped_idx] < params$alpha_fdr
    is_concordant <- purity[mapped_idx] >= params$purity_thr
    is_gap <- gap12[mapped_idx] >= params$gap_thr
    stage_category[mapped_idx] <- ifelse(
      is_sig & is_concordant & is_gap,
      "Stage-specific",
      ifelse(is_sig, "Pan-stage", "Intermediate")
    )
  }

  names(pred) <- names(mapped_stage) <- names(stage_category) <- rownames(sim)
  names(top1) <- names(gap12) <- names(purity) <- names(p_raw) <- names(fdr) <- rownames(sim)

  query$ProjectionName <- query_name
  query$PredictedStage <- pred[colnames(query)]
  query$MappedStage <- mapped_stage[colnames(query)]
  query$StageCategory <- stage_category[colnames(query)]
  query$SimilarityTop1 <- top1[colnames(query)]
  query$SimilarityGap12 <- gap12[colnames(query)]
  query$NeighborPurityK <- purity[colnames(query)]
  query$SimilarityP <- p_raw[colnames(query)]
  query$SimilarityFDR <- fdr[colnames(query)]

  list(query = query, dc_query = dc_query, similarity = sim)
}

run_seurat_label_transfer <- function(reference, query, config) {
  DefaultAssay(reference) <- "integrated"
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    dims = 1:30,
    reference.reduction = "pca"
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference@meta.data[[config$ldlr_celltype_col]],
    dims = 1:30
  )
  query <- AddMetaData(query, metadata = predictions)
  query <- MapQuery(
    anchorset = anchors,
    reference = reference,
    query = query,
    refdata = list(celltype = config$ldlr_celltype_col),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
  query
}

summarise_target_celltype_by_stage <- function(object, celltype_col, stage_col, target_celltype) {
  object@meta.data %>%
    mutate(.is_target = .data[[celltype_col]] == target_celltype) %>%
    group_by(stage = .data[[stage_col]]) %>%
    summarise(
      n_cells = n(),
      n_target = sum(.is_target, na.rm = TRUE),
      target_fraction = n_target / n_cells,
      .groups = "drop"
    )
}

write_mapping_outputs <- function(mapped, name, output_dir) {
  ensure_dir(output_dir)
  saveRDS(mapped$query, file.path(output_dir, paste0(name, "_mapped.rds")))
  write.csv(
    mapped$query@meta.data,
    file.path(output_dir, paste0(name, "_metadata.csv")),
    quote = TRUE
  )
  write.csv(
    mapped$dc_query,
    file.path(output_dir, paste0(name, "_diffusion_coordinates.csv")),
    quote = FALSE
  )
}

main <- function(config) {
  set.seed(config$seed)
  ensure_dir(config$output_dir)

  ldlr <- read_seurat(path_here(config$input_dir, config$ldlr_rds))

  genes_path <- file.path(config$output_dir, "genes_use.rds")
  gene_selection <- select_ldlr_timeline_genes(
    ldlr = ldlr,
    stage_col = config$ldlr_stage_col,
    assay = config$assay,
    slot = config$data_slot,
    hvg_n = config$hvg_n,
    target_n = config$geneset_n,
    output_path = genes_path
  )

  reference <- build_ldlr_diffusion_reference(
    ldlr = gene_selection$ldlr,
    genes_use = gene_selection$genes_use,
    config = config
  )
  params <- calibrate_mapping_parameters(reference, config)

  saveRDS(reference, file.path(config$output_dir, "ldlr_diffusion_reference.rds"))
  saveRDS(params, file.path(config$output_dir, "mapping_parameters.rds"))

  apoe <- read_seurat(path_here(config$input_dir, config$apoe_rds))
  apoe_mapped <- map_query_stage(apoe, reference, params, config, query_name = "APOE")
  write_mapping_outputs(apoe_mapped, "apoe", config$output_dir)

  human_mapped <- lapply(names(config$human_rds), function(name) {
    query <- read_seurat(path_here(config$input_dir, config$human_rds[[name]]))
    mapped <- map_query_stage(query, reference, params, config, query_name = name)
    write_mapping_outputs(mapped, name, config$output_dir)
    mapped$query
  })
  names(human_mapped) <- names(config$human_rds)

  human_all <- Reduce(function(x, y) merge(x, y), human_mapped)
  saveRDS(human_all, file.path(config$output_dir, "human_queries_mapped_merged.rds"))
  write.csv(human_all@meta.data, file.path(config$output_dir, "human_queries_mapped_metadata.csv"))

  target_summary <- summarise_target_celltype_by_stage(
    human_all,
    celltype_col = config$query_celltype_col,
    stage_col = "MappedStage",
    target_celltype = config$target_celltype
  )
  write.csv(target_summary, file.path(config$output_dir, "target_celltype_by_mapped_stage.csv"), row.names = FALSE)

  invisible(list(
    reference = reference,
    params = params,
    apoe = apoe_mapped$query,
    human = human_all,
    target_summary = target_summary
  ))
}

if (identical(environment(), globalenv()) && !interactive()) {
  main(config)
}
