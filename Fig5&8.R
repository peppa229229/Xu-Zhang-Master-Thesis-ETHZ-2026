## =========================================================
## TUMOR (tumor_sc==TRUE) — PATIENT FILTER (Pt108) + FUNCTIONAL-STATE-ONLY FILTER
## + RE-CLUSTER (PCA/Neighbors/FindClusters/UMAP) ON FILTERED CELLS
## + FIX: cluster ordering numeric (c1,c2,c3,...,c10...) everywhere
## + FIX: drop weird prefixes like cg/c/C/cluster_ -> ALWAYS normalize to integer
## + FIX: one typo bug in UMAP section (mdx -> um2 already resolved)
## + INCLUDE: NeoTCR scoring + UMAP projection + clonotype ranking
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(openxlsx)
  library(ggrepel)
  library(scales)
  library(grid)
  library(Matrix)
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Missing package: patchwork. Run install.packages('patchwork') and retry.")
  }
  library(patchwork)
  
  if (!requireNamespace("ragg", quietly = TRUE)) {
    message("[warn] ragg not installed; PNG export will fallback to default device.")
  }
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    message("[warn] GSVA not installed; ssGSEA heatmap section will fail unless you install it.")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("[warn] ComplexHeatmap not installed; heatmap sections will fail unless you install it.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    message("[warn] circlize not installed; heatmap sections will fail unless you install it.")
  }
})

## =========================================================
## User-tunable parameters
## =========================================================
deep_red <- "#8B0000"

label_every <- 1L
state_max_chars <- 14
label_font_size <- 1.9
heat_label_heights <- c(9.2, 1.25)

png_dpi <- 600
options(ragg.max_dim = 100000)

npcs_run <- 50
dims_use <- 1:30
var_nfeatures <- 2000
umap_name <- "umap"

cluster_resolution <- 0.6

expr_layer <- "data"
min_present_genes <- 5
ssgsea_norm <- TRUE
drop_const_genes <- TRUE

match_key_for_tumor <- "vjaa"        # "vjaa" or "CDR3b"
match_mode_for_tumor <- "contains"   # "contains" or "exact"

## =========================================================
## Patient filter (Pt108)
## =========================================================
patient_col <- "Patient"
patient_keep <- "Pt108"

## =========================================================
## 1) Paths
## =========================================================
xlsx_path <- "C:/Users/peppa/Downloads/sc_data_availibility.xlsx"
rds_path  <- "C:/Users/peppa/Downloads/MEL_scRNA_TCR_objects/Melanoma_TUMOR_scRNA_072324_scRep_TCRs_Reactivity_Categories_scCloneIDs.Rds"
sig_merged_xlsx <- "C:/Users/peppa/Downloads/signature_xu_merged.xlsx"

stopifnot(file.exists(xlsx_path), file.exists(rds_path))
if (!file.exists(sig_merged_xlsx)) {
  message("[warn] signature_xu_merged.xlsx not found; ssGSEA/signature heatmaps will fail unless you set sig_merged_xlsx correctly.")
}

out_dir <- "C:/Users/peppa/Downloads/"
out_prefix <- file.path(
  out_dir,
  paste0("tumor_Pt108_funStateOnly_RECLUSTER_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1.5) Save helpers
## =========================================================
save_plot_best <- function(filename_base, plot, width, height, dpi = 1200) {
  tryCatch({
    ggsave(
      paste0(filename_base, ".pdf"),
      plot, width = width, height = height, device = grDevices::cairo_pdf
    )
  }, error = function(e) {
    message("[save] cairo_pdf failed, fallback to pdf(): ", conditionMessage(e))
    ggsave(
      paste0(filename_base, ".pdf"),
      plot, width = width, height = height, device = grDevices::pdf
    )
  })
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      paste0(filename_base, ".png"),
      plot, width = width, height = height, device = ragg::agg_png, dpi = dpi
    )
  } else {
    tryCatch({
      ggsave(
        paste0(filename_base, ".png"),
        plot, width = width, height = height, device = grDevices::cairo_png, dpi = dpi
      )
    }, error = function(e) {
      message("[save] cairo_png failed, fallback to default png(): ", conditionMessage(e))
      ggsave(
        paste0(filename_base, ".png"),
        plot, width = width, height = height, dpi = dpi
      )
    })
  }
  
  tryCatch({
    ggsave(
      paste0(filename_base, ".svg"),
      plot, width = width, height = height, device = grDevices::svg
    )
  }, error = function(e) {
    message("[save] svg failed (skip): ", conditionMessage(e))
  })
}

save_ht_best <- function(filename_base, ht, width, height, dpi = 600) {
  stopifnot(requireNamespace("ComplexHeatmap", quietly = TRUE))
  
  tryCatch({
    grDevices::cairo_pdf(paste0(filename_base, ".pdf"), width = width, height = height)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    grDevices::dev.off()
  }, error = function(e) {
    grDevices::pdf(paste0(filename_base, ".pdf"), width = width, height = height)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    grDevices::dev.off()
  })
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      paste0(filename_base, ".png"),
      width = width, height = height, units = "in", res = dpi
    )
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    grDevices::dev.off()
  } else {
    grDevices::png(
      paste0(filename_base, ".png"),
      width = width, height = height, units = "in", res = dpi
    )
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    grDevices::dev.off()
  }
}

save_plot_best_local <- function(filename_base, plot, width, height, dpi = 600) {
  if (exists("save_plot_best", mode = "function")) {
    save_plot_best(filename_base, plot, width = width, height = height, dpi = dpi)
  } else {
    ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height)
    ggsave(paste0(filename_base, ".png"), plot, width = width, height = height, dpi = dpi)
  }
}

## =========================================================
## 1.6) Cluster helpers
## =========================================================
normalize_cluster_id <- function(x, keep_na_last = TRUE) {
  x0 <- as.character(x)
  x0[is.na(x0) | trimws(x0) == ""] <- "NA"
  
  num <- suppressWarnings(as.integer(stringr::str_extract(x0, "\\d+")))
  id <- ifelse(!is.na(num) & x0 != "NA", as.character(num), "NA")
  
  uniq <- unique(id)
  uniq <- uniq[!is.na(uniq) & trimws(uniq) != ""]
  nums <- suppressWarnings(as.integer(uniq))
  num_levels <- sort(unique(nums[!is.na(nums)]))
  lev <- as.character(num_levels)
  
  if ("NA" %in% uniq) {
    if (keep_na_last) {
      lev <- c(lev, "NA")
    } else {
      lev <- c("NA", lev)
    }
  }
  
  list(id = id, num = num, levels = lev)
}

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  x2 <- toupper(trimws(as.character(x)))
  x2 %in% c("TRUE","T","1","YES","Y")
}

GetLayerDataSafe <- function(seu, assay = "RNA", layer = "data") {
  DefaultAssay(seu) <- assay
  GetAssayData(seu, assay = assay, layer = layer)
}

mode_chr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

short_state <- function(x, n = state_max_chars) {
  x <- ifelse(is.na(x) | trimws(x) == "", "Other", as.character(x))
  stringr::str_trunc(x, width = n, side = "right")
}

pick_first_existing <- function(cands, pool) {
  hit <- cands[cands %in% pool]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

clean_key <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x
}

canon_key <- function(x) {
  x <- toupper(as.character(x))
  gsub("[^A-Z0-9]+", "", x)
}

find_feature <- function(seu, assay, candidates) {
  if (!(assay %in% Assays(seu))) return(NA_character_)
  feats <- rownames(seu[[assay]])
  if (length(feats) == 0) return(NA_character_)
  fk <- canon_key(feats)
  ck <- canon_key(candidates)
  
  idx <- match(ck, fk)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(feats[idx[1]])
  
  for (k in ck) {
    hit <- which(grepl(k, fk, fixed = TRUE))
    if (length(hit) > 0) return(feats[hit[1]])
  }
  NA_character_
}

clean_gene_vec <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "NONE", "0")] <- NA_character_
  unique(x[!is.na(x)])
}

map_genes_to_features <- function(seu, genes_upper, assay = "RNA") {
  feats <- rownames(seu[[assay]])
  feats_u <- toupper(feats)
  idx <- match(genes_upper, feats_u)
  feats[idx[!is.na(idx)]]
}

get_gene_vec <- function(seu, gene) {
  feats <- rownames(seu[["RNA"]])
  fu <- toupper(feats)
  idx <- match(toupper(gene), fu)
  if (is.na(idx)) return(rep(0, ncol(seu)))
  vdf <- tryCatch(FetchData(seu, vars = feats[idx], layer = "data"), error = function(e) NULL)
  if (is.null(vdf)) vdf <- FetchData(seu, vars = feats[idx])
  v <- as.numeric(vdf[, 1])
  names(v) <- rownames(vdf)
  v[colnames(seu)]
}

get_seurat <- function(x) {
  if (inherits(x, "Seurat")) return(x)
  if (is.list(x)) {
    idx <- which(sapply(x, inherits, "Seurat"))
    if (length(idx) >= 1) return(x[[idx[1]]])
  }
  stop("No Seurat object found inside the RDS.")
}

resolve_to_rna_gene <- function(key) {
  map <- list(
    "CD4"   = c("CD4"),
    "CD8a"  = c("CD8A","CD8"),
    "IL7Aa" = c("IL7R","IL7RA"),
    "CD28"  = c("CD28"),
    "CCR7"  = c("CCR7"),
    "CD62L" = c("SELL"),
    "CD95"  = c("FAS"),
    "CD27"  = c("CD27"),
    "TCF7"  = c("TCF7"),
    "PD1"   = c("PDCD1"),
    "TIGIT" = c("TIGIT"),
    "CD39"  = c("ENTPD1"),
    "ICOS"  = c("ICOS"),
    "CD25"  = c("IL2RA"),
    "CD137" = c("TNFRSF9"),
    "PRF1"  = c("PRF1"),
    "IFNG"  = c("IFNG"),
    "GZMB"  = c("GZMB"),
    "GZMK"  = c("GZMK"),
    "FASLG" = c("FASLG"),
    "GZMH"  = c("GZMH"),
    "GZMM"  = c("GZMM"),
    "MKI67" = c("MKI67"),
    "TOX"   = c("TOX"),
    "FOXP3" = c("FOXP3"),
    "PRDM1" = c("PRDM1"),
    "NR4A1" = c("NR4A1")
  )
  map[[key]]
}

## =========================================================
## 2) Read Excel and filter tumor_sc == TRUE
## =========================================================
dat <- readxl::read_xlsx(xlsx_path, .name_repair = "minimal")

need_cols <- c(
  "CDR3b","TCR","guide_id","vjaa","Beta_clonotype","PD","tumor_vjaa_counts",
  "Post3rd_skin","Pre3rd_skin","skin_vjaa_counts","Postvax_blood_bulk_TCR",
  "ifng_vjaa_counts_blood_scTCR","total_vjaa_counts","IMP-expanded","tumor_sc",
  "skin_sc","Postvax_blood","in_vitro_exp_sc_CR","VAL"
)
missing_cols <- setdiff(need_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("Excel is missing columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat %>%
  mutate(
    tumor_sc = to_logical(tumor_sc),
    VAL = to_logical(VAL)
  )

tumor_tab <- dat %>%
  filter(tumor_sc %in% TRUE) %>%
  mutate(
    TCR = as.character(TCR),
    vjaa = as.character(vjaa),
    CDR3b = as.character(CDR3b)
  ) %>%
  filter(!is.na(TCR), TCR != "")

if (nrow(tumor_tab) == 0) stop("No rows found where tumor_sc == TRUE.")

tumor_tab <- tumor_tab %>%
  mutate(
    tcr_num = str_extract(TCR, "\\d+"),
    TCR_label = ifelse(!is.na(tcr_num), paste0("TCR", tcr_num), paste0("TCR", TCR))
  )

## =========================================================
## 3) Load Seurat + patient filter + functional-state filter + recluster + UMAP
## =========================================================
obj <- readRDS(rds_path)
seu0 <- get_seurat(obj)

DefaultAssay(seu0) <- "RNA"

if (!(patient_col %in% colnames(seu0@meta.data))) {
  stop(
    "meta.data missing patient column '", patient_col, "'. Available: ",
    paste(colnames(seu0@meta.data), collapse = ", ")
  )
}

pat_vec <- as.character(seu0@meta.data[[patient_col]])
pat_vec[is.na(pat_vec)] <- ""
pat_target <- toupper(trimws(as.character(patient_keep)))

pat_keep_mask <- toupper(trimws(pat_vec)) %in% c(
  pat_target,
  gsub("^PT", "P", pat_target),
  gsub("^P", "PT", pat_target)
)

patient_cells <- rownames(seu0@meta.data)[pat_keep_mask]
if (length(patient_cells) == 0) {
  stop(
    "No cells matched patient=", patient_keep, " using meta column ", patient_col,
    ". Unique Patient values (first 30): ",
    paste(head(unique(pat_vec), 30), collapse = ", ")
  )
}

seu0p <- subset(seu0, cells = patient_cells)
message("[info] patient filter kept cells: ", ncol(seu0p), " | Patient=", patient_keep, " | col=", patient_col)

state_col <- "functional.cluster"
if (!(state_col %in% colnames(seu0p@meta.data))) {
  stop("meta.data missing functional.cluster; cannot filter.")
}

st_vec <- as.character(seu0p@meta.data[[state_col]])
keep_cells <- rownames(seu0p@meta.data)[!is.na(st_vec) & trimws(st_vec) != ""]
if (length(keep_cells) == 0) {
  stop("After filtering by functional state (within patient), 0 cells remain.")
}

seu <- subset(seu0p, cells = keep_cells)
message("[info] kept cells with functional state (within patient): ", ncol(seu))

seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
seu <- FindVariableFeatures(seu, assay = "RNA", nfeatures = var_nfeatures, verbose = FALSE)

if (!is.null(seu[["RNA"]]@layers) && "scale.data" %in% names(seu[["RNA"]]@layers)) {
  seu[["RNA"]]@layers[["scale.data"]] <- NULL
}

seu <- ScaleData(seu, assay = "RNA", features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, assay = "RNA", features = VariableFeatures(seu), npcs = npcs_run, verbose = FALSE)

seu <- FindNeighbors(seu, reduction = "pca", dims = dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = dims_use, reduction.name = umap_name, verbose = FALSE)

stopifnot(umap_name %in% names(seu@reductions))
stopifnot("seurat_clusters" %in% colnames(seu@meta.data))

cln <- normalize_cluster_id(seu$seurat_clusters, keep_na_last = TRUE)
seu$seurat_clusters <- factor(cln$id, levels = cln$levels)
seu$seurat_clusters <- droplevels(seu$seurat_clusters)

seu_rds_out <- paste0(
  out_prefix,
  "_Seurat_PATIENT_",
  gsub("[^A-Za-z0-9]+","_",patient_keep),
  "_filtered_funStateOnly_RECLUSTERED.rds"
)
saveRDS(seu, seu_rds_out)
message("[saved] ", seu_rds_out)

## =========================================================
## 4) Marker modules
## =========================================================
genes_tcell  <- c("CD3E","CD4","CD8A","CD8B")
genes_memory <- c("CCR7","TCF7","SELL","CD28","CD27","IL7R","FAS","LEF1")
genes_effector <- c("NKG7","CST7","FASLG","PRF1","GZMA","GZMB","GZMH","GZMK","GZMM","GNLY",
                    "IL21","IL2","TNF","IFNG","CCL4","CCL5")
genes_pheno <- c("CD38","ENTPD1","ITGAE","CD69","KLRG1","CD40LG","TNFRSF4","TNFRSF9","ICOS",
                 "ITGAL","ITGA1","ITGB1","CX3CR1","CXCR6","IL2RA")
genes_inhib <- c("PDCD1","CTLA4","TIGIT","LAG3","HAVCR2","CD244","CD160","BTLA","VTCN1","TOX")
genes_tf <- c("MKI67","TBX21","EOMES","ZNF683","ID2","ID3","PRDM1","GATA3","FOXP3")

genes_all <- c(genes_tcell, genes_memory, genes_effector, genes_pheno, genes_inhib, genes_tf)
gene_group <- c(
  rep("T cell", length(genes_tcell)),
  rep("Memory/naive", length(genes_memory)),
  rep("Effector function", length(genes_effector)),
  rep("Phenotype", length(genes_pheno)),
  rep("Inhibitory/exhaustion", length(genes_inhib)),
  rep("Transcription factors", length(genes_tf))
)
group_levels <- c(
  "T cell","Memory/naive","Effector function","Phenotype",
  "Inhibitory/exhaustion","Transcription factors"
)

feat_names <- rownames(seu[["RNA"]])
present_mask <- !is.na(match(genes_all, feat_names))
present_genes  <- genes_all[present_mask]
present_groups <- gene_group[present_mask]
if (length(present_genes) < 10) {
  stop("Too few marker genes found in the RNA assay after filtering.")
}

## =========================================================
## 5) Match tumor_sc rows -> Seurat cells
## =========================================================
meta <- seu@meta.data
stopifnot("seurat_clusters" %in% colnames(meta))

if (!(match_key_for_tumor %in% colnames(tumor_tab))) {
  stop("tumor_tab missing key column: ", match_key_for_tumor)
}
if (!(match_key_for_tumor %in% colnames(meta))) {
  stop("Seurat meta missing key column: ", match_key_for_tumor)
}

meta_key <- as.character(meta[[match_key_for_tumor]])
meta_key[is.na(meta_key)] <- ""

tumor_tab <- tumor_tab %>%
  mutate(
    row_id = row_number(),
    TCR_unique = make.unique(TCR_label, sep = "_dup")
  )

map_list <- lapply(seq_len(nrow(tumor_tab)), function(i) {
  keyi <- as.character(tumor_tab[[match_key_for_tumor]][i])
  if (is.na(keyi) || trimws(keyi) == "") return(NULL)
  
  idx <- if (match_mode_for_tumor == "exact") {
    which(meta_key == keyi)
  } else {
    which(grepl(keyi, meta_key, fixed = TRUE))
  }
  
  cells <- rownames(meta)[idx]
  rep_cell <- if (length(cells) >= 1) cells[1] else NA_character_
  
  cl_raw <- if (length(cells) >= 1) {
    mode_chr(as.character(meta[cells, "seurat_clusters"]))
  } else {
    NA_character_
  }
  cl_norm <- normalize_cluster_id(cl_raw, keep_na_last = TRUE)$id[1]
  
  st <- if (length(cells) >= 1 && (state_col %in% colnames(meta))) {
    mode_chr(as.character(meta[cells, state_col]))
  } else {
    NA_character_
  }
  
  list(
    row_id = tumor_tab$row_id[i],
    TCR_unique = tumor_tab$TCR_unique[i],
    TCR_label  = tumor_tab$TCR_label[i],
    TCR_excel  = as.character(tumor_tab$TCR[i]),
    VAL        = as.logical(tumor_tab$VAL[i]),
    key        = keyi,
    n_cells    = length(cells),
    cells      = cells,
    rep_cell   = rep_cell,
    cluster    = cl_norm,
    state      = st
  )
})
map_list <- Filter(Negate(is.null), map_list)

map_df <- bind_rows(lapply(map_list, function(x) {
  tibble(
    row_id     = x$row_id,
    TCR_unique = x$TCR_unique,
    TCR_label  = x$TCR_label,
    TCR_excel  = x$TCR_excel,
    VAL        = x$VAL,
    key        = x$key,
    n_cells    = x$n_cells,
    rep_cell   = x$rep_cell,
    cluster    = x$cluster,
    state      = x$state
  )
}))

missing_df <- map_df %>% filter(is.na(rep_cell) | n_cells == 0)
keep_map   <- map_df %>% filter(!is.na(rep_cell), n_cells > 0)
keep_list  <- map_list[sapply(map_list, function(x) !is.na(x$rep_cell) && x$n_cells > 0)]
if (nrow(keep_map) == 0) {
  stop("No tumor_sc matches found in filtered Seurat using key=", match_key_for_tumor)
}

expr_layer_mat <- GetLayerDataSafe(seu, assay = "RNA", layer = "data")
expr_layer_mat <- expr_layer_mat[present_genes, , drop = FALSE]

expr_mat <- matrix(NA_real_, nrow = length(present_genes), ncol = nrow(keep_map))
rownames(expr_mat) <- present_genes
colnames(expr_mat) <- keep_map$TCR_unique

for (k in seq_along(keep_list)) {
  cells <- keep_list[[k]]$cells
  cells <- intersect(cells, colnames(expr_layer_mat))
  if (length(cells) == 0) next
  
  X <- expr_layer_mat[, cells, drop = FALSE]
  avg <- Matrix::rowMeans(X)
  expr_mat[, k] <- pmin(pmax(as.numeric(avg), 0), 2)
}

## =========================================================
## 6) Main table + module means
## =========================================================
module_mean <- function(gset) {
  gset <- intersect(gset, rownames(expr_mat))
  if (length(gset) == 0) return(rep(NA_real_, ncol(expr_mat)))
  colMeans(expr_mat[gset, , drop = FALSE])
}

mod_df <- tibble(
  TCR_unique    = colnames(expr_mat),
  mean_Tcell    = module_mean(genes_tcell),
  mean_Memory   = module_mean(genes_memory),
  mean_Effector = module_mean(genes_effector),
  mean_Pheno    = module_mean(genes_pheno),
  mean_Inhib    = module_mean(genes_inhib),
  mean_TF       = module_mean(genes_tf)
)

big_table <- keep_map %>%
  left_join(mod_df, by = "TCR_unique") %>%
  mutate(
    TCR_index = row_number(),
    cluster_clean = ifelse(is.na(cluster) | trimws(cluster) == "", "NA", as.character(cluster)),
    cluster_label = paste0("c", cluster_clean),
    state_clean   = short_state(state)
  )

## =========================================================
## 7) Marker-module heatmap + bottom labels
## =========================================================
x_df <- big_table %>%
  transmute(
    TCR_unique,
    TCR_index,
    TCR_label,
    TCR_excel,
    VAL,
    cluster_label,
    state_clean
  )

n_cols <- nrow(x_df)
show_label <- rep(TRUE, n_cols)

df_long <- as.data.frame(expr_mat) %>%
  mutate(gene = rownames(expr_mat)) %>%
  pivot_longer(-gene, names_to = "TCR_unique", values_to = "value") %>%
  left_join(x_df, by = "TCR_unique") %>%
  mutate(
    group = present_groups[match(gene, present_genes)],
    group = factor(group, levels = group_levels),
    gene  = factor(gene, levels = rev(present_genes)),
    TCR_unique = factor(TCR_unique, levels = x_df$TCR_unique)
  )

vline_df <- if (n_cols >= 2) {
  tibble(x = seq(1.5, n_cols - 0.5, by = 1))
} else {
  tibble(x = numeric(0))
}

rect_df <- x_df %>%
  filter(VAL %in% TRUE) %>%
  transmute(
    xmin = TCR_index - 0.5,
    xmax = TCR_index + 0.5,
    ymin = -Inf,
    ymax = Inf
  )

p_expr <- ggplot(df_long, aes(x = TCR_unique, y = gene, fill = value)) +
  geom_tile() +
  { if (nrow(vline_df) > 0)
    geom_vline(
      data = vline_df, aes(xintercept = x), inherit.aes = FALSE,
      color = "white", linewidth = 0.25
    )
  } +
  { if (nrow(rect_df) > 0)
    geom_rect(
      data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = NA, color = "black",
      linetype = "dashed", linewidth = 0.9
    )
  } +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low = "lightgrey", high = "red", limits = c(0, 2), breaks = c(0, 1, 2)) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0), limits = x_df$TCR_unique) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = paste0("Tumor (tumor_sc == TRUE) | Patient=", patient_keep, " | Marker-module heatmap (FILTERED + RE-CLUSTERED)"),
    subtitle = paste0("Bottom labels: State / NewCluster(c#) / TCR | match_key=", match_key_for_tumor),
    x = NULL, y = NULL, fill = "Expr (clamped)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.spacing.y = unit(0.10, "lines"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(4, 6, 0, 6)
  )

label_panel_df <- x_df %>%
  mutate(
    show = show_label,
    label_color = ifelse(VAL %in% TRUE, deep_red, "black"),
    TCR_unique = factor(TCR_unique, levels = x_df$TCR_unique)
  )

df_lab <- bind_rows(
  label_panel_df %>% transmute(TCR_unique, row = "State",   txt = state_clean,   col = label_color, show),
  label_panel_df %>% transmute(TCR_unique, row = "Cluster", txt = cluster_label, col = label_color, show),
  label_panel_df %>% transmute(TCR_unique, row = "TCR",     txt = TCR_label,     col = label_color, show)
) %>%
  filter(show) %>%
  mutate(row = factor(row, levels = c("State","Cluster","TCR")))

p_labels <- ggplot() +
  geom_text(
    data = df_lab,
    aes(x = TCR_unique, y = row, label = txt, color = col),
    size = label_font_size, fontface = "bold", lineheight = 0.95
  ) +
  scale_color_identity() +
  scale_x_discrete(drop = FALSE, expand = c(0, 0), limits = x_df$TCR_unique) +
  scale_y_discrete(drop = FALSE, expand = expansion(mult = c(0.10, 0.10))) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 11) +
  theme(plot.margin = margin(0, 10, 2, 10))

p_heat <- p_expr / p_labels + plot_layout(heights = heat_label_heights)

heat_w <- max(14, 0.34 * n_cols + 4)
heat_h <- max(10, 0.12 * length(present_genes) + 6)
save_plot_best(
  paste0(out_prefix, "_heatmap_markerPartitions_labelsAll_NEWclusters"),
  p_heat, width = heat_w, height = heat_h, dpi = png_dpi
)

## =========================================================
## 8) UMAP #1: label all selected TCRs; VAL red
## =========================================================
um <- as.data.frame(Embeddings(seu, reduction = umap_name))[, 1:2, drop = FALSE]
colnames(um) <- c("UMAP_1","UMAP_2")
um$cell <- rownames(um)

sel_all <- big_table %>%
  select(TCR_unique, TCR_excel, VAL, rep_cell) %>%
  left_join(um, by = c("rep_cell" = "cell")) %>%
  filter(!is.na(UMAP_1), !is.na(UMAP_2)) %>%
  mutate(
    tcr_num     = str_extract(as.character(TCR_excel), "\\d+"),
    TCR_label2  = ifelse(!is.na(tcr_num), paste0("TCR", tcr_num), paste0("TCR", as.character(TCR_excel))),
    label_color = ifelse(VAL %in% TRUE, deep_red, "black")
  )

sel_val_um1 <- sel_all %>% filter(VAL %in% TRUE)

p_umap_selected <- ggplot(um, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey85", size = 0.22) +
  geom_point(data = sel_all, aes(UMAP_1, UMAP_2), color = "black", size = 1.4) +
  { if (nrow(sel_val_um1) > 0)
    geom_point(data = sel_val_um1, aes(UMAP_1, UMAP_2), color = deep_red, size = 2.6)
  } +
  ggrepel::geom_text_repel(
    data = sel_all,
    aes(UMAP_1, UMAP_2, label = TCR_label2, color = label_color),
    size = 3.2,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.20,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.35,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  labs(
    title = paste0("UMAP (FILTERED + RE-CLUSTERED): Patient=", patient_keep, " | all selected TCRs; VAL in red"),
    x = "UMAP_1", y = "UMAP_2"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

save_plot_best(
  paste0(out_prefix, "_umap_selectedTCR_allLabels_NEWUMAP_NEWclusters"),
  p_umap_selected, width = 8.5, height = 6.5, dpi = png_dpi
)

## =========================================================
## 9) UMAP #2: functional.cluster + VAL labels red
## =========================================================
um$state <- as.character(seu@meta.data[um$cell, state_col])
um$state <- ifelse(is.na(um$state) | trimws(um$state) == "", "NA", um$state)

state_levels <- c(sort(unique(um$state[um$state != "NA"])), "NA")
um$state <- factor(um$state, levels = state_levels)

non_na_levels <- state_levels[state_levels != "NA"]
pal_state <- if (length(non_na_levels) > 0) {
  c(
    setNames(scales::hue_pal()(length(non_na_levels)), non_na_levels),
    setNames("grey80", "NA")
  )
} else {
  c(setNames("grey80", "NA"))
}

sel_val_um2 <- big_table %>%
  filter(VAL %in% TRUE) %>%
  select(TCR_unique, TCR_excel, rep_cell) %>%
  left_join(um, by = c("rep_cell" = "cell")) %>%
  filter(!is.na(UMAP_1), !is.na(UMAP_2)) %>%
  mutate(
    tcr_num    = str_extract(as.character(TCR_excel), "\\d+"),
    TCR_label2 = ifelse(!is.na(tcr_num), paste0("TCR", tcr_num), paste0("TCR", as.character(TCR_excel)))
  )

p_umap_state_val <- ggplot(um, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point(size = 0.28) +
  scale_color_manual(values = pal_state, drop = FALSE) +
  { if (nrow(sel_val_um2) > 0)
    geom_point(
      data = sel_val_um2, aes(UMAP_1, UMAP_2),
      inherit.aes = FALSE, color = deep_red, size = 2.6, show.legend = FALSE
    )
  } +
  { if (nrow(sel_val_um2) > 0)
    ggrepel::geom_text_repel(
      data = sel_val_um2,
      aes(UMAP_1, UMAP_2, label = TCR_label2),
      inherit.aes = FALSE,
      color = deep_red,
      size = 3.4,
      fontface = "bold",
      box.padding = 0.40,
      point.padding = 0.25,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = deep_red,
      segment.size = 0.45,
      show.legend = FALSE
    )
  } +
  labs(
    title = paste0("UMAP (FILTERED + RE-CLUSTERED): Patient=", patient_keep, " | functional.cluster + VAL labels (red)"),
    x = "UMAP_1", y = "UMAP_2", color = "Functional state"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4.8))) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_plot_best(
  paste0(out_prefix, "_umap_functionalStateClusters_VALlabelsRed_NEWUMAP"),
  p_umap_state_val, width = 8.5, height = 6.5, dpi = png_dpi
)

## =========================================================
## 10) Export tables
## =========================================================
wb_out1 <- paste0(out_prefix, "_outputs_markerModules.xlsx")
openxlsx::write.xlsx(
  list(
    big_table = big_table,
    heatmap_matrix = as.data.frame(expr_mat) %>% tibble::rownames_to_column("gene"),
    missing_matches = missing_df
  ),
  file = wb_out1
)

## =========================================================
## 11) UMAP background colored by new clusters
##     + overlay tumor_sc TRUE (black) and VAL (red)
## =========================================================
col_nonval <- "black"
label_tcr_on_umap <- TRUE
label_only_val <- FALSE
label_font_nonval <- 2.4
label_font_val <- 3.0
label_max_overlaps <- Inf

key_col <- match_key_for_tumor
tab_keys <- tumor_tab %>%
  mutate(key = as.character(.data[[key_col]])) %>%
  filter(!is.na(key), trimws(key) != "") %>%
  select(key, TCR_label, VAL) %>%
  distinct()

md <- seu@meta.data
meta_key2 <- as.character(md[[key_col]])
meta_key2[is.na(meta_key2)] <- ""

hits_all <- bind_rows(lapply(seq_len(nrow(tab_keys)), function(i) {
  k <- as.character(tab_keys$key[i])
  idx <- if (match_mode_for_tumor == "exact") {
    which(meta_key2 == k)
  } else {
    which(grepl(k, meta_key2, fixed = TRUE))
  }
  if (length(idx) == 0) return(NULL)
  tibble(
    cell = rownames(md)[idx],
    TCR_label = tab_keys$TCR_label[i],
    VAL = as.logical(tab_keys$VAL[i]),
    key = k
  )
}))

if (is.null(hits_all) || nrow(hits_all) == 0) {
  stop("No tumor_sc hits found in this filtered Seurat object using key=", key_col)
}
hits_one <- hits_all %>% group_by(cell) %>% slice(1) %>% ungroup()

um2 <- um %>% select(cell, UMAP_1, UMAP_2)
um2$cluster_raw <- as.character(seu@meta.data[um2$cell, "seurat_clusters"])
um2$cluster_raw[is.na(um2$cluster_raw) | trimws(um2$cluster_raw) == ""] <- "NA"

clu2 <- normalize_cluster_id(um2$cluster_raw, keep_na_last = TRUE)
um2$cluster <- factor(clu2$id, levels = clu2$levels)

pal_cluster <- setNames(scales::hue_pal()(length(clu2$levels)), clu2$levels)

um2$tumor_hit <- FALSE
um2$val_hit <- FALSE
um2$TCR_label <- NA_character_

idxm <- match(hits_one$cell, um2$cell)
ok <- which(!is.na(idxm))
um2$tumor_hit[idxm[ok]] <- TRUE
um2$val_hit[idxm[ok]] <- hits_one$VAL[ok] %in% TRUE
um2$TCR_label[idxm[ok]] <- hits_one$TCR_label[ok]

df_nonval <- um2 %>% filter(tumor_hit, !val_hit)
df_val    <- um2 %>% filter(tumor_hit,  val_hit)

df_lab2 <- um2 %>%
  filter(tumor_hit) %>%
  mutate(
    label_group = ifelse(val_hit, "VAL", "nonVAL"),
    label = TCR_label
  )

if (isTRUE(label_only_val)) {
  df_lab2 <- df_lab2 %>% filter(val_hit)
}

p_umap_cluster_bg_tumor_val <- ggplot(um2, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cluster), size = 0.25, alpha = 0.95) +
  scale_color_manual(values = pal_cluster, drop = FALSE, labels = function(z) paste0("c", z)) +
  guides(color = guide_legend(override.aes = list(size = 3.6))) +
  { if (nrow(df_nonval) > 0)
    geom_point(
      data = df_nonval, inherit.aes = FALSE,
      aes(UMAP_1, UMAP_2), color = col_nonval, size = 1.25, alpha = 0.95
    )
  } +
  { if (nrow(df_val) > 0)
    geom_point(
      data = df_val, inherit.aes = FALSE,
      aes(UMAP_1, UMAP_2), color = deep_red, size = 2.35, alpha = 0.98
    )
  } +
  { if (isTRUE(label_tcr_on_umap) && nrow(df_lab2) > 0)
    ggrepel::geom_text_repel(
      data = df_lab2,
      aes(UMAP_1, UMAP_2, label = label),
      inherit.aes = FALSE,
      size = ifelse(df_lab2$label_group == "VAL", label_font_val, label_font_nonval),
      color = ifelse(df_lab2$label_group == "VAL", deep_red, col_nonval),
      fontface = "bold",
      box.padding = 0.22,
      point.padding = 0.18,
      max.overlaps = label_max_overlaps,
      min.segment.length = 0,
      segment.alpha = 0.7,
      segment.size = 0.25
    )
  } +
  labs(
    title = paste0("NEW UMAP (FILTERED + RE-CLUSTERED): Patient=", patient_keep,
                   " | background=new clusters; tumor_sc black; VAL red"),
    subtitle = paste0("match_key=", key_col,
                      " | tumor_sc matched cells=", sum(um2$tumor_hit),
                      " | VAL cells=", sum(um2$val_hit),
                      " | resolution=", cluster_resolution),
    x = "UMAP_1", y = "UMAP_2", color = "New cluster"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

save_plot_best(
  paste0(out_prefix, "_umap_NEWclusters_BGcolor_tumorBlack_VALred"),
  p_umap_cluster_bg_tumor_val, width = 9.5, height = 6.8, dpi = png_dpi
)

hits_one2 <- hits_one %>%
  mutate(cluster_raw = as.character(seu@meta.data[cell, "seurat_clusters"])) %>%
  mutate(cluster_raw = ifelse(is.na(cluster_raw) | trimws(cluster_raw) == "", "NA", cluster_raw))

clu_bar <- normalize_cluster_id(hits_one2$cluster_raw, keep_na_last = TRUE)
hits_one2 <- hits_one2 %>%
  mutate(
    cluster = factor(clu_bar$id, levels = clu_bar$levels),
    class = ifelse(VAL %in% TRUE, "VAL", "tumor_sc_nonVAL")
  ) %>%
  count(cluster, class, name = "n") %>%
  tidyr::complete(
    cluster = factor(clu_bar$levels, levels = clu_bar$levels),
    class = c("tumor_sc_nonVAL", "VAL"),
    fill = list(n = 0)
  ) %>%
  mutate(class = factor(class, levels = c("tumor_sc_nonVAL", "VAL")))

tot_df <- hits_one2 %>%
  group_by(cluster) %>%
  summarise(n_tot = sum(n), .groups = "drop")

pal_rb <- c(tumor_sc_nonVAL = col_nonval, VAL = deep_red)

p_bar_cluster_val_nonval <- ggplot(hits_one2, aes(x = cluster, y = n, fill = class)) +
  geom_col(width = 0.85) +
  geom_text(
    data = tot_df,
    aes(x = cluster, y = n_tot, label = n_tot),
    inherit.aes = FALSE,
    vjust = -0.30,
    fontface = "bold",
    size = 3.2,
    color = "black"
  ) +
  scale_x_discrete(labels = function(z) paste0("c", z)) +
  scale_fill_manual(
    values = pal_rb,
    breaks = c("VAL", "tumor_sc_nonVAL"),
    labels = c("VAL (red)", "tumor_sc non-VAL (black)")
  ) +
  labs(
    title = paste0("Per NEW cluster (FILTERED + RE-CLUSTERED): Patient=", patient_keep,
                   " | tumor_sc hits = non-VAL black + VAL red"),
    subtitle = paste0("resolution=", cluster_resolution, " | match_key=", key_col),
    x = "New cluster",
    y = "# tumor_sc matched cells",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold", size = 9),
    legend.position = "top",
    legend.text = element_text(face = "bold")
  ) +
  coord_cartesian(clip = "off")

bar_w <- max(8.5, 0.30 * length(clu_bar$levels) + 4)
save_plot_best(
  paste0(out_prefix, "_BAR_NEWcluster_tumorBlack_VALred_stacked"),
  p_bar_cluster_val_nonval, width = bar_w, height = 5.8, dpi = png_dpi
)

out_xlsx2 <- paste0(out_prefix, "_NEWcluster_tumorSC_VAL_counts.xlsx")
openxlsx::write.xlsx(
  list(
    tumor_tab_used = tumor_tab,
    hits_onePerCell = hits_one %>% mutate(cluster_raw = as.character(seu@meta.data[cell, "seurat_clusters"])),
    bar_table = hits_one2,
    totals = tot_df
  ),
  file = out_xlsx2
)

## =========================================================
## 12) Marker panels per category
## =========================================================
out_merged <- paste0(out_prefix, "_marker_panels_perCategory_RNA")
dir.create(out_merged, showWarnings = FALSE, recursive = TRUE)

LEG_H_MM <- 5.0
LEG_W_MM <- 2.0
LEG_TXT  <- 7

spec_merged <- list(
  "Lineage and Naive/Memory markers" = list(items = list(
    list(key="CD4",   title="CD4 (RNA)"),
    list(key="CD8a",  title="CD8A (RNA)"),
    list(key="IL7Aa", title="IL7R (RNA)"),
    list(key="CCR7",  title="CCR7 (RNA)"),
    list(key="TCF7",  title="TCF7 (RNA)"),
    list(key="CD28",  title="CD28 (RNA)"),
    list(key="CD27",  title="CD27 (RNA)"),
    list(key="CD62L", title="SELL (RNA)"),
    list(key="CD95",  title="FAS (RNA)")
  )),
  "Exhaustion/activation markers" = list(items = list(
    list(key="PD1",   title="PDCD1 (RNA)"),
    list(key="TIGIT", title="TIGIT (RNA)"),
    list(key="CD39",  title="ENTPD1 (RNA)"),
    list(key="ICOS",  title="ICOS (RNA)"),
    list(key="CD25",  title="IL2RA (RNA)"),
    list(key="CD137", title="TNFRSF9 (RNA)"),
    list(key="TOX",   title="TOX (RNA)")
  )),
  "Effector molecules" = list(items = list(
    list(key="PRF1",  title="PRF1 (RNA)"),
    list(key="IFNG",  title="IFNG (RNA)"),
    list(key="GZMB",  title="GZMB (RNA)"),
    list(key="GZMK",  title="GZMK (RNA)"),
    list(key="FASLG", title="FASLG (RNA)"),
    list(key="GZMH",  title="GZMH (RNA)"),
    list(key="GZMM",  title="GZMM (RNA)")
  )),
  "Transcription factors / Proliferation" = list(items = list(
    list(key="MKI67", title="MKI67 (RNA)"),
    list(key="TOX",   title="TOX (RNA)"),
    list(key="PRDM1", title="PRDM1 (RNA)"),
    list(key="NR4A1", title="NR4A1 (RNA)"),
    list(key="FOXP3", title="FOXP3 (RNA)"),
    list(key="TCF7",  title="TCF7 (RNA)")
  ))
)

dedup_items <- function(items) {
  keys <- vapply(items, `[[`, "", "key")
  items[!duplicated(keys)]
}
spec_merged <- lapply(spec_merged, function(x) {
  x$items <- dedup_items(x$items)
  x
})

plot_one_category <- function(seu, category_name, items,
                              max_cols = 5,
                              pt_size = 0.18,
                              min_q = 0.01, max_q = 0.99) {
  
  stopifnot(umap_name %in% Reductions(seu))
  DefaultAssay(seu) <- "RNA"
  
  emb <- as.data.frame(Embeddings(seu, umap_name))[, 1:2, drop = FALSE]
  colnames(emb) <- c("UMAP_1","UMAP_2")
  emb$cell <- rownames(emb)
  
  titles <- vapply(items, `[[`, "", "title")
  vals_list <- vector("list", length(items))
  
  for (i in seq_along(items)) {
    gene_cands <- resolve_to_rna_gene(items[[i]]$key)
    if (is.null(gene_cands)) {
      vals_list[[i]] <- rep(NA_real_, nrow(emb))
      next
    }
    feat <- find_feature(seu, "RNA", gene_cands)
    if (is.na(feat)) {
      vals_list[[i]] <- rep(NA_real_, nrow(emb))
      next
    }
    
    vdf <- tryCatch(FetchData(seu, vars = feat, layer = "data"), error = function(e) NULL)
    if (is.null(vdf)) vdf <- tryCatch(FetchData(seu, vars = feat), error = function(e) NULL)
    if (is.null(vdf) || nrow(vdf) == 0) {
      vals_list[[i]] <- rep(NA_real_, nrow(emb))
      next
    }
    
    v <- as.numeric(vdf[, 1])
    names(v) <- rownames(vdf)
    vals_list[[i]] <- v[match(emb$cell, names(v))]
  }
  
  allv <- unlist(vals_list, use.names = FALSE)
  allv <- allv[is.finite(allv)]
  if (length(allv) < 10) {
    lo <- 0
    hi <- 1
  } else {
    lo <- as.numeric(quantile(allv, probs = min_q, na.rm = TRUE))
    hi <- as.numeric(quantile(allv, probs = max_q, na.rm = TRUE))
    if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
      lo <- min(allv)
      hi <- max(allv)
      if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
        lo <- 0
        hi <- 1
      }
    }
  }
  
  plist <- vector("list", length(items))
  for (i in seq_along(items)) {
    df <- emb
    df$val <- vals_list[[i]]
    
    if (!any(is.finite(df$val))) {
      plist[[i]] <- ggplot() + theme_void() +
        labs(title = paste0(titles[i], "\n[MISSING]")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
      next
    }
    
    plist[[i]] <- ggplot(df, aes(UMAP_1, UMAP_2, color = val)) +
      geom_point(size = pt_size) +
      coord_equal() +
      scale_color_gradient(
        low = "grey90", high = "#2C5BFF",
        limits = c(lo, hi), oob = scales::squish
      ) +
      labs(title = titles[i], color = NULL) +
      theme_void(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.key.height = unit(LEG_H_MM, "mm"),
        legend.key.width  = unit(LEG_W_MM, "mm"),
        legend.text = element_text(size = LEG_TXT),
        plot.margin = margin(2, 2, 2, 2)
      )
  }
  
  combined <- wrap_plots(plist, ncol = min(max_cols, length(plist))) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = category_name,
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
    )
  
  combined & theme(
    legend.position = "right",
    legend.box.just = "center",
    plot.margin = margin(3, 3, 3, 3)
  )
}

make_dims <- function(n_panels, max_cols = 5) {
  ncol <- min(max_cols, n_panels)
  nrow <- ceiling(n_panels / ncol)
  width  <- max(6.8, ncol * 2.55 + 1.4)
  height <- max(4.2, nrow * 2.15 + 0.9)
  list(w = width, h = height)
}

for (catn in names(spec_merged)) {
  items <- spec_merged[[catn]]$items
  p <- plot_one_category(seu, catn, items, max_cols = 5)
  dims <- make_dims(length(items), max_cols = 5)
  fname <- paste0("MERGED_CD4CD8_", gsub("[^A-Za-z0-9]+", "_", catn), "_RNA")
  save_plot_best(file.path(out_merged, fname), p, dims$w, dims$h, dpi = png_dpi)
}

## =========================================================
## 13) per-cell ssGSEA heatmap (ComplexHeatmap)
## =========================================================
if (requireNamespace("GSVA", quietly = TRUE) &&
    requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE) &&
    file.exists(sig_merged_xlsx)) {
  
  suppressPackageStartupMessages({
    library(GSVA)
    library(ComplexHeatmap)
    library(circlize)
  })
  
  safe_str <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x[is.na(x)] <- ""
    x
  }
  
  clean_gene <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x[x %in% c("", "NA", "N/A", "NULL", "NONE", "0")] <- NA_character_
    x
  }
  
  read_signatures_with_order <- function(sig_xlsx) {
    sig_wide <- suppressMessages(readxl::read_xlsx(sig_xlsx, col_names = FALSE))
    if (nrow(sig_wide) < 3) stop("Signature file must have >=3 rows (paper/state/genes...)")
    
    paper_row <- as.character(unlist(sig_wide[1, ])); paper_row[is.na(paper_row)] <- ""
    state_row <- as.character(unlist(sig_wide[2, ])); state_row[is.na(state_row)] <- ""
    
    keep_cols <- which(trimws(paper_row) != "" & trimws(state_row) != "")
    if (length(keep_cols) == 0) stop("No valid signature columns (paper+state missing).")
    
    genes_block <- sig_wide[3:nrow(sig_wide), keep_cols, drop = FALSE]
    gene_list_raw <- lapply(seq_along(keep_cols), function(j) {
      g <- clean_gene(unlist(genes_block[, j]))
      g <- g[!is.na(g)]
      unique(g)
    })
    
    tibble::tibble(
      col_i   = keep_cols,
      paper   = paper_row[keep_cols],
      state   = state_row[keep_cols],
      n_genes = vapply(gene_list_raw, length, 1L),
      genes   = gene_list_raw
    ) %>%
      filter(n_genes > 0) %>%
      mutate(sig_name = make.unique(paste0(paper, " | ", state), sep = " #")) %>%
      arrange(col_i)
  }
  
  drop_constant_genes_dense <- function(expr_mat) {
    if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) return(expr_mat)
    v <- apply(expr_mat, 1, function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) return(0)
      stats::var(x)
    })
    keep <- which(is.finite(v) & v > 0)
    expr_mat[keep, , drop = FALSE]
  }
  
  key_col_excel <- match_key_for_tumor
  key_col_meta  <- match_key_for_tumor
  
  tumor_keys <- tumor_tab %>%
    transmute(key = safe_str(.data[[key_col_excel]])) %>%
    filter(key != "") %>%
    distinct()
  
  md3 <- seu@meta.data
  meta_key3 <- safe_str(md3[[key_col_meta]])
  
  tumor_cells <- unique(unlist(lapply(seq_len(nrow(tumor_keys)), function(i) {
    k <- tumor_keys$key[i]
    idx <- if (match_mode_for_tumor == "exact") {
      which(meta_key3 == k)
    } else {
      which(grepl(k, meta_key3, fixed = TRUE))
    }
    rownames(md3)[idx]
  })))
  tumor_cells <- tumor_cells[!is.na(tumor_cells)]
  if (length(tumor_cells) == 0) stop("ssGSEA: no matched tumor_sc cells in filtered object.")
  
  expr <- GetLayerDataSafe(seu, assay = "RNA", layer = expr_layer)
  expr <- expr[, intersect(tumor_cells, colnames(expr)), drop = FALSE]
  rownames(expr) <- toupper(rownames(expr))
  
  sig_tbl <- read_signatures_with_order(sig_merged_xlsx)
  feat <- rownames(expr)
  
  sig_tbl2 <- sig_tbl %>%
    mutate(
      genes_present = lapply(genes, function(g) intersect(g, feat)),
      n_present = vapply(genes_present, length, 1L)
    ) %>%
    filter(n_present >= min_present_genes)
  
  if (nrow(sig_tbl2) == 0) stop("ssGSEA: no signatures pass min_present_genes=", min_present_genes)
  
  gene_sets <- sig_tbl2$genes_present
  names(gene_sets) <- sig_tbl2$sig_name
  
  expr_mat_ss <- as.matrix(expr)
  if (drop_const_genes) {
    expr2 <- drop_constant_genes_dense(expr_mat_ss)
    if (nrow(expr2) >= 200) expr_mat_ss <- expr2
  }
  
  param <- tryCatch(
    GSVA::ssgseaParam(expr = expr_mat_ss, geneSets = gene_sets, normalize = ssgsea_norm),
    error = function(e1) {
      GSVA::ssgseaParam(expr = expr_mat_ss, gset.idx.list = gene_sets, normalize = ssgsea_norm)
    }
  )
  ssg <- GSVA::gsva(param, verbose = FALSE)
  
  mat <- as.matrix(ssg)
  mat <- mat[sig_tbl2$sig_name, , drop = FALSE]
  
  mdx <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    filter(cell %in% colnames(mat))
  
  umx <- as.data.frame(Embeddings(seu, umap_name))[, 1:2, drop = FALSE]
  colnames(umx) <- c("UMAP_1","UMAP_2")
  umx$cell <- rownames(umx)
  
  mdx <- mdx %>%
    left_join(umx, by = "cell") %>%
    mutate(
      UMAP_1 = ifelse(is.na(UMAP_1), 0, UMAP_1),
      UMAP_2 = ifelse(is.na(UMAP_2), 0, UMAP_2),
      cluster_raw = as.character(seurat_clusters)
    )
  
  clu_mdx <- normalize_cluster_id(mdx$cluster_raw, keep_na_last = TRUE)
  mdx$cluster <- factor(clu_mdx$id, levels = clu_mdx$levels)
  
  mdx <- mdx %>% arrange(cluster, UMAP_1, UMAP_2, cell)
  cell_order <- mdx$cell
  mat <- mat[, cell_order, drop = FALSE]
  
  col_split <- mdx$cluster
  names(col_split) <- mdx$cell
  
  cell2tcr <- rep("", length(cell_order))
  names(cell2tcr) <- cell_order
  if (nrow(hits_one) > 0) {
    m2 <- hits_one %>% select(cell, TCR_label)
    idx2 <- match(m2$cell, names(cell2tcr))
    ok2 <- which(!is.na(idx2))
    cell2tcr[idx2[ok2]] <- m2$TCR_label[ok2]
  }
  
  cluster_text <- paste0("c", as.character(col_split))
  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = ComplexHeatmap::anno_text(cluster_text, rot = 90, just = "right",
                                        gp = gpar(fontsize = 7, fontface = "bold")),
    TCR = ComplexHeatmap::anno_text(cell2tcr[cell_order], rot = 90, just = "right",
                                    gp = gpar(fontsize = 7)),
    show_annotation_name = FALSE,
    show_legend = c(Cluster = FALSE, TCR = FALSE),
    annotation_height = unit.c(unit(10, "mm"), unit(10, "mm"))
  )
  
  paper_levels <- unique(sig_tbl2$paper)
  row_split <- factor(sig_tbl2$paper, levels = paper_levels)
  
  row_state <- sig_tbl2$state
  names(row_state) <- sig_tbl2$sig_name
  
  ra_left <- ComplexHeatmap::rowAnnotation(
    Type = ComplexHeatmap::anno_text(row_state, just = "left", gp = gpar(fontsize = 9)),
    width = unit(55, "mm"),
    show_annotation_name = FALSE
  )
  
  ra_right <- ComplexHeatmap::rowAnnotation(
    Paper = ComplexHeatmap::anno_block(
      gp = gpar(fill = NA, col = "black", lwd = 1.0),
      labels = levels(row_split),
      labels_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_rot = 0,
      labels_just = "center"
    ),
    width = unit(34, "mm"),
    show_annotation_name = FALSE
  )
  
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  
  zlim <- 2
  col_fun_z <- circlize::colorRamp2(c(-zlim, 0, zlim), c("#2C7BB6", "white", "#D7191C"))
  
  ht <- ComplexHeatmap::Heatmap(
    mat_z,
    name = "z",
    col = col_fun_z,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    row_split = row_split,
    row_gap = unit(2.2, "mm"),
    column_split = col_split,
    column_gap = unit(0.9, "mm"),
    top_annotation = top_ha,
    left_annotation = ra_left,
    right_annotation = ra_right,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      title = "z",
      legend_direction = "horizontal",
      title_position = "topcenter",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      grid_height = unit(6, "mm"),
      grid_width  = unit(2.2, "mm")
    ),
    use_raster = TRUE,
    raster_quality = 2
  )
  
  nC <- ncol(mat_z)
  nR <- nrow(mat_z)
  w_in <- min(40, max(10, 6 + nC * 0.012))
  h_in <- min(26, max(7,  4 + nR * 0.22))
  
  out_ssg <- paste0(out_prefix, "_TUMORSC_allCells_ssGSEA_perCell_heatmap_Z_NEWclusters_NOclusterBAR")
  save_ht_best(out_ssg, ht, width = w_in, height = h_in, dpi = 600)
  
} else {
  message("[skip] ssGSEA heatmap: missing GSVA/ComplexHeatmap/circlize or signature file path.")
}

## =========================================================
## 14) New clusters x signatures heatmap (NO ssGSEA)
## =========================================================
if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE) &&
    file.exists(sig_merged_xlsx)) {
  
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
  })
  
  clean_gene2 <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x[x %in% c("", "NA", "N/A", "NULL", "NONE", "0")] <- NA_character_
    x
  }
  
  read_signatures_with_order2 <- function(sig_xlsx) {
    sig_wide <- suppressMessages(readxl::read_xlsx(sig_xlsx, col_names = FALSE))
    if (nrow(sig_wide) < 3) stop("Signature file must have >=3 rows (paper/state/genes...)")
    
    paper_row <- as.character(unlist(sig_wide[1, ])); paper_row[is.na(paper_row)] <- ""
    state_row <- as.character(unlist(sig_wide[2, ])); state_row[is.na(state_row)] <- ""
    
    keep_cols <- which(trimws(paper_row) != "" & trimws(state_row) != "")
    if (length(keep_cols) == 0) stop("No valid signature columns (paper+state missing).")
    
    genes_block <- sig_wide[3:nrow(sig_wide), keep_cols, drop = FALSE]
    gene_list_raw <- lapply(seq_along(keep_cols), function(j) {
      g <- clean_gene2(unlist(genes_block[, j]))
      g <- g[!is.na(g)]
      unique(g)
    })
    
    tibble::tibble(
      col_i   = keep_cols,
      paper   = paper_row[keep_cols],
      state   = state_row[keep_cols],
      n_genes = vapply(gene_list_raw, length, 1L),
      genes   = gene_list_raw
    ) %>%
      filter(n_genes > 0) %>%
      mutate(sig_name = make.unique(paste0(paper, " | ", state), sep = " #")) %>%
      arrange(col_i)
  }
  
  seu$cluster_for_heatmap_raw <- as.character(seu$seurat_clusters)
  clu_hm <- normalize_cluster_id(seu$cluster_for_heatmap_raw, keep_na_last = TRUE)
  seu$cluster_for_heatmap <- factor(clu_hm$id, levels = clu_hm$levels)
  seu$cluster_for_heatmap <- droplevels(seu$cluster_for_heatmap)
  
  avg_list <- Seurat::AverageExpression(
    seu,
    assays = "RNA",
    layer = expr_layer,
    group.by = "cluster_for_heatmap",
    verbose = FALSE,
    return.seurat = FALSE
  )
  avg_expr <- avg_list$RNA
  if (is.null(avg_expr) || nrow(avg_expr) == 0) {
    stop("AverageExpression returned empty matrix.")
  }
  rownames(avg_expr) <- toupper(rownames(avg_expr))
  
  colnames(avg_expr) <- as.character(colnames(avg_expr))
  col_num <- suppressWarnings(as.integer(stringr::str_extract(colnames(avg_expr), "\\d+")))
  ord <- order(is.na(col_num), col_num, colnames(avg_expr))
  avg_expr <- avg_expr[, ord, drop = FALSE]
  
  sig_tbl <- read_signatures_with_order2(sig_merged_xlsx)
  feat <- rownames(avg_expr)
  
  sig_tbl2 <- sig_tbl %>%
    mutate(
      genes_present = lapply(genes, function(g) intersect(g, feat)),
      n_present = vapply(genes_present, length, 1L)
    ) %>%
    filter(n_present >= min_present_genes)
  
  if (nrow(sig_tbl2) == 0) stop("No signatures pass min_present_genes=", min_present_genes)
  
  mat <- vapply(seq_len(nrow(sig_tbl2)), function(i) {
    g <- sig_tbl2$genes_present[[i]]
    colMeans(avg_expr[g, , drop = FALSE])
  }, FUN.VALUE = numeric(ncol(avg_expr)))
  
  mat <- t(mat)
  rownames(mat) <- sig_tbl2$sig_name
  colnames(mat) <- colnames(avg_expr)
  
  paper_levels <- unique(sig_tbl2$paper)
  row_split <- factor(sig_tbl2$paper, levels = paper_levels)
  
  row_state <- sig_tbl2$state
  names(row_state) <- sig_tbl2$sig_name
  
  col_split <- factor(colnames(mat), levels = colnames(mat))
  
  cluster_text <- paste0("c", as.character(col_split))
  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    Cluster = ComplexHeatmap::anno_text(cluster_text, rot = 90, just = "right",
                                        gp = gpar(fontsize = 7, fontface = "bold")),
    show_annotation_name = FALSE,
    show_legend = c(Cluster = FALSE),
    annotation_height = unit(10, "mm")
  )
  
  ra_right <- ComplexHeatmap::rowAnnotation(
    Paper = ComplexHeatmap::anno_block(
      gp = gpar(fill = NA, col = "black", lwd = 1.0),
      labels = levels(row_split),
      labels_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_rot = 0,
      labels_just = "center"
    ),
    width = unit(34, "mm"),
    show_annotation_name = FALSE
  )
  
  ra_left <- ComplexHeatmap::rowAnnotation(
    Type = ComplexHeatmap::anno_text(row_state, just = "left", gp = gpar(fontsize = 9)),
    width = unit(55, "mm"),
    show_annotation_name = FALSE
  )
  
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  
  zlim <- 2
  col_fun_z <- circlize::colorRamp2(c(-zlim, 0, zlim), c("#2C7BB6", "white", "#D7191C"))
  
  ht2 <- ComplexHeatmap::Heatmap(
    mat_z,
    name = "z",
    col = col_fun_z,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    row_split = row_split,
    row_gap = unit(2.2, "mm"),
    column_split = col_split,
    column_gap = unit(0.9, "mm"),
    top_annotation = top_ha,
    left_annotation = ra_left,
    right_annotation = ra_right,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      title = "z",
      legend_direction = "horizontal",
      title_position = "topcenter",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      grid_height = unit(6, "mm"),
      grid_width  = unit(2.2, "mm")
    ),
    use_raster = TRUE,
    raster_quality = 2
  )
  
  nC <- ncol(mat_z)
  nR <- nrow(mat_z)
  w_in <- min(20, max(8, 5 + nC * 0.35))
  h_in <- min(26, max(7, 4 + nR * 0.22))
  
  out_sig <- paste0(out_prefix, "_NEWclusters_signatureList_NOssGSEA_heatmap_Z_NOclusterBAR")
  save_ht_best(out_sig, ht2, width = w_in, height = h_in, dpi = 600)
  
} else {
  message("[skip] NEW clusters x signatures heatmap: missing ComplexHeatmap/circlize or signature file path.")
}

## =========================================================
## 15) NeoTCR scoring + UMAP projection + clonotype ranking
## =========================================================
if (!all(c("NeoTCR4_score", "NeoTCR8_score") %in% colnames(seu@meta.data))) {
  
  score_sig_xlsx <- "C:/Users/peppa/Downloads/scoresignature_xu_merged.xlsx"
  stopifnot(file.exists(score_sig_xlsx))
  
  sheets <- readxl::excel_sheets(score_sig_xlsx)
  stopifnot(all(c("NeoTCR4","NeoTCR8") %in% sheets))
  
  neo4_genes_u <- clean_gene_vec(readxl::read_xlsx(score_sig_xlsx, sheet = "NeoTCR4", col_names = TRUE)[[1]])
  neo8_genes_u <- clean_gene_vec(readxl::read_xlsx(score_sig_xlsx, sheet = "NeoTCR8", col_names = TRUE)[[1]])
  
  neo4_feats <- map_genes_to_features(seu, neo4_genes_u, "RNA")
  neo8_feats <- map_genes_to_features(seu, neo8_genes_u, "RNA")
  
  if (length(neo4_feats) < 8)  stop("Too few NeoTCR4 genes present: ", length(neo4_feats))
  if (length(neo8_feats) < 15) stop("Too few NeoTCR8 genes present: ", length(neo8_feats))
  
  DefaultAssay(seu) <- "RNA"
  
  seu <- AddModuleScore(seu, features = list(NeoTCR4 = neo4_feats), name = "NeoTCR4", assay = "RNA", search = FALSE)
  seu <- AddModuleScore(seu, features = list(NeoTCR8 = neo8_feats), name = "NeoTCR8", assay = "RNA", search = FALSE)
  
  col4 <- grep("^NeoTCR4(_)?1$", colnames(seu@meta.data), value = TRUE)
  col8 <- grep("^NeoTCR8(_)?1$", colnames(seu@meta.data), value = TRUE)
  if (length(col4) != 1 || length(col8) != 1) {
    stop(
      "Cannot uniquely identify NeoTCR4/8 AddModuleScore columns. Found: ",
      paste(grep("NeoTCR4|NeoTCR8", colnames(seu@meta.data), value = TRUE), collapse = ", ")
    )
  }
  seu$NeoTCR4_score <- as.numeric(seu@meta.data[[col4]])
  seu$NeoTCR8_score <- as.numeric(seu@meta.data[[col8]])
}

meta_cols <- colnames(seu@meta.data)

key_priority <- c("vjaa","CDR3b","CDR3B","Beta_clonotype","guide_id","TCR")
key_candidates <- key_priority[key_priority %in% colnames(tumor_tab)]
key_in_meta <- key_candidates[key_candidates %in% meta_cols]

if (length(key_in_meta) == 0) {
  overlap <- intersect(colnames(tumor_tab), meta_cols)
  overlap <- setdiff(overlap, c("VAL","tumor_sc","skin_sc","in_vitro_exp_sc_CR"))
  if (length(overlap) == 0) {
    stop(
      "No shared key column between tumor_tab (Excel) and seu@meta.data.\n",
      "You must add one TCR identifier column into seu@meta.data (e.g., vjaa or CDR3b) OR rerun the matching step on the pre-filter object."
    )
  }
  key_use <- overlap[1]
} else {
  key_use <- key_in_meta[1]
}
message("[key] using shared clonotype key = ", key_use)

tab_keys2 <- tumor_tab %>%
  mutate(key = as.character(.data[[key_use]])) %>%
  filter(!is.na(key), trimws(key) != "") %>%
  transmute(key, VAL = as.logical(VAL)) %>%
  group_by(key) %>%
  summarise(
    tumor_hit = TRUE,
    val_hit = any(VAL %in% TRUE),
    .groups = "drop"
  )

md_rank <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(key = as.character(.data[[key_use]]))

md_rank <- md_rank %>%
  left_join(tab_keys2, by = "key") %>%
  mutate(
    tumor_hit = ifelse(is.na(tumor_hit), FALSE, tumor_hit),
    val_hit   = ifelse(is.na(val_hit),   FALSE, val_hit)
  )

seu$tumor_hit <- md_rank$tumor_hit[match(rownames(seu@meta.data), md_rank$cell)]
seu$val_hit   <- md_rank$val_hit[match(rownames(seu@meta.data), md_rank$cell)]

um_rank <- as.data.frame(Embeddings(seu, reduction = umap_name))[,1:2, drop=FALSE]
colnames(um_rank) <- c("UMAP_1","UMAP_2")
um_rank$cell <- rownames(um_rank)

um_rank$NeoTCR4_score <- seu@meta.data[um_rank$cell, "NeoTCR4_score"]
um_rank$NeoTCR8_score <- seu@meta.data[um_rank$cell, "NeoTCR8_score"]
um_rank$tumor_hit <- seu@meta.data[um_rank$cell, "tumor_hit"]
um_rank$val_hit   <- seu@meta.data[um_rank$cell, "val_hit"]

df_nonval2 <- um_rank %>% filter(tumor_hit, !val_hit)
df_val2    <- um_rank %>% filter(tumor_hit,  val_hit)

pt_bg <- 2
pt_hit_nonval <- 1.5
pt_hit_val <- 2.5

p_um4 <- ggplot(um_rank, aes(UMAP_1, UMAP_2, color = NeoTCR4_score)) +
  geom_point(size = pt_bg) +
  { if (nrow(df_nonval2)>0) geom_point(
    data = df_nonval2, aes(UMAP_1,UMAP_2),
    inherit.aes = FALSE, color = "black", size = pt_hit_nonval
  ) } +
  { if (nrow(df_val2)>0) geom_point(
    data = df_val2, aes(UMAP_1,UMAP_2),
    inherit.aes = FALSE, color = deep_red, size = pt_hit_val
  ) } +
  labs(title = "NeoTCR4 score", color = "NeoTCR4") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_um8 <- ggplot(um_rank, aes(UMAP_1, UMAP_2, color = NeoTCR8_score)) +
  geom_point(size = pt_bg) +
  { if (nrow(df_nonval2)>0) geom_point(
    data = df_nonval2, aes(UMAP_1,UMAP_2),
    inherit.aes = FALSE, color = "black", size = pt_hit_nonval
  ) } +
  { if (nrow(df_val2)>0) geom_point(
    data = df_val2, aes(UMAP_1,UMAP_2),
    inherit.aes = FALSE, color = deep_red, size = pt_hit_val
  ) } +
  labs(title = "NeoTCR8 score", color = "NeoTCR8") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_um_both <- (p_um4 + p_um8) +
  plot_layout(guides = "keep") +
  plot_annotation(
    title = paste0("NEW UMAP: NeoTCR4 & NeoTCR8 | VAL red; tumor_sc black | key=", key_use),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

save_plot_best_local(
  paste0(out_prefix, "_UMAP_NeoTCR4_NeoTCR8_scores_overlay_VALred_tumorBlack_key_", key_use),
  p_um_both, width = 13.5, height = 6.2, dpi = png_dpi
)

cd8a <- get_gene_vec(seu, "CD8A")
cd8b <- get_gene_vec(seu, "CD8B")
is_cd8 <- (cd8a > 0) | (cd8b > 0)

rank_df <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    clonotype = as.character(.data[[key_use]]),
    clonotype = ifelse(is.na(clonotype) | trimws(clonotype)=="", NA_character_, clonotype),
    lineage = ifelse(is_cd8[match(cell, colnames(seu))], "CD8", "CD4"),
    NeoTCR4 = as.numeric(NeoTCR4_score),
    NeoTCR8 = as.numeric(NeoTCR8_score),
    VAL_cell = as.logical(val_hit)
  ) %>%
  filter(!is.na(clonotype))

clono_tab <- rank_df %>%
  group_by(lineage, clonotype) %>%
  summarise(
    n_cells = n(),
    NeoTCR4_max = max(NeoTCR4, na.rm = TRUE),
    NeoTCR8_max = max(NeoTCR8, na.rm = TRUE),
    VAL_any = any(VAL_cell %in% TRUE),
    .groups = "drop"
  ) %>%
  mutate(NeoScore = ifelse(lineage=="CD4", NeoTCR4_max, NeoTCR8_max))

p90_cd4 <- quantile(clono_tab$NeoScore[clono_tab$lineage=="CD4"], probs = 0.9, na.rm = TRUE)
p90_cd8 <- quantile(clono_tab$NeoScore[clono_tab$lineage=="CD8"], probs = 0.9, na.rm = TRUE)

clono_tab <- clono_tab %>%
  mutate(in_top90 = case_when(
    lineage=="CD4" ~ NeoScore >= p90_cd4,
    lineage=="CD8" ~ NeoScore >= p90_cd8,
    TRUE ~ FALSE
  ))

pie_df <- clono_tab %>%
  filter(in_top90) %>%
  group_by(lineage) %>%
  summarise(n_top90 = n(), n_val = sum(VAL_any), n_nonval = n_top90 - n_val, .groups = "drop") %>%
  tidyr::pivot_longer(cols = c(n_val, n_nonval), names_to = "class", values_to = "n") %>%
  mutate(class = recode(class, n_val = "Tumor-specific (VAL)", n_nonval = "Non-tumor-specific"))

make_pie <- function(df, title_txt) {
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = title_txt))
  ggplot(df, aes(x = "", y = n, fill = class)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Tumor-specific (VAL)" = deep_red, "Non-tumor-specific" = "grey60")) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(face = "bold")
    ) +
    labs(title = title_txt)
}

rank_plot <- function(df, lineage, top_n = 80, p90_cd4, p90_cd8) {
  df2 <- df %>%
    dplyr::filter(lineage == !!lineage) %>%
    dplyr::arrange(dplyr::desc(NeoScore)) %>%
    dplyr::mutate(rank = dplyr::row_number())
  
  n_take <- min(top_n, nrow(df2))
  if (n_take <= 0) {
    return(ggplot() + theme_void() + labs(title = paste0(lineage, " ranking (empty)")))
  }
  
  df2 <- df2 %>%
    dplyr::slice(1:n_take) %>%
    dplyr::mutate(pt_col = ifelse(VAL_any, deep_red, "black"))
  
  thr <- if (lineage == "CD4") p90_cd4 else p90_cd8
  
  ggplot(df2, aes(x = rank, y = NeoScore)) +
    geom_hline(yintercept = thr, linetype = "dashed", linewidth = 0.7, color = "black") +
    annotate(
      "text",
      x = max(df2$rank), y = thr,
      label = "90th", hjust = 1.05, vjust = -0.4,
      fontface = "bold", size = 3.2, color = "black"
    ) +
    geom_line(color = "grey75", linewidth = 0.5) +
    geom_point(aes(color = pt_col), size = 2.2) +
    scale_color_identity() +
    labs(
      title = ifelse(
        lineage == "CD4",
        "Ranking of CD4+ TIL clonotypes by NeoTCR4 (max per clonotype)",
        "Ranking of CD8+ TIL clonotypes by NeoTCR8 (max per clonotype)"
      ),
      subtitle = "VAL highlighted red | 90th percentile dashed line",
      x = "Clonotype rank (desc score)",
      y = "Clonotype score (max over cells)"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

p_rank_cd4 <- rank_plot(clono_tab, "CD4", top_n = 80, p90_cd4 = p90_cd4, p90_cd8 = p90_cd8)
p_rank_cd8 <- rank_plot(clono_tab, "CD8", top_n = 80, p90_cd4 = p90_cd4, p90_cd8 = p90_cd8)
pie_cd4 <- make_pie(pie_df %>% filter(lineage=="CD4"), "CD4 top90: VAL fraction")
pie_cd8 <- make_pie(pie_df %>% filter(lineage=="CD8"), "CD8 top90: VAL fraction")

p_rank_all <- (p_rank_cd4 | pie_cd4) / (p_rank_cd8 | pie_cd8) +
  plot_layout(widths = c(2.8,1.2), heights = c(1,1)) +
  plot_annotation(
    title = paste0("HC25 clonotype ranking + top90 tumor-specific fraction (VAL red) | key=", key_use),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

save_plot_best_local(
  paste0(out_prefix, "_HC25_clonotypeRanking_NeoTCR4_CD4_NeoTCR8_CD8_withPie_top90_VALred_key_", key_use),
  p_rank_all, width = 13.5, height = 8.2, dpi = png_dpi
)

out_xlsx <- paste0(out_prefix, "_NeoTCR_scores_clonotypeRanking_key_", key_use, ".xlsx")
openxlsx::write.xlsx(
  list(
    per_cell = seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      select(cell, all_of(key_use), NeoTCR4_score, NeoTCR8_score, tumor_hit, val_hit),
    clonotype = clono_tab,
    pie_top90 = pie_df
  ),
  file = out_xlsx
)
message("[saved] ", out_xlsx)

## =========================================================
## 16) Final print
## =========================================================
cat("\nDONE.\n")
cat("Patient filter:\n")
cat(" - patient_col: ", patient_col, "\n")
cat(" - patient_keep: ", patient_keep, "\n")
cat("Saved key outputs:\n")
cat(" - Seurat (patient+functional filtered+reclustered): ", seu_rds_out, "\n")
cat(" - Marker heatmap: ", paste0(out_prefix, "_heatmap_markerPartitions_labelsAll_NEWclusters.pdf/.png/.svg"), "\n")
cat(" - UMAP selected TCR: ", paste0(out_prefix, "_umap_selectedTCR_allLabels_NEWUMAP_NEWclusters.pdf/.png/.svg"), "\n")
cat(" - UMAP functional.cluster: ", paste0(out_prefix, "_umap_functionalStateClusters_VALlabelsRed_NEWUMAP.pdf/.png/.svg"), "\n")
cat(" - UMAP new clusters + overlays: ", paste0(out_prefix, "_umap_NEWclusters_BGcolor_tumorBlack_VALred.pdf/.png/.svg"), "\n")
cat(" - Barplot: ", paste0(out_prefix, "_BAR_NEWcluster_tumorBlack_VALred_stacked.pdf/.png/.svg"), "\n")
cat(" - Excel (marker modules): ", wb_out1, "\n")
cat(" - Excel (cluster counts): ", out_xlsx2, "\n")
cat(" - Marker panels dir: ", out_merged, "\n")
cat(" - NeoTCR UMAP overlay: ", paste0(out_prefix, "_UMAP_NeoTCR4_NeoTCR8_scores_overlay_VALred_tumorBlack_key_", key_use, ".pdf/.png/.svg"), "\n")
cat(" - NeoTCR ranking Excel: ", out_xlsx, "\n")
cat("Resolution used for FindClusters: ", cluster_resolution, "\n")