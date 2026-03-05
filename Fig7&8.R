############################################################
## INVITRO EXP (in_vitro_exp_sc_CR==TRUE) — ONE MASTER SCRIPT
## Seurat v5 pipeline (CD4 + CD8) + downstream plots/tables
##
## Goals:
## - Single consolidated script (deduplicated helpers)
## - English only
## - Modular run switches
##
## Major outputs (under out_root):
## - CD4/ , CD8/ , shared_tables/
############################################################

## =========================
## 0) RUN SWITCHES
## =========================
RUN <- list(
  # Core
  cluster_cd4_cd8                = TRUE,   # gate + recluster + UMAP + save clustered RDS
  rawtype_vote_and_clusterTreeHM = TRUE,   # AddModuleScore + raw-type vote + ComplexHeatmap (cluster dendrogram)
  overlay_invitro_VALlegend      = TRUE,   # UMAP overlay: nonVAL single color; VAL colored by TCR label (legend) + export hits tables
  deg_export                     = TRUE,   # per-cluster DEG tables (FindAllMarkers)
  
  # Add-ons based on existing run folder outputs
  pie_VAL_perTCR_clusterComp     = TRUE,   # per TCR pie: cluster composition (VAL only)
  bar_VALcells_perCluster        = TRUE,   # per cluster bar: # VAL cells
  marker_panels_perCategory_RNA  = TRUE,   # per category marker panels (RNA only)
  percell_ssGSEA_heatmaps        = FALSE,  # per-cell ssGSEA heatmaps (large). Turn on if needed.
  qc_VAL_sanity                  = TRUE,   # QC sanity check for VAL cells (plots + tables)
  
  # NeoTCR (choose ONE flavor; you can enable multiple if you want)
  neotcr_simple_scores_only      = FALSE,  # compute NeoTCR4/8 and plot projections (no overlays)
  neotcr_overlay_pair_sc_plus_nv = FALSE,  # overlay sc_availability + nonvax_hits by (CDR3a,CDR3b) pair on UMAP
  cd8_3d_pillars_neotcr_pair     = FALSE   # CD8-only 3D pillars (heavy; can be slow)
)

## =========================
## 1) PATHS (EDIT HERE)
## =========================
rdata_path     <- "C:/Users/peppa/Downloads/20241107_Sorts_10X/Analysis_v3/seurat_filtered.RData"
xlsx_scAvail   <- "C:/Users/peppa/Downloads/sc_data_availibility.xlsx"
sig_cd4_xlsx   <- "C:/Users/peppa/Downloads/Xu_sig_cd4.xlsx"
sig_cd8_xlsx   <- "C:/Users/peppa/Downloads/xu_sig_cd8.xlsx"

## NeoTCR sources (optional modules)
sig_merged_xlsx     <- "C:/Users/peppa/Downloads/signature_xu_merged.xlsx"      # contains NeoTCR4/NeoTCR8 in row2
sig_merged_neo_xlsx <- "C:/Users/peppa/Downloads/signature_xu_merged_neo.xlsx"  # sheets NeoTCR4/NeoTCR8 (one gene column)

## Optional: extra overlay input (pair-match)
xlsx_nonvax_hits <- "C:/Users/peppa/Downloads/nonvax_hits_summary.xlsx"

## Output root (new folder per run)
out_dir <- "C:/Users/peppa/Downloads/"
run_id  <- paste0("INVITRO_CD4CD8_MASTER_", format(Sys.time(), "%Y%m%d_%H%M%S"))
out_root <- file.path(out_dir, run_id)

dir_cd4    <- file.path(out_root, "CD4")
dir_cd8    <- file.path(out_root, "CD8")
dir_shared <- file.path(out_root, "shared_tables")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_cd4, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_cd8, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_shared, showWarnings = FALSE, recursive = TRUE)

## If you want to run add-ons from an EXISTING run folder, set this (otherwise it uses out_root)
existing_run_root <- ""  # e.g. "C:/Users/peppa/Downloads/INVITRO_CD4CD8_20260222_171428"

## =========================
## 2) GLOBAL OPTIONS
## =========================
use_harmony_if_multi_patient <- TRUE
cluster_resolution <- 0.6
dims_use <- 1:30
png_dpi  <- 1200

deg_min_cells_cluster <- 100
deg_min_pct  <- 0.25
deg_logfc    <- 1
deg_test     <- "wilcox"

overlay_match_mode <- "contains"   # "contains" or "exact"
pair_match_mode    <- "exact"      # for (CDR3a,CDR3b) pair matching

## Overlay colors
col_hit_nonval <- "#1F2D3A"  # invitro hits but non-VAL
deep_red       <- "#8B0000"  # VAL

## Raw-type string scoring weights
W_NEOTCR     <- 6
W_NEOTCR_CD4 <- 7

## =========================
## 3) LIBRARIES (ONCE)
## =========================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(scales)
  library(Matrix)
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Missing package: patchwork. Install with install.packages('patchwork')")
  }
  library(patchwork)
  
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Missing package: ComplexHeatmap. Install with BiocManager::install('ComplexHeatmap')")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Missing package: circlize. Install with install.packages('circlize')")
  }
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## =========================
## 4) SHARED HELPERS (DEDUPED)
## =========================
save_plot_best <- function(filename_base, plot, width, height, dpi = 1200) {
  tryCatch({
    ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height,
           device = grDevices::cairo_pdf)
  }, error = function(e) {
    ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height,
           device = grDevices::pdf)
  })
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
           device = ragg::agg_png, dpi = dpi)
  } else {
    ggsave(paste0(filename_base, ".png"), plot, width = width, height = height, dpi = dpi)
  }
  
  tryCatch({
    ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height,
           device = grDevices::svg)
  }, error = function(e) {})
}

save_ht_best <- function(filename_base, ht, width, height, dpi = 1200) {
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
    ragg::agg_png(paste0(filename_base, ".png"), width = width, height = height, units = "in", res = dpi)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    grDevices::dev.off()
  } else {
    grDevices::png(paste0(filename_base, ".png"), width = width, height = height, units = "in", res = dpi)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    grDevices::dev.off()
  }
}

pick_first_existing <- function(cands, x) {
  hit <- cands[which(cands %in% x)][1]
  if (length(hit) == 0) NA_character_ else hit
}

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  x2 <- toupper(trimws(as.character(x)))
  x2 %in% c("TRUE","T","1","YES","Y")
}

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

auto_cut <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 50) return(NA_real_)
  med <- median(x); madv <- mad(x)
  max(med + 2 * madv, as.numeric(quantile(x, 0.80)))
}

get_adt_vec <- function(seu, feat, layer = "data") {
  if (!("ADT" %in% Assays(seu))) return(NULL)
  if (!(feat %in% rownames(seu[["ADT"]]))) return(NULL)
  DefaultAssay(seu) <- "ADT"
  FetchData(seu, vars = feat, layer = layer)[, 1]
}

get_first_seurat_in_env <- function() {
  objs <- mget(ls(envir = .GlobalEnv), envir = .GlobalEnv)
  idx <- which(sapply(objs, inherits, "Seurat"))
  if (length(idx) >= 1) return(objs[[idx[1]]])
  stop("No Seurat object found in the loaded RData.")
}

make_tcr_label <- function(x) {
  x <- as.character(x)
  n <- stringr::str_extract(x, "\\d+")
  ifelse(!is.na(n), paste0("TCR", n), paste0("TCR", x))
}

## =========================
## 5) LOAD INPUTS (Seurat + Excel)
## =========================
stopifnot(file.exists(rdata_path), file.exists(xlsx_scAvail))
load(rdata_path)

if (exists("filtered_seurat", envir = .GlobalEnv) && inherits(filtered_seurat, "Seurat")) {
  seu_all <- filtered_seurat
} else {
  seu_all <- get_first_seurat_in_env()
}
if (any(duplicated(colnames(seu_all)))) seu_all <- RenameCells(seu_all, new.names = make.unique(colnames(seu_all)))
DefaultAssay(seu_all) <- if ("RNA" %in% Assays(seu_all)) "RNA" else Assays(seu_all)[1]

patient_col <- pick_first_existing(
  c("patient","Patient","patient_id","PatientID","donor","Donor","sample","Sample","orig.ident"),
  colnames(seu_all@meta.data)
)

cdr3b_col_meta <- pick_first_existing(
  c("CDR3b","CDR3B","VDJ_cdr3_TRB","TRB_cdr3","cdr3b"),
  colnames(seu_all@meta.data)
)
if (is.na(cdr3b_col_meta)) {
  stop("Cannot find TRB CDR3 column in Seurat meta.data. Tried: CDR3b/CDR3B/VDJ_cdr3_TRB/TRB_cdr3/cdr3b")
}

dat_sc <- readxl::read_xlsx(xlsx_scAvail, .name_repair = "minimal")
need_cols <- c("TCR","VAL","in_vitro_exp_sc_CR","CDR3b")
miss <- setdiff(need_cols, colnames(dat_sc))
if (length(miss) > 0) stop("Excel missing columns: ", paste(miss, collapse=", "))

invitro_tab <- dat_sc %>%
  mutate(
    in_vitro_exp_sc_CR = to_logical(in_vitro_exp_sc_CR),
    VAL  = to_logical(VAL),
    TCR  = as.character(TCR),
    CDR3b = as.character(CDR3b)
  ) %>%
  filter(in_vitro_exp_sc_CR %in% TRUE) %>%
  filter(!is.na(TCR), TCR != "", !is.na(CDR3b), CDR3b != "") %>%
  mutate(
    TCR_label = make_tcr_label(TCR),
    CDR3b_u   = safe_str(CDR3b)
  )

write.xlsx(list(invitro_tab_used = invitro_tab),
           file = file.path(dir_shared, "invitro_tab_used.xlsx"))

message("[info] Seurat meta CDR3b column = ", cdr3b_col_meta)
message("[info] invitro_tab rows = ", nrow(invitro_tab), " | unique CDR3b = ", dplyr::n_distinct(invitro_tab$CDR3b_u))

## =========================
## 6) CD4/CD8 GATING + CLUSTERING
## =========================
gate_cd4_cd8 <- function(parent, want = c("CD4","CD8")) {
  want <- match.arg(want)
  have_adt <- "ADT" %in% Assays(parent)
  
  cd4_pos <- rep(FALSE, ncol(parent))
  cd8_pos <- rep(FALSE, ncol(parent))
  gate_source <- "ADT"
  
  if (have_adt) {
    adt_feats <- rownames(parent[["ADT"]])
    cd4_feat <- pick_first_existing(c("CD4","cd4","CD4_TotalSeqB","CD4-TotalSeqB","CD4ADT"), adt_feats)
    cd8_feat <- pick_first_existing(c("CD8","CD8A","cd8","CD8_TotalSeqB","CD8-TotalSeqB","CD8ADT","CD8a","CD8a_TotalSeqB"), adt_feats)
    
    if (!is.na(cd4_feat) && !is.na(cd8_feat)) {
      adt_cd4 <- get_adt_vec(parent, cd4_feat, layer="data")
      adt_cd8 <- get_adt_vec(parent, cd8_feat, layer="data")
      cd4_cut <- auto_cut(adt_cd4)
      cd8_cut <- auto_cut(adt_cd8)
      cd4_pos <- is.finite(cd4_cut) & (adt_cd4 > cd4_cut)
      cd8_pos <- is.finite(cd8_cut) & (adt_cd8 > cd8_cut)
    } else {
      have_adt <- FALSE
    }
  }
  
  if (!have_adt) {
    gate_source <- "RNA"
    if (!("CD4" %in% rownames(parent[["RNA"]]))) stop("RNA missing CD4 -> cannot gate")
    if (!("CD8A" %in% rownames(parent[["RNA"]])) && !("CD8B" %in% rownames(parent[["RNA"]]))) stop("RNA missing CD8A/B -> cannot gate")
    
    r_cd4 <- FetchData(parent, vars="CD4", layer="data")[,1]
    r_cd8 <- 0
    if ("CD8A" %in% rownames(parent[["RNA"]])) r_cd8 <- FetchData(parent, vars="CD8A", layer="data")[,1]
    if ("CD8B" %in% rownames(parent[["RNA"]])) r_cd8 <- pmax(r_cd8, FetchData(parent, vars="CD8B", layer="data")[,1])
    
    cd4_pos <- r_cd4 > 0
    cd8_pos <- r_cd8 > 0
  }
  
  cells <- if (want == "CD4") colnames(parent)[cd4_pos & !cd8_pos] else colnames(parent)[cd8_pos & !cd4_pos]
  list(cells = cells, gate_source = gate_source)
}

run_clustering <- function(obj, label, out_subdir) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, assay="RNA", verbose=FALSE)
  obj <- FindVariableFeatures(obj, assay="RNA", nfeatures=2000, verbose=FALSE)
  try({ obj[["RNA"]]@scale.data <- matrix(numeric(0), 0, 0) }, silent = TRUE)
  obj <- ScaleData(obj, assay="RNA", features=VariableFeatures(obj), verbose=FALSE)
  obj <- RunPCA(obj, assay="RNA", features=VariableFeatures(obj), npcs=50, verbose=FALSE)
  
  red <- "pca"
  if (use_harmony_if_multi_patient && !is.na(patient_col)) {
    n_pat <- dplyr::n_distinct(obj@meta.data[[patient_col]])
    if (n_pat >= 2 && requireNamespace("harmony", quietly = TRUE)) {
      suppressPackageStartupMessages(library(harmony))
      obj <- RunHarmony(obj, group.by.vars = patient_col, reduction="pca", assay.use="RNA", verbose=FALSE)
      red <- "harmony"
    }
  }
  
  obj <- FindNeighbors(obj, reduction=red, dims=dims_use, verbose=FALSE)
  obj <- FindClusters(obj, resolution=cluster_resolution, verbose=FALSE)
  obj <- RunUMAP(obj, reduction=red, dims=dims_use, verbose=FALSE)
  
  p <- DimPlot(obj, reduction="umap", group.by="seurat_clusters", label=TRUE, repel=TRUE) +
    ggtitle(paste0(label, " UMAP clusters (res=", cluster_resolution, ") | reduction=", red))
  save_plot_best(file.path(out_subdir, paste0(label, "_UMAP_clusters")), p, 8.8, 6.6, dpi = png_dpi)
  
  obj
}

## Clustered objects (will be assigned/loaded)
seu_cd4 <- NULL
seu_cd8 <- NULL

if (RUN$cluster_cd4_cd8) {
  gate_cd4 <- gate_cd4_cd8(seu_all, "CD4")
  seu_cd4 <- subset(seu_all, cells = gate_cd4$cells)
  message("[info] CD4 cells: ", ncol(seu_cd4), " | gate=", gate_cd4$gate_source)
  seu_cd4 <- run_clustering(seu_cd4, "CD4", dir_cd4)
  
  gate_cd8 <- gate_cd4_cd8(seu_all, "CD8")
  seu_cd8 <- subset(seu_all, cells = gate_cd8$cells)
  message("[info] CD8 cells: ", ncol(seu_cd8), " | gate=", gate_cd8$gate_source)
  seu_cd8 <- run_clustering(seu_cd8, "CD8", dir_cd8)
  
  saveRDS(seu_cd4, file = file.path(dir_cd4, "seu_cd4_clustered.rds"))
  saveRDS(seu_cd8, file = file.path(dir_cd8, "seu_cd8_clustered.rds"))
} else {
  ## Use existing run folder RDS if provided
  use_root <- if (nzchar(existing_run_root)) existing_run_root else out_root
  rds_cd4 <- file.path(use_root, "CD4", "seu_cd4_clustered.rds")
  rds_cd8 <- file.path(use_root, "CD8", "seu_cd8_clustered.rds")
  stopifnot(file.exists(rds_cd4), file.exists(rds_cd8))
  seu_cd4 <- readRDS(rds_cd4)
  seu_cd8 <- readRDS(rds_cd8)
}

## =========================
## 7) SIGNATURE SCORING + RAW TYPE VOTE + CLUSTER-TREE HEATMAP
## =========================
canon_type <- function(x) {
  s <- tolower(as.character(x))
  s <- gsub("[_/\\|]+", " ", s)
  s <- gsub("\\+", " + ", s)
  s <- gsub("[^a-z0-9\\+ ]+", " ", s)
  s <- gsub("\\s+", " ", s)
  trimws(s)
}

score_type_string_raw <- function(type_string) {
  s <- canon_type(type_string)
  toks <- unlist(strsplit(s, " ", fixed = TRUE))
  toks <- toks[toks != ""]
  sc <- 0
  add <- function(hit, w) { if (isTRUE(hit)) sc <<- sc + w }
  
  add(grepl("neotcr", s), W_NEOTCR)
  add(grepl("neotcr", s) && grepl("cd4", s), W_NEOTCR_CD4)
  
  add(("naive" %in% toks) || grepl("\\bnaive\\b", s), 3)
  add(("cm" %in% toks) || grepl("central memory", s), 3)
  add(("em" %in% toks) || grepl("\\bem\\b", s), 3)
  
  add(grepl("memory", s) || grepl("\\btcm\\b", s) || grepl("\\btem\\b", s), 2)
  add(grepl("memory effector", s) || grepl("memory/effector", s), 3)
  add(grepl("\\beffector\\b", s) || grepl("terminal", s) || grepl("cytotoxic", s), 3)
  add(("trm" %in% toks) || grepl("\\btrm\\b", s), 3)
  add(grepl("stem", s), 3)
  
  add(grepl("activated", s), 3)
  add(grepl("early", s) && grepl("activ", s), 3)
  
  add(("tex" %in% toks) || grepl("\\btex\\b", s), 3)
  add(grepl("exhaust", s), 3)
  add(grepl("progenitor", s) && grepl("exhaust", s), 3)
  add(("dys" %in% toks) || grepl("\\bdys\\b", s), 2)
  add(("tox" %in% toks) || grepl("\\btox\\b", s), 2)
  add(grepl("entpd1", s), 3)
  add(grepl("cxcl13", s), 3)
  
  add(grepl("prolif", s), 3)
  add(grepl("mitosis", s), 3)
  add(grepl("cell", s) && grepl("cycle", s), 3)
  add(grepl("hsp", s) || grepl("heat shock", s), 3)
  add(grepl("heat", s) && grepl("shock", s), 3)
  add(grepl("apopt", s), 3)
  
  add(grepl("\\btreg\\b", s), 3)
  add(grepl("treg", s) && grepl("like", s), 3)
  add(grepl("nr4a1", s) && grepl("treg", s), 3)
  add(grepl("il2rahi", s), 3)
  add(grepl("il2ralo", s), 3)
  add(grepl("\\btfh\\b", s), 3)
  add(grepl("\\bth\\b", s), 2)
  
  add(grepl("\\bmait\\b", s), 3)
  add(grepl("nk", s) && grepl("like", s), 3)
  add(grepl("innate", s) && grepl("like", s), 3)
  add(grepl("gamma", s) || grepl("delta", s) || grepl("gd", s), 3)
  add(grepl("klrb1", s), 3)
  
  add(grepl("gzmk", s), 3)
  add(grepl("gzmb", s), 3)
  add(grepl("fgfbp2", s), 3)
  add(grepl("\\bxcl\\b", s) || grepl("xcl1|xcl2", s), 3)
  add(grepl("mhc", s) && grepl("ii", s), 3)
  
  add(grepl("\\brpl\\b", s) || grepl("rpl32", s), 3)
  add(grepl("\\bmito\\b", s), 3)
  add(grepl("\\bchrom\\b", s), 3)
  add(grepl("\\btcf7\\b", s), 3)
  add(grepl("\\bfos\\b", s), 3)
  add(grepl("\\bil6st\\b", s), 3)
  
  if (nchar(trimws(as.character(type_string))) <= 10) sc <- sc + 2
  sc
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

make_type_order_from_signature <- function(sig_tbl_ordered) {
  sig_tbl_ordered %>%
    select(paper, state, col_i) %>%
    distinct() %>%
    arrange(col_i) %>%
    group_by(paper) %>%
    mutate(type_order_in_paper = row_number()) %>%
    ungroup() %>%
    mutate(
      paper = as.character(paper),
      state = as.character(state),
      paper_id = factor(paper, levels = unique(paper)),
      row_id = paste0(paper, "||", state)
    )
}

build_complex_type_heatmap <- function(prefix, out_subdir, cluster_top_perpaper_type_maxZ, type_order_df) {
  df <- cluster_top_perpaper_type_maxZ %>%
    mutate(
      cluster = as.character(cluster),
      paper = as.character(paper),
      state = as.character(state),
      row_id = paste0(paper, "||", state)
    ) %>%
    group_by(cluster, row_id) %>%
    summarise(z = max(z, na.rm = TRUE), .groups = "drop") %>%
    semi_join(type_order_df %>% select(row_id), by = "row_id")
  
  clus <- sort(unique(df$cluster))
  suppressWarnings({
    clus_num <- as.integer(clus)
    if (all(is.finite(clus_num))) clus <- as.character(sort(clus_num))
  })
  
  mat_wide <- df %>%
    mutate(cluster = factor(cluster, levels = clus)) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = z) %>%
    as.data.frame()
  
  rn <- mat_wide$row_id
  mat <- as.matrix(mat_wide[, setdiff(colnames(mat_wide), "row_id"), drop = FALSE])
  rownames(mat) <- rn
  colnames(mat) <- paste0("c", colnames(mat))
  
  row_ord <- type_order_df %>%
    filter(row_id %in% rownames(mat)) %>%
    arrange(paper_id, type_order_in_paper, row_id)
  
  mat <- mat[row_ord$row_id, , drop = FALSE]
  row_split <- factor(row_ord$paper, levels = levels(row_ord$paper_id))
  
  left_type_labels <- row_ord$state
  names(left_type_labels) <- row_ord$row_id
  
  ra_left <- rowAnnotation(
    Type = anno_text(left_type_labels, just = "left", gp = gpar(fontsize = 9)),
    width = unit(45, "mm"),
    show_annotation_name = FALSE
  )
  
  ra_right <- rowAnnotation(
    Paper = anno_block(
      gp = gpar(fill = NA, col = "black", lwd = 1.1),
      labels = levels(row_split),
      labels_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_rot = 0,
      labels_just = "center"
    ),
    width = unit(38, "mm"),
    show_annotation_name = FALSE
  )
  
  q <- as.numeric(stats::quantile(as.vector(mat), probs = c(0.02, 0.98), na.rm = TRUE))
  zlim <- if (all(is.finite(q)) && q[1] < q[2]) q else c(-2, 2)
  col_fun <- circlize::colorRamp2(c(zlim[1], 0, zlim[2]), c("#2C7BB6", "white", "#D7191C"))
  
  ht <- Heatmap(
    mat,
    name = "z",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = TRUE,
    column_dend_height = unit(25, "mm"),
    row_split = row_split,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    top_annotation = HeatmapAnnotation(
      Cluster = anno_text(colnames(mat), gp = gpar(fontsize = 9)),
      show_annotation_name = FALSE,
      annotation_height = unit(6, "mm")
    ),
    left_annotation = ra_left,
    right_annotation = ra_right,
    heatmap_legend_param = list(
      title = "z",
      legend_direction = "horizontal",
      title_position = "topcenter"
    )
  )
  
  w_in <- max(9, min(22, 3 + 0.35 * ncol(mat)))
  h_in <- max(7, min(28, 3 + 0.22 * nrow(mat)))
  
  save_ht_best(file.path(out_subdir, paste0(prefix, "_ComplexHeatmap_clusterTree_typesByPaperOrder")),
               ht, width = w_in, height = h_in, dpi = png_dpi)
  
  invisible(list(ht = ht, mat = mat, row_order = row_ord, row_split = row_split))
}

cross_label_invitro_rawVote <- function(obj, sig_xlsx, prefix, out_subdir) {
  sig_tbl_ordered <- read_signatures_with_order(sig_xlsx)
  
  DefaultAssay(obj) <- "RNA"
  feat <- toupper(rownames(obj[["RNA"]]))
  
  sig_tbl <- sig_tbl_ordered %>%
    mutate(
      genes_present = lapply(genes, function(g) intersect(g, feat)),
      n_present = vapply(genes_present, length, 1L),
      frac_present = ifelse(n_genes > 0, n_present/n_genes, 0)
    )
  
  min_present <- 8
  sig_use <- sig_tbl %>% filter(n_present >= min_present)
  if (nrow(sig_use) == 0) stop(prefix, ": no signatures with >= ", min_present, " genes present.")
  
  gene_list <- sig_use$genes_present
  names(gene_list) <- sig_use$sig_name
  
  score_name <- paste0("XuSig_", prefix, "_")
  obj <- AddModuleScore(obj, features = gene_list, assay="RNA", name = score_name)
  
  score_cols <- grep(paste0("^", score_name, "[0-9]+$"), colnames(obj@meta.data), value=TRUE)
  stopifnot(length(score_cols) == length(gene_list))
  
  score_map <- tibble(score_col = score_cols, sig_name = names(gene_list)) %>%
    left_join(sig_use %>% select(sig_name, paper, state, n_genes, n_present, frac_present), by="sig_name")
  
  md <- obj@meta.data %>%
    mutate(cluster=as.character(seurat_clusters)) %>%
    select(cluster, all_of(score_cols))
  
  cluster_sig_mean <- md %>%
    group_by(cluster) %>%
    summarise(across(all_of(score_cols), mean, na.rm=TRUE), .groups="drop")
  
  suppressWarnings({
    cc <- as.integer(cluster_sig_mean$cluster)
    if (all(is.finite(cc))) cluster_sig_mean <- cluster_sig_mean %>% arrange(cc)
  })
  
  cluster_sig_long <- cluster_sig_mean %>%
    pivot_longer(-cluster, names_to="score_col", values_to="mean_score") %>%
    left_join(score_map, by="score_col") %>%
    group_by(sig_name) %>%
    mutate(z = as.numeric(scale(mean_score))) %>%
    ungroup()
  
  ## same (cluster,paper,state): keep max z
  cluster_top_perpaper_type <- cluster_sig_long %>%
    group_by(cluster, paper, state) %>%
    slice_max(order_by=z, n=1, with_ties=FALSE) %>%
    ungroup()
  
  score_rows <- cluster_top_perpaper_type %>%
    mutate(type_raw = as.character(state)) %>%
    rowwise() %>%
    mutate(type_score = score_type_string_raw(type_raw)) %>%
    ungroup()
  
  cluster_type_vote <- score_rows %>%
    group_by(cluster, type_raw) %>%
    summarise(
      votes = n(),
      sum_score = sum(type_score, na.rm=TRUE),
      max_z = max(z, na.rm=TRUE),
      sum_z = sum(z, na.rm=TRUE),
      .groups="drop"
    ) %>%
    group_by(cluster) %>%
    arrange(desc(votes), desc(sum_score), desc(max_z), desc(sum_z), type_raw) %>%
    slice(1) %>%
    ungroup() %>%
    rename(final_type_label = type_raw)
  
  obj$final_type_label <- cluster_type_vote$final_type_label[
    match(as.character(obj$seurat_clusters), cluster_type_vote$cluster)
  ]
  
  ## cell best signature
  score_mat <- as.matrix(obj@meta.data[, score_cols, drop=FALSE])
  best_idx <- max.col(score_mat, ties.method="first")
  obj$cell_best_signature <- score_map$sig_name[match(score_cols[best_idx], score_map$score_col)]
  obj$cell_best_state_raw <- score_map$state[match(score_cols[best_idx], score_map$score_col)]
  
  ## heatmap order: paper order + within paper order from Excel column order
  type_order_df <- make_type_order_from_signature(sig_tbl_ordered)
  keep_pairs <- sig_use %>% distinct(paper, state)
  type_order_df <- type_order_df %>%
    inner_join(keep_pairs, by = c("paper","state")) %>%
    distinct(row_id, .keep_all = TRUE)
  
  ht_res <- build_complex_type_heatmap(
    prefix = prefix,
    out_subdir = out_subdir,
    cluster_top_perpaper_type_maxZ = cluster_top_perpaper_type %>% select(cluster, paper, state, z),
    type_order_df = type_order_df
  )
  
  wb <- file.path(out_subdir, paste0(prefix, "_signature_RAWtypeVote.xlsx"))
  write.xlsx(
    list(
      signatures_used = sig_use %>%
        transmute(paper, state, sig_name, n_genes, n_present,
                  frac_present = round(frac_present, 3),
                  genes_present = vapply(genes_present, function(g) paste(g, collapse=","), "")),
      score_map = score_map,
      cluster_top_perpaper_type_maxZ = cluster_top_perpaper_type,
      perpaper_type_scoring = score_rows %>% select(cluster, paper, state, z, type_score),
      cluster_type_vote_final = cluster_type_vote,
      type_order_from_excel = type_order_df %>% select(paper, state, col_i, type_order_in_paper),
      heatmap_matrix_used = as.data.frame(ht_res$mat) %>% tibble::rownames_to_column("paper||type")
    ),
    file = wb
  )
  
  p_final <- DimPlot(obj, reduction="umap", group.by="final_type_label", label=TRUE, repel=TRUE) +
    ggtitle(paste0(prefix, " UMAP: final_type_label (raw-type vote)"))
  save_plot_best(file.path(out_subdir, paste0(prefix, "_UMAP_finalTypeLabel_RAWvote")),
                 p_final, 9.0, 6.7, dpi=png_dpi)
  
  list(obj=obj, wb=wb)
}

if (RUN$rawtype_vote_and_clusterTreeHM) {
  stopifnot(file.exists(sig_cd4_xlsx), file.exists(sig_cd8_xlsx))
  res_cd4 <- cross_label_invitro_rawVote(seu_cd4, sig_cd4_xlsx, "CD4", dir_cd4)
  seu_cd4 <- res_cd4$obj
  res_cd8 <- cross_label_invitro_rawVote(seu_cd8, sig_cd8_xlsx, "CD8", dir_cd8)
  seu_cd8 <- res_cd8$obj
  
  saveRDS(seu_cd4, file=file.path(dir_cd4, "seu_cd4_withTypeLabels.rds"))
  saveRDS(seu_cd8, file=file.path(dir_cd8, "seu_cd8_withTypeLabels.rds"))
}

## =========================
## 8) OVERLAY: invitro hits + VAL legend (CDR3b match)
## =========================
invitro_overlay_VALlegend <- function(obj, prefix, out_subdir, invitro_tab, seurat_cdr3b_col,
                                      match_mode = c("contains","exact")) {
  match_mode <- match.arg(match_mode)
  md <- obj@meta.data
  if (!(seurat_cdr3b_col %in% colnames(md))) stop(prefix, ": meta missing TRB CDR3 column: ", seurat_cdr3b_col)
  
  meta_cdr3b <- safe_str(md[[seurat_cdr3b_col]])
  cells <- rownames(md)
  
  keys <- invitro_tab %>% select(TCR_label, VAL, CDR3b_u) %>% distinct()
  
  hits_all <- bind_rows(lapply(seq_len(nrow(keys)), function(i){
    b <- keys$CDR3b_u[i]
    idx <- if (match_mode == "exact") which(meta_cdr3b == b) else which(grepl(b, meta_cdr3b, fixed = TRUE))
    if (length(idx)==0) return(NULL)
    tibble(
      cell = cells[idx],
      TCR_label = keys$TCR_label[i],
      CDR3b = b,
      VAL = as.logical(keys$VAL[i])
    )
  })) %>% filter(!is.na(cell), cell!="")
  
  hits_one <- hits_all %>% group_by(cell) %>% slice(1) %>% ungroup()
  
  um <- as.data.frame(Embeddings(obj, "umap"))
  um <- um[,1:2,drop=FALSE]; colnames(um) <- c("UMAP_1","UMAP_2")
  um$cell <- rownames(um)
  um$cluster <- as.character(obj$seurat_clusters[match(um$cell, colnames(obj))])
  
  um$hit <- FALSE
  um$val_hit <- FALSE
  um$TCR_label <- NA_character_
  
  if (nrow(hits_one) > 0) {
    idxm <- match(hits_one$cell, um$cell)
    ok <- which(!is.na(idxm))
    um$hit[idxm[ok]] <- TRUE
    um$val_hit[idxm[ok]] <- hits_one$VAL[ok] %in% TRUE
    um$TCR_label[idxm[ok]] <- as.character(hits_one$TCR_label[ok])
  }
  
  df_nonval <- um %>% filter(hit, !val_hit)
  df_val    <- um %>% filter(hit, val_hit, !is.na(TCR_label), TCR_label!="")
  
  cent_cluster <- um %>%
    group_by(cluster) %>%
    summarise(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2), n=n(), .groups="drop") %>%
    mutate(lbl=paste0("c", cluster))
  
  val_levels <- df_val %>% distinct(TCR_label) %>% arrange(TCR_label) %>% pull(TCR_label)
  pal_val <- setNames(scales::hue_pal()(max(length(val_levels), 1)), val_levels)
  
  p <- ggplot(um, aes(UMAP_1, UMAP_2)) +
    geom_point(color="grey88", size=0.20) +
    geom_text(data=cent_cluster, aes(UMAP_1,UMAP_2,label=lbl),
              inherit.aes=FALSE, fontface="bold", size=4.0, color="black") +
    { if (nrow(df_nonval)>0) geom_point(data=df_nonval, inherit.aes=FALSE,
                                        aes(UMAP_1,UMAP_2),
                                        color=col_hit_nonval, size=1.10, alpha=0.95) } +
    { if (nrow(df_val)>0) geom_point(data=df_val, inherit.aes=FALSE,
                                     aes(UMAP_1,UMAP_2, color=TCR_label),
                                     size=1.55, alpha=0.98) } +
    { if (nrow(df_val)>0) scale_color_manual(values = pal_val, drop=FALSE) else NULL } +
    { if (nrow(df_val)>0) guides(color = guide_legend(override.aes = list(size = 3.5))) else guides(color="none") } +
    labs(
      title=paste0(prefix, " UMAP: invitro hits + VAL colored by TCR (legend)"),
      subtitle=paste0("Match=", match_mode,
                      " | meta CDR3b col=", seurat_cdr3b_col,
                      " | matched cells=", sum(um$hit),
                      " | VAL cells=", sum(um$val_hit)),
      x="UMAP_1", y="UMAP_2", color="VAL TCR"
    ) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"),
          legend.title=element_text(face="bold"))
  
  save_plot_best(file.path(out_subdir, paste0(prefix, "_UMAP_invitroHits_VALlegend")),
                 p, 9.2, 6.6, dpi=png_dpi)
  
  wb <- file.path(out_subdir, paste0(prefix, "_invitro_hits_VALlegend.xlsx"))
  write.xlsx(
    list(
      hits_all = hits_all,
      hits_onePerCell = hits_one,
      summary = tibble(
        prefix = prefix,
        seurat_cdr3b_col = seurat_cdr3b_col,
        match_mode = match_mode,
        n_cells_total = nrow(um),
        n_cells_hit = sum(um$hit),
        n_cells_val = sum(um$val_hit),
        n_val_tcr = dplyr::n_distinct(df_val$TCR_label)
      )
    ),
    file = wb
  )
  
  invisible(list(plot=p, hits_one=hits_one, wb=wb))
}

hits_cd4_xlsx <- file.path(dir_cd4, "CD4_invitro_hits_VALlegend.xlsx")
hits_cd8_xlsx <- file.path(dir_cd8, "CD8_invitro_hits_VALlegend.xlsx")

if (RUN$overlay_invitro_VALlegend) {
  stopifnot("umap" %in% Reductions(seu_cd4), "umap" %in% Reductions(seu_cd8))
  invitro_overlay_VALlegend(seu_cd4, "CD4", dir_cd4, invitro_tab, cdr3b_col_meta, match_mode = overlay_match_mode)
  invitro_overlay_VALlegend(seu_cd8, "CD8", dir_cd8, invitro_tab, cdr3b_col_meta, match_mode = overlay_match_mode)
}

## =========================
## 9) DEG EXPORT
## =========================
deg_export <- function(obj, prefix, out_subdir) {
  DefaultAssay(obj) <- "RNA"
  label_col <- if ("final_type_label" %in% colnames(obj@meta.data)) "final_type_label" else "seurat_clusters"
  
  clu_sizes <- table(obj$seurat_clusters)
  keep_clusters <- names(clu_sizes[clu_sizes >= deg_min_cells_cluster])
  
  Idents(obj) <- "seurat_clusters"
  markers_all <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    test.use = deg_test,
    min.pct = deg_min_pct,
    logfc.threshold = deg_logfc,
    verbose = FALSE
  ) %>%
    mutate(cluster = as.character(cluster)) %>%
    left_join(tibble(cluster=names(clu_sizes), n_cells=as.integer(clu_sizes)), by="cluster")
  
  markers_ge <- markers_all %>% filter(cluster %in% keep_clusters)
  
  cluster2label <- obj@meta.data %>%
    mutate(cluster=as.character(seurat_clusters), label=as.character(.data[[label_col]])) %>%
    distinct(cluster, label)
  
  markers_ge <- markers_ge %>% left_join(cluster2label, by="cluster")
  
  label_gene_support <- markers_ge %>%
    group_by(label, gene) %>%
    summarise(
      n_clusters_support = n_distinct(cluster),
      max_log2FC = max(avg_log2FC, na.rm=TRUE),
      min_p_adj = min(p_val_adj, na.rm=TRUE),
      .groups="drop"
    ) %>%
    arrange(label, desc(n_clusters_support), desc(max_log2FC))
  
  wb <- createWorkbook()
  addWorksheet(wb, "cluster_sizes")
  writeData(wb, "cluster_sizes",
            cluster2label %>%
              left_join(tibble(cluster=names(clu_sizes), n_cells=as.integer(clu_sizes)), by="cluster") %>%
              arrange(suppressWarnings(as.integer(cluster))))
  
  addWorksheet(wb, "markers_all_ge")
  writeData(wb, "markers_all_ge", markers_ge)
  
  addWorksheet(wb, "label_gene_support")
  writeData(wb, "label_gene_support", label_gene_support)
  
  clusters_sorted <- keep_clusters[order(suppressWarnings(as.integer(keep_clusters)))]
  for (cl in clusters_sorted) {
    df_one <- markers_ge %>% filter(cluster==cl) %>% arrange(p_val_adj, desc(avg_log2FC))
    lab <- unique(cluster2label$label[cluster2label$cluster==cl])
    lab <- ifelse(length(lab)==0, "NA", lab[1])
    sheet <- paste0("c", cl, "_", substr(lab, 1, 18))
    sheet <- str_replace_all(sheet, "[\\[\\]\\*\\?/\\\\:]", "_")
    sheet <- substr(sheet, 1, 31)
    addWorksheet(wb, sheet)
    writeData(wb, sheet, df_one)
  }
  
  out_xlsx <- file.path(out_subdir, paste0(prefix, "_perCluster_DEG.xlsx"))
  saveWorkbook(wb, out_xlsx, overwrite=TRUE)
  out_xlsx
}

if (RUN$deg_export) {
  deg_cd4 <- deg_export(seu_cd4, "CD4", dir_cd4)
  deg_cd8 <- deg_export(seu_cd8, "CD8", dir_cd8)
  message("[DEG] CD4: ", deg_cd4)
  message("[DEG] CD8: ", deg_cd8)
}

## =====================================================================
## FROM HERE: add-ons that typically use EXISTING run outputs
## If existing_run_root is set, use it; otherwise use out_root
## =====================================================================
use_root <- if (nzchar(existing_run_root)) existing_run_root else out_root
dir_cd4_use <- file.path(use_root, "CD4")
dir_cd8_use <- file.path(use_root, "CD8")
dir_shared_use <- file.path(use_root, "shared_tables")

rds_cd4_type <- file.path(dir_cd4_use, "seu_cd4_withTypeLabels.rds")
rds_cd8_type <- file.path(dir_cd8_use, "seu_cd8_withTypeLabels.rds")
hits_cd4_xlsx_use <- file.path(dir_cd4_use, "CD4_invitro_hits_VALlegend.xlsx")
hits_cd8_xlsx_use <- file.path(dir_cd8_use, "CD8_invitro_hits_VALlegend.xlsx")

## =========================
## 10) PIE: VAL-only per TCR cluster composition
## =========================
read_hits_val <- function(xlsx_path) {
  hits_one <- readxl::read_xlsx(xlsx_path, sheet = "hits_onePerCell") %>%
    mutate(cell = as.character(cell), TCR_label = as.character(TCR_label), VAL = as.logical(VAL))
  hits_one %>% filter(VAL %in% TRUE, !is.na(cell), cell != "", !is.na(TCR_label), TCR_label != "")
}

get_cluster_palette <- function(seu) {
  cl <- as.character(seu$seurat_clusters)
  cl[is.na(cl)] <- "NA"
  lv <- sort(unique(cl), na.last = TRUE)
  suppressWarnings({
    lv_num <- as.integer(lv)
    if (all(is.finite(lv_num))) lv <- as.character(sort(lv_num))
  })
  cluster_labels <- paste0("c", lv)
  pal <- setNames(scales::hue_pal()(max(length(cluster_labels), 1)), cluster_labels)
  list(levels = cluster_labels, pal = pal)
}

make_pie_df <- function(seu, hits_val, cluster_palette) {
  md <- seu@meta.data
  md$cell <- rownames(md)
  md$cluster <- as.character(md$seurat_clusters)
  md$cluster[is.na(md$cluster) | trimws(md$cluster) == ""] <- "NA"
  md$cluster_label <- paste0("c", md$cluster)
  
  df <- hits_val %>%
    select(cell, TCR_label) %>%
    distinct(cell, .keep_all = TRUE) %>%
    inner_join(md %>% select(cell, cluster_label), by = "cell") %>%
    mutate(cluster_label = factor(cluster_label, levels = cluster_palette$levels),
           TCR_label = as.character(TCR_label)) %>%
    count(TCR_label, cluster_label, name = "n") %>%
    group_by(TCR_label) %>%
    mutate(total = sum(n), prop = ifelse(total > 0, n / total, 0), pct = round(prop * 100)) %>%
    ungroup()
  
  df %>%
    tidyr::complete(
      TCR_label,
      cluster_label = factor(cluster_palette$levels, levels = cluster_palette$levels),
      fill = list(n = 0, total = NA_real_, prop = 0, pct = 0)
    ) %>%
    group_by(TCR_label) %>%
    mutate(total = ifelse(all(is.na(total)), sum(n), max(total, na.rm = TRUE)),
           prop = ifelse(total > 0, n / total, 0),
           pct = round(prop * 100)) %>%
    ungroup()
}

plot_pies <- function(pie_df, cluster_palette, title_prefix,
                      facet_ncol = 6, label_min_prop = 0.10) {
  if (nrow(pie_df) == 0) return(ggplot() + theme_void() + ggtitle(paste0(title_prefix, ": no VAL hits")))
  
  tcr_totals <- pie_df %>%
    group_by(TCR_label) %>%
    summarise(total = sum(n), .groups = "drop") %>%
    mutate(TCR_facet = paste0(TCR_label, "\n(n=", total, ")"))
  
  df <- pie_df %>%
    left_join(tcr_totals %>% select(TCR_label, TCR_facet), by = "TCR_label") %>%
    mutate(
      label_txt = dplyr::case_when(
        prop >= label_min_prop & n > 0 ~ paste0(as.character(cluster_label), "\n", "n=", n, "\n", pct, "%"),
        TRUE ~ ""
      )
    )
  
  tcr_levels <- tcr_totals %>%
    mutate(tcr_num = suppressWarnings(as.integer(str_extract(TCR_label, "\\d+")))) %>%
    arrange(ifelse(is.na(tcr_num), 1e9, tcr_num), TCR_label) %>%
    pull(TCR_label)
  
  df$TCR_label <- factor(df$TCR_label, levels = tcr_levels)
  df$TCR_facet <- factor(df$TCR_facet, levels = tcr_totals$TCR_facet[match(tcr_levels, tcr_totals$TCR_label)])
  
  ggplot(df %>% filter(prop > 0), aes(x = 1, y = prop, fill = cluster_label)) +
    geom_col(width = 1, color = "white", linewidth = 0.25) +
    geom_text(aes(label = label_txt), position = position_stack(vjust = 0.5),
              size = 2.9, fontface = "bold", color = "black") +
    coord_polar(theta = "y") +
    facet_wrap(~ TCR_facet, ncol = facet_ncol) +
    scale_fill_manual(values = cluster_palette$pal, drop = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = paste0(title_prefix, " | VAL-only per-TCR cluster composition (pie)"),
      subtitle = paste0("Each TCR normalized to 100%. Label if >= ", round(label_min_prop * 100), "%."),
      x = NULL, y = NULL, fill = "Cluster"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 9),
      panel.spacing = unit(0.35, "lines"),
      legend.position = "right"
    )
}

pie_size <- function(pie_df, facet_ncol = 6) {
  n_tcr <- pie_df %>% distinct(TCR_label) %>% nrow()
  w <- max(10, 1.6 * min(n_tcr, 18))
  h <- max(7, 0.95 * ceiling(max(n_tcr, 1) / facet_ncol) * 3.0)
  list(w = w, h = h)
}

if (RUN$pie_VAL_perTCR_clusterComp) {
  stopifnot(file.exists(rds_cd4_type), file.exists(rds_cd8_type),
            file.exists(hits_cd4_xlsx_use), file.exists(hits_cd8_xlsx_use))
  
  seu4 <- readRDS(rds_cd4_type)
  seu8 <- readRDS(rds_cd8_type)
  
  hits4 <- read_hits_val(hits_cd4_xlsx_use)
  hits8 <- read_hits_val(hits_cd8_xlsx_use)
  
  pal4 <- get_cluster_palette(seu4)
  pal8 <- get_cluster_palette(seu8)
  
  pie4 <- make_pie_df(seu4, hits4, pal4)
  pie8 <- make_pie_df(seu8, hits8, pal8)
  
  p_pie4 <- plot_pies(pie4, pal4, "CD4", facet_ncol = 6, label_min_prop = 0.10)
  p_pie8 <- plot_pies(pie8, pal8, "CD8", facet_ncol = 6, label_min_prop = 0.10)
  
  sz4 <- pie_size(pie4, facet_ncol = 6)
  sz8 <- pie_size(pie8, facet_ncol = 6)
  
  save_plot_best(file.path(dir_cd4_use, "CD4_PIE_VALonly_TCRclusterComposition"), p_pie4, sz4$w, sz4$h, dpi = png_dpi)
  save_plot_best(file.path(dir_cd8_use, "CD8_PIE_VALonly_TCRclusterComposition"), p_pie8, sz8$w, sz8$h, dpi = png_dpi)
  
  dir.create(file.path(use_root, "shared_tables"), showWarnings = FALSE, recursive = TRUE)
  write.xlsx(list(CD4_pie_table = pie4, CD8_pie_table = pie8),
             file = file.path(use_root, "shared_tables", "VALonly_TCR_cluster_pie_tables.xlsx"))
}

## =========================
## 11) BAR: # VAL cells per cluster (CD4/CD8)
## =========================
read_hits_val_onecell <- function(xlsx_path) {
  readxl::read_xlsx(xlsx_path, sheet = "hits_onePerCell") %>%
    mutate(cell = as.character(cell), VAL = as.logical(VAL)) %>%
    filter(VAL %in% TRUE, !is.na(cell), cell != "") %>%
    distinct(cell, .keep_all = TRUE)
}

cluster_levels_from_seu <- function(seu) {
  cl <- as.character(seu$seurat_clusters)
  cl[is.na(cl) | trimws(cl) == ""] <- "NA"
  cl_lab <- paste0("c", cl)
  
  ord <- tibble(cluster_label = cl_lab) %>%
    distinct() %>%
    mutate(
      cluster_num = suppressWarnings(as.integer(str_extract(cluster_label, "\\d+"))),
      cluster_num = ifelse(is.na(cluster_num), 1e9, cluster_num)
    ) %>%
    arrange(cluster_num, cluster_label) %>%
    pull(cluster_label)
  ord
}

make_cluster_palette <- function(cluster_levels) {
  setNames(scales::hue_pal()(max(length(cluster_levels), 1)), cluster_levels)
}

count_val_cells_by_cluster <- function(seu, hits_val) {
  md <- seu@meta.data
  md$cell <- rownames(md)
  md$cluster <- as.character(md$seurat_clusters)
  md$cluster[is.na(md$cluster) | trimws(md$cluster) == ""] <- "NA"
  md$cluster_label <- paste0("c", md$cluster)
  
  hits_val %>%
    select(cell) %>%
    inner_join(md %>% select(cell, cluster_label), by = "cell") %>%
    count(cluster_label, name = "n_val_cells")
}

make_bar_plot <- function(df_counts, cluster_levels, pal_cluster, title_prefix) {
  df_counts <- df_counts %>%
    mutate(cluster_label = factor(cluster_label, levels = cluster_levels)) %>%
    tidyr::complete(cluster_label = factor(cluster_levels, levels = cluster_levels),
                    fill = list(n_val_cells = 0)) %>%
    arrange(cluster_label)
  
  ggplot(df_counts, aes(x = cluster_label, y = n_val_cells, fill = cluster_label)) +
    geom_col(width = 0.85) +
    geom_text(aes(label = n_val_cells), vjust = -0.25, fontface = "bold", size = 3.0) +
    scale_fill_manual(values = pal_cluster, drop = FALSE, guide = "none") +
    labs(title = paste0(title_prefix, ": # VAL cells per cluster"),
         x = "Cluster", y = "# VAL cells") +
    theme_classic(base_size=11) +
    theme(plot.title = element_text(face="bold"),
          axis.text.x = element_text(size=8, face="bold")) +
    coord_cartesian(clip = "off")
}

if (RUN$bar_VALcells_perCluster) {
  stopifnot(file.exists(rds_cd4_type), file.exists(rds_cd8_type),
            file.exists(hits_cd4_xlsx_use), file.exists(hits_cd8_xlsx_use))
  
  seu4 <- readRDS(rds_cd4_type)
  seu8 <- readRDS(rds_cd8_type)
  hits4 <- read_hits_val_onecell(hits_cd4_xlsx_use)
  hits8 <- read_hits_val_onecell(hits_cd8_xlsx_use)
  
  lv4 <- cluster_levels_from_seu(seu4); pal4 <- make_cluster_palette(lv4)
  lv8 <- cluster_levels_from_seu(seu8); pal8 <- make_cluster_palette(lv8)
  
  c4 <- count_val_cells_by_cluster(seu4, hits4)
  c8 <- count_val_cells_by_cluster(seu8, hits8)
  
  p4 <- make_bar_plot(c4, lv4, pal4, "CD4")
  p8 <- make_bar_plot(c8, lv8, pal8, "CD8")
  
  save_plot_best(file.path(dir_cd4_use, "CD4_BAR_VALcells_perCluster_clusterColors"),
                 p4, width = max(8.5, 0.32 * length(lv4) + 4), height = 5.5, dpi = png_dpi)
  save_plot_best(file.path(dir_cd8_use, "CD8_BAR_VALcells_perCluster_clusterColors"),
                 p8, width = max(8.5, 0.32 * length(lv8) + 4), height = 5.5, dpi = png_dpi)
  
  dir.create(file.path(use_root, "shared_tables"), showWarnings = FALSE, recursive = TRUE)
  openxlsx::write.xlsx(
    list(
      CD4_VAL_cells_per_cluster = c4,
      CD4_cluster_palette = tibble(cluster_label = names(pal4), color = unname(pal4)),
      CD8_VAL_cells_per_cluster = c8,
      CD8_cluster_palette = tibble(cluster_label = names(pal8), color = unname(pal8))
    ),
    file = file.path(use_root, "shared_tables", "VAL_cells_per_cluster_with_palette.xlsx")
  )
}

## =========================
## 12) MARKER PANELS PER CATEGORY (RNA ONLY)
## =========================
if (RUN$marker_panels_perCategory_RNA) {
  suppressPackageStartupMessages({
    library(grid)
  })
  
  out_cd4 <- file.path(dir_cd4_use, "marker_panels_perCategory_RNA_only_noCD45_smallLegend")
  out_cd8 <- file.path(dir_cd8_use, "marker_panels_perCategory_RNA_only_noCD45_smallLegend")
  dir.create(out_cd4, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_cd8, showWarnings = FALSE, recursive = TRUE)
  
  LEG_H_MM <- 6
  LEG_W_MM <- 2.3
  LEG_TXT  <- 7
  
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
      "TCF7"  = c("TCF7"),
      "TOX"   = c("TOX"),
      "FOXP3" = c("FOXP3"),
      "PRDM1" = c("PRDM1"),
      "NR4A1" = c("NR4A1")
    )
    map[[key]]
  }
  
  spec_cd4 <- list(
    "Lineage and Naive/Memory markers" = list(items = list(
      list(key="CD4",  title="CD4 (RNA)"),
      list(key="CD8a", title="CD8A (RNA)")
    )),
    "Exhaustion/activation markers" = list(items = list(
      list(key="PD1",  title="PDCD1 (RNA)"),
      list(key="CD39", title="ENTPD1 (RNA)"),
      list(key="ICOS", title="ICOS (RNA)"),
      list(key="CD25", title="IL2RA (RNA)")
    )),
    "Effector molecules" = list(items = list(
      list(key="PRF1", title="PRF1 (RNA)"),
      list(key="IFNG", title="IFNG (RNA)"),
      list(key="GZMB", title="GZMB (RNA)"),
      list(key="GZMK", title="GZMK (RNA)")
    )),
    "Transcription factors" = list(items = list(
      list(key="MKI67", title="MKI67 (RNA)"),
      list(key="TCF7",  title="TCF7 (RNA)"),
      list(key="TOX",   title="TOX (RNA)"),
      list(key="FOXP3", title="FOXP3 (RNA)")
    ))
  )
  
  spec_cd8 <- list(
    "Lineage and Naive/Memory markers" = list(items = list(
      list(key="CD8a",  title="CD8A (RNA)"),
      list(key="IL7Aa", title="IL7R (RNA)"),
      list(key="CD28",  title="CD28 (RNA)"),
      list(key="CD4",   title="CD4 (RNA)"),
      list(key="CCR7",  title="CCR7 (RNA)"),
      list(key="CD62L", title="SELL (RNA)"),
      list(key="CD95",  title="FAS (RNA)"),
      list(key="CD27",  title="CD27 (RNA)")
    )),
    "Exhaustion/activation markers" = list(items = list(
      list(key="PD1",   title="PDCD1 (RNA)"),
      list(key="TIGIT", title="TIGIT (RNA)"),
      list(key="CD39",  title="ENTPD1 (RNA)"),
      list(key="ICOS",  title="ICOS (RNA)"),
      list(key="CD137", title="TNFRSF9 (RNA)")
    )),
    "Effector molecules" = list(items = list(
      list(key="PRF1",  title="PRF1 (RNA)"),
      list(key="FASLG", title="FASLG (RNA)"),
      list(key="GZMB",  title="GZMB (RNA)"),
      list(key="GZMH",  title="GZMH (RNA)"),
      list(key="GZMM",  title="GZMM (RNA)")
    )),
    "Transcription factors" = list(items = list(
      list(key="MKI67", title="MKI67 (RNA)"),
      list(key="TCF7",  title="TCF7 (RNA)"),
      list(key="TOX",   title="TOX (RNA)"),
      list(key="PRDM1", title="PRDM1 (RNA)"),
      list(key="NR4A1", title="NR4A1 (RNA)")
    ))
  )
  
  plot_one_category <- function(seu, category_name, items,
                                max_cols = 5,
                                pt_size = 0.18,
                                min_q = 0.01, max_q = 0.99) {
    stopifnot("umap" %in% Reductions(seu))
    DefaultAssay(seu) <- "RNA"
    
    emb <- as.data.frame(Embeddings(seu, "umap"))[, 1:2, drop = FALSE]
    colnames(emb) <- c("UMAP_1","UMAP_2")
    emb$cell <- rownames(emb)
    
    titles <- vapply(items, `[[`, "", "title")
    vals_list <- vector("list", length(items))
    
    for (i in seq_along(items)) {
      gene_cands <- resolve_to_rna_gene(items[[i]]$key)
      if (is.null(gene_cands)) {
        vals_list[[i]] <- rep(NA_real_, nrow(emb)); next
      }
      feat <- find_feature(seu, "RNA", gene_cands)
      if (is.na(feat)) {
        vals_list[[i]] <- rep(NA_real_, nrow(emb)); next
      }
      vdf <- tryCatch(FetchData(seu, vars = feat, layer = "data"), error = function(e) NULL)
      if (is.null(vdf)) vdf <- tryCatch(FetchData(seu, vars = feat), error = function(e) NULL)
      if (is.null(vdf) || nrow(vdf) == 0) {
        vals_list[[i]] <- rep(NA_real_, nrow(emb)); next
      }
      v <- as.numeric(vdf[, 1]); names(v) <- rownames(vdf)
      vals_list[[i]] <- v[match(emb$cell, names(v))]
    }
    
    allv <- unlist(vals_list, use.names = FALSE)
    allv <- allv[is.finite(allv)]
    if (length(allv) < 10) {
      lo <- 0; hi <- 1
    } else {
      lo <- as.numeric(quantile(allv, probs = min_q, na.rm = TRUE))
      hi <- as.numeric(quantile(allv, probs = max_q, na.rm = TRUE))
      if (!is.finite(lo) || !is.finite(hi) || lo >= hi) { lo <- 0; hi <- 1 }
    }
    
    plist <- vector("list", length(items))
    for (i in seq_along(items)) {
      df <- emb
      df$val <- vals_list[[i]]
      
      if (!any(is.finite(df$val))) {
        plist[[i]] <- ggplot() + theme_void() +
          labs(title = paste0(titles[i], "\n[MISSING]")) +
          theme(plot.title = element_text(hjust=0.5, face="bold", size=11))
        next
      }
      
      plist[[i]] <- ggplot(df, aes(UMAP_1, UMAP_2, color = val)) +
        geom_point(size = pt_size) +
        coord_equal() +
        scale_color_gradient(low = "grey90", high = "#2C5BFF",
                             limits = c(lo, hi), oob = scales::squish) +
        labs(title = titles[i], color = NULL) +
        theme_void(base_size = 10) +
        theme(
          plot.title = element_text(hjust=0.5, face="bold", size=11),
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
        theme = theme(plot.title = element_text(color = "#1E88E5", face="bold", size=14, hjust=0))
      )
    
    combined & theme(legend.position = "right", legend.box.just = "center", plot.margin = margin(3, 3, 3, 3))
  }
  
  make_dims <- function(n_panels, max_cols = 5) {
    ncol <- min(max_cols, n_panels)
    nrow <- ceiling(n_panels / ncol)
    width  <- max(6.8, ncol * 2.55 + 1.4)
    height <- max(4.2, nrow * 2.15 + 0.9)
    list(w=width, h=height)
  }
  
  run_all_panels <- function(seu, spec, out_dir, prefix) {
    for (catn in names(spec)) {
      items <- spec[[catn]]$items
      p <- plot_one_category(seu, catn, items, max_cols = 5)
      dims <- make_dims(length(items), max_cols = 5)
      fname <- paste0(prefix, "_", gsub("[^A-Za-z0-9]+", "_", catn), "_RNA")
      save_plot_best(file.path(out_dir, fname), p, dims$w, dims$h, dpi = png_dpi)
    }
  }
  
  stopifnot(file.exists(rds_cd4_type), file.exists(rds_cd8_type))
  seu4 <- readRDS(rds_cd4_type)
  seu8 <- readRDS(rds_cd8_type)
  stopifnot("umap" %in% Reductions(seu4), "umap" %in% Reductions(seu8))
  
  run_all_panels(seu4, spec_cd4, out_cd4, "CD4")
  run_all_panels(seu8, spec_cd8, out_cd8, "CD8")
}

## =========================
## 13) QC SANITY: VAL CELLS
## =========================
if (RUN$qc_VAL_sanity) {
  suppressPackageStartupMessages({
    library(patchwork)
  })
  
  add_basic_qc_metrics_if_missing <- function(seu) {
    DefaultAssay(seu) <- "RNA"
    if (!all(c("nCount_RNA","nFeature_RNA") %in% colnames(seu@meta.data))) {
      counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
      seu$nCount_RNA   <- Matrix::colSums(counts)
      seu$nFeature_RNA <- Matrix::colSums(counts > 0)
    }
    if (!("percent.mt" %in% colnames(seu@meta.data))) {
      feats <- rownames(seu[["RNA"]])
      mt <- grep("^MT-", feats, value = TRUE)
      if (length(mt) > 0) seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt)
    }
    if (!("percent.ribo" %in% colnames(seu@meta.data))) {
      feats <- rownames(seu[["RNA"]])
      rb <- grep("^RP[LS]", feats, value = TRUE)
      if (length(rb) > 0) seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = rb)
    }
    seu
  }
  
  make_val_subset <- function(seu_all, hits_val) {
    md <- seu_all@meta.data %>% tibble::rownames_to_column("cell")
    map <- hits_val %>% select(cell, TCR_label) %>% distinct()
    
    md2 <- md %>%
      left_join(map, by = "cell") %>%
      mutate(VAL_cell = !is.na(TCR_label) & TCR_label != "",
             VAL_TCR  = ifelse(is.na(TCR_label), NA_character_, TCR_label))
    
    seu_all$VAL_cell <- md2$VAL_cell[match(rownames(seu_all@meta.data), md2$cell)]
    seu_all$VAL_TCR  <- md2$VAL_TCR[match(rownames(seu_all@meta.data), md2$cell)]
    
    cells_val <- rownames(seu_all@meta.data)[which(seu_all$VAL_cell %in% TRUE)]
    seu_val <- subset(seu_all, cells = cells_val)
    
    list(seu_all = seu_all, seu_val = seu_val, md2 = md2)
  }
  
  qc_summary_tables <- function(seu_val, prefix) {
    md <- seu_val@meta.data %>% tibble::rownames_to_column("cell") %>%
      mutate(cluster = as.character(seurat_clusters), VAL_TCR = as.character(VAL_TCR))
    
    by_tcr <- md %>%
      group_by(VAL_TCR) %>%
      summarise(
        n_cells = n(),
        nCount_median = median(nCount_RNA, na.rm=TRUE),
        nCount_IQR = IQR(nCount_RNA, na.rm=TRUE),
        nFeature_median = median(nFeature_RNA, na.rm=TRUE),
        nFeature_IQR = IQR(nFeature_RNA, na.rm=TRUE),
        percent_mt_median = if ("percent.mt" %in% colnames(md)) median(percent.mt, na.rm=TRUE) else NA_real_,
        percent_ribo_median = if ("percent.ribo" %in% colnames(md)) median(percent.ribo, na.rm=TRUE) else NA_real_,
        .groups = "drop"
      ) %>% arrange(desc(n_cells), VAL_TCR)
    
    by_cluster <- md %>%
      group_by(cluster) %>%
      summarise(
        n_cells = n(),
        nCount_median = median(nCount_RNA, na.rm=TRUE),
        nFeature_median = median(nFeature_RNA, na.rm=TRUE),
        percent_mt_median = if ("percent.mt" %in% colnames(md)) median(percent.mt, na.rm=TRUE) else NA_real_,
        percent_ribo_median = if ("percent.ribo" %in% colnames(md)) median(percent.ribo, na.rm=TRUE) else NA_real_,
        .groups = "drop"
      ) %>% arrange(suppressWarnings(as.integer(cluster)))
    
    tcr_cluster <- md %>%
      count(VAL_TCR, cluster, name="n") %>%
      group_by(VAL_TCR) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup() %>%
      arrange(VAL_TCR, desc(n), suppressWarnings(as.integer(cluster)))
    
    val_cells_list <- md %>%
      select(cell, seurat_clusters, VAL_TCR, nCount_RNA, nFeature_RNA,
             dplyr::any_of(c("percent.mt","percent.ribo"))) %>%
      arrange(VAL_TCR, suppressWarnings(as.integer(seurat_clusters)), desc(nCount_RNA))
    
    list(by_tcr=by_tcr, by_cluster=by_cluster, tcr_cluster=tcr_cluster, val_cells_list=val_cells_list)
  }
  
  plot_qc_panels <- function(seu_all, seu_val, prefix) {
    DefaultAssay(seu_all) <- "RNA"
    DefaultAssay(seu_val) <- "RNA"
    
    vln_feats <- c("nCount_RNA","nFeature_RNA")
    if ("percent.mt" %in% colnames(seu_val@meta.data)) vln_feats <- c(vln_feats, "percent.mt")
    if ("percent.ribo" %in% colnames(seu_val@meta.data)) vln_feats <- c(vln_feats, "percent.ribo")
    
    p_vln_cluster <- VlnPlot(seu_val, features = vln_feats, group.by = "seurat_clusters",
                             pt.size = 0.05, ncol = min(4, length(vln_feats))) +
      plot_annotation(title = paste0(prefix, " VAL cells QC: violin by cluster"))
    
    p_vln_tcr <- VlnPlot(seu_val, features = c("nCount_RNA","nFeature_RNA"),
                         group.by = "VAL_TCR", pt.size = 0.05, ncol = 2) +
      plot_annotation(title = paste0(prefix, " VAL cells QC: violin by VAL_TCR")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    p_sc <- FeatureScatter(seu_val, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      ggtitle(paste0(prefix, " VAL cells: nCount_RNA vs nFeature_RNA"))
    
    if (!("umap" %in% Reductions(seu_all))) {
      p_umaps <- ggplot() + theme_void() + ggtitle(paste0(prefix, ": no UMAP found"))
    } else {
      if (!("VAL_cell" %in% colnames(seu_all@meta.data))) seu_all$VAL_cell <- FALSE
      seu_all$VAL_cell <- ifelse(is.na(seu_all$VAL_cell), FALSE, seu_all$VAL_cell)
      
      p_umap_bg <- DimPlot(seu_all, reduction="umap", group.by="seurat_clusters",
                           pt.size=0.1, label=FALSE) +
        ggtitle(paste0(prefix, " UMAP clusters (background)"))
      
      p_umap_val_count <- FeaturePlot(seu_val, features = "nCount_RNA", reduction="umap", pt.size=0.25) +
        ggtitle(paste0(prefix, " VAL UMAP: nCount_RNA"))
      
      p_umap_val_feat <- FeaturePlot(seu_val, features = "nFeature_RNA", reduction="umap", pt.size=0.25) +
        ggtitle(paste0(prefix, " VAL UMAP: nFeature_RNA"))
      
      p_umap_val_mt <- NULL
      if ("percent.mt" %in% colnames(seu_val@meta.data)) {
        p_umap_val_mt <- FeaturePlot(seu_val, features = "percent.mt", reduction="umap", pt.size=0.25) +
          ggtitle(paste0(prefix, " VAL UMAP: percent.mt"))
      }
      
      p_umaps <- if (is.null(p_umap_val_mt)) {
        (p_umap_bg / (p_umap_val_count | p_umap_val_feat))
      } else {
        (p_umap_bg / (p_umap_val_count | p_umap_val_feat | p_umap_val_mt))
      }
    }
    
    list(p_vln_cluster=p_vln_cluster, p_vln_tcr=p_vln_tcr, p_sc=p_sc, p_umaps=p_umaps)
  }
  
  run_val_qc <- function(rds_path, hits_xlsx, prefix, out_dir, dpi = 1200) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    seu_all <- readRDS(rds_path)
    hits_val <- read_hits_val_onecell(hits_xlsx)
    seu_all <- add_basic_qc_metrics_if_missing(seu_all)
    
    tmp <- make_val_subset(seu_all, hits_val)
    seu_all <- tmp$seu_all
    seu_val <- tmp$seu_val
    
    if (ncol(seu_val) == 0) stop(prefix, ": no VAL cells after matching. Check hits table.")
    
    qs <- qc_summary_tables(seu_val, prefix)
    ps <- plot_qc_panels(seu_all, seu_val, prefix)
    
    n_tcr <- dplyr::n_distinct(seu_val$VAL_TCR)
    w_tcr <- max(10, min(40, 6 + 0.25 * n_tcr))
    
    save_plot_best(file.path(out_dir, paste0(prefix, "_QC_VALcells_violin_byCluster")),
                   ps$p_vln_cluster, width = 12, height = 7, dpi = dpi)
    save_plot_best(file.path(out_dir, paste0(prefix, "_QC_VALcells_violin_byTCR")),
                   ps$p_vln_tcr, width = w_tcr, height = 7, dpi = dpi)
    save_plot_best(file.path(out_dir, paste0(prefix, "_QC_VALcells_scatter_nCount_vs_nFeature")),
                   ps$p_sc, width = 7.2, height = 6.2, dpi = dpi)
    save_plot_best(file.path(out_dir, paste0(prefix, "_QC_VALcells_UMAP_QCmaps")),
                   ps$p_umaps, width = 14, height = 10, dpi = dpi)
    
    out_xlsx <- file.path(out_dir, paste0(prefix, "_QC_VALcells_summary.xlsx"))
    openxlsx::write.xlsx(
      list(
        VAL_hits_onePerCell = hits_val,
        VAL_cells_list = qs$val_cells_list,
        VAL_TCR_x_cluster = qs$tcr_cluster,
        QC_by_TCR = qs$by_tcr,
        QC_by_cluster = qs$by_cluster
      ),
      file = out_xlsx
    )
    out_xlsx
  }
  
  stopifnot(file.exists(rds_cd4_type), file.exists(rds_cd8_type),
            file.exists(hits_cd4_xlsx_use), file.exists(hits_cd8_xlsx_use))
  
  out_cd4_qc <- file.path(dir_cd4_use, "QC_VALcells_sanitycheck")
  out_cd8_qc <- file.path(dir_cd8_use, "QC_VALcells_sanitycheck")
  
  qc_cd4_xlsx <- run_val_qc(rds_cd4_type, hits_cd4_xlsx_use, "CD4", out_cd4_qc, dpi = png_dpi)
  qc_cd8_xlsx <- run_val_qc(rds_cd8_type, hits_cd8_xlsx_use, "CD8", out_cd8_qc, dpi = png_dpi)
  
  message("[QC] CD4: ", qc_cd4_xlsx)
  message("[QC] CD8: ", qc_cd8_xlsx)
}

############################################################
## FINAL SUMMARY
############################################################
cat("\nDONE.\n")
cat("Output root:\n - ", out_root, "\n", sep = "")
cat("Primary outputs:\n")
cat(" - CD4 clustered RDS: ", file.path(dir_cd4, "seu_cd4_clustered.rds"), "\n", sep = "")
cat(" - CD8 clustered RDS: ", file.path(dir_cd8, "seu_cd8_clustered.rds"), "\n", sep = "")
if (file.exists(file.path(dir_cd4, "seu_cd4_withTypeLabels.rds"))) {
  cat(" - CD4 type labels RDS: ", file.path(dir_cd4, "seu_cd4_withTypeLabels.rds"), "\n", sep = "")
}
if (file.exists(file.path(dir_cd8, "seu_cd8_withTypeLabels.rds"))) {
  cat(" - CD8 type labels RDS: ", file.path(dir_cd8, "seu_cd8_withTypeLabels.rds"), "\n", sep = "")
}
cat(" - invitro_tab_used.xlsx: ", file.path(dir_shared, "invitro_tab_used.xlsx"), "\n", sep = "")
cat(" - Overlay meta CDR3b col: ", cdr3b_col_meta, " | match mode: ", overlay_match_mode, "\n", sep = "")