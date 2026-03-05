suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  library(readxl)
  library(scales)
  library(ragg)
  library(plotly)
  library(htmlwidgets)
})

setwd("C:/Users/peppa/Downloads/NonVax_2026")

###############################################
## GLOBAL FIGURE SIZE
###############################################
FIG_W_MM <- 190
FIG_H_MM <- 140
FIG_DPI  <- 600

###############################################
## PATHS
###############################################
base_dir <- "C:/Users/peppa/Downloads/NonVax_2026"

cd4_xlsx <- file.path(base_dir, "cd4.xlsx")
cd8_xlsx <- file.path(base_dir, "cd8.xlsx")

nonvax_hits_path <- "C:/Users/peppa/Downloads/nonvax_hits.xlsx"

sc_avail_path <- "C:/Users/peppa/Downloads/sc_data_availibility.xlsx"
tcr_map_path  <- "C:/Users/peppa/Downloads/tcr_to_guide_id_matches.xlsx"

###############################################
## COLORS
###############################################
A_GREEN <- "#2E7D32"
B_RED   <- "#C62828"

BG_A <- "#DFF0D8"
BG_B <- "#F8D7DA"

BASE_GREY <- "grey80"
SIG_GOLD  <- scales::alpha("gold", 0.35)

VAX_FILL_GREY <- "grey75"
VAX_NONVAL_BORDER <- "black"
VAX_VAL_BORDER    <- B_RED

###############################################
## HELPERS
###############################################
clean_names_simple <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", "_", nms)
  nms <- gsub("^_|_$", "", nms)
  tolower(nms)
}

clean_names_unique <- function(nms) {
  nms2 <- clean_names_simple(nms)
  make.unique(nms2, sep = "_dup")
}

to_logical <- function(v) {
  if (is.logical(v)) return(v)
  if (is.numeric(v)) return(v != 0)
  if (is.factor(v))  return(to_logical(as.character(v)))
  if (is.character(v)) {
    x <- toupper(trimws(v))
    return(x %in% c("TRUE", "T", "1", "YES", "Y"))
  }
  as.logical(v)
}

pick_first <- function(nms, pattern_vec) {
  for (pat in pattern_vec) {
    idx <- which(grepl(pat, nms, ignore.case = TRUE))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

make_unique_names <- function(x) {
  make.unique(x, sep = "_dup")
}

safe_num_max_abs <- function(x, fallback = 1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(fallback)
  max(abs(x), na.rm = TRUE)
}

safe_num_max <- function(x, fallback = 1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(fallback)
  max(x, na.rm = TRUE)
}

nice_limit_sym <- function(x, min_limit = 1) {
  v <- safe_num_max_abs(x, fallback = min_limit)
  v <- max(v, min_limit)
  ceiling(v * 10) / 10
}

nice_limit_pos <- function(x, min_limit = 1) {
  v <- safe_num_max(x, fallback = min_limit)
  v <- max(v, min_limit)
  ceiling(v * 10) / 10
}

save_gg <- function(p, outdir, fname,
                    width_mm = FIG_W_MM,
                    height_mm = FIG_H_MM,
                    dpi = FIG_DPI) {
  ggsave(
    filename = file.path(outdir, paste0(fname, ".pdf")),
    plot = p,
    width = width_mm, height = height_mm, units = "mm",
    device = grDevices::cairo_pdf
  )
  ggsave(
    filename = file.path(outdir, paste0(fname, ".png")),
    plot = p,
    width = width_mm, height = height_mm, units = "mm",
    dpi = dpi,
    device = ragg::agg_png,
    bg = "white"
  )
}

save_plotly <- function(p_gg, outdir, fname) {
  p_int <- plotly::ggplotly(p_gg, tooltip = "text") %>%
    plotly::layout(hoverlabel = list(align = "left"))
  htmlwidgets::saveWidget(
    p_int,
    file = file.path(outdir, paste0(fname, ".html")),
    selfcontained = TRUE
  )
  invisible(p_int)
}

###############################################
## BACKGROUND LAYERS
###############################################
ma_bg_layers <- function(x_mid, y_lim = c(-2, 2),
                         fill_pos = BG_A, fill_neg = BG_B,
                         rect_alpha = 0.45,
                         text_alpha = 0.92,
                         text_size = 34) {
  list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
             fill = fill_pos, alpha = rect_alpha),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = fill_neg, alpha = rect_alpha),
    annotate("text", x = x_mid, y = y_lim[2] * 0.65, label = "A",
             color = "white", fontface = "bold", size = text_size, alpha = text_alpha),
    annotate("text", x = x_mid, y = y_lim[1] * 0.65, label = "B",
             color = "white", fontface = "bold", size = text_size, alpha = text_alpha)
  )
}

volcano_bg_layers <- function(x_lim = c(-3, 3), y_lim = c(0, 20),
                              fill_pos = BG_A, fill_neg = BG_B,
                              rect_alpha = 0.45,
                              text_alpha = 0.92,
                              text_size = 34) {
  xA <- (0 + x_lim[2]) / 2
  xB <- (x_lim[1] + 0) / 2
  yT <- y_lim[2] * 0.70
  
  list(
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = fill_pos, alpha = rect_alpha),
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
             fill = fill_neg, alpha = rect_alpha),
    annotate("text", x = xA, y = yT, label = "A",
             color = "white", fontface = "bold", size = text_size, alpha = text_alpha),
    annotate("text", x = xB, y = yT, label = "B",
             color = "white", fontface = "bold", size = text_size, alpha = text_alpha)
  )
}

###############################################
## READ COUNTS
###############################################
read_count_xlsx <- function(path) {
  df <- read_excel(path, col_names = TRUE) %>% as.data.frame()
  colnames(df) <- make_unique_names(colnames(df))
  colnames(df)[1] <- "guide_id"
  
  df <- df %>%
    filter(!is.na(guide_id), guide_id != "")
  
  if (anyDuplicated(df$guide_id) > 0) {
    dup_ids <- unique(df$guide_id[duplicated(df$guide_id)])
    warning(
      "Duplicated guide_id found in ", basename(path), ": ",
      paste(head(dup_ids, 10), collapse = ", "),
      if (length(dup_ids) > 10) " ..." else ""
    )
    
    df <- df %>%
      group_by(guide_id) %>%
      summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  }
  
  sample_cols <- setdiff(colnames(df), "guide_id")
  df[, sample_cols] <- lapply(df[, sample_cols, drop = FALSE], as.numeric)
  df[, sample_cols] <- lapply(df[, sample_cols, drop = FALSE], function(x) {
    x[is.na(x)] <- 0
    x
  })
  
  df
}

make_count_and_meta <- function(df, label_for_error = "dataset") {
  count_data <- df %>%
    column_to_rownames("guide_id") %>%
    as.matrix()
  
  count_data <- round(count_data)
  storage.mode(count_data) <- "integer"
  
  sample_names <- colnames(count_data)
  
  condition <- case_when(
    grepl("_A_", sample_names) ~ "A",
    grepl("_B_", sample_names) ~ "B",
    TRUE ~ NA_character_
  )
  
  if (any(is.na(condition))) {
    stop(
      "Cannot assign A/B condition for these columns in ", label_for_error, ": ",
      paste(sample_names[is.na(condition)], collapse = ", ")
    )
  }
  
  ## positive log2FC = A enrichment
  meta <- data.frame(
    row.names = sample_names,
    condition = factor(condition, levels = c("B", "A"))
  )
  
  stopifnot(all(colnames(count_data) == rownames(meta)))
  
  list(count_data = count_data, meta = meta)
}

###############################################
## READ OLD VAX GUIDES
###############################################
read_old_vax_guides <- function(cell = c("cd4", "cd8"),
                                sc_avail_path,
                                tcr_map_path) {
  cell <- match.arg(cell)
  
  sc_avail <- read_excel(sc_avail_path) %>% as.data.frame()
  colnames(sc_avail) <- clean_names_unique(colnames(sc_avail))
  
  tcr_map <- read_excel(tcr_map_path) %>% as.data.frame()
  colnames(tcr_map) <- clean_names_unique(colnames(tcr_map))
  
  sc_cell_col <- if (cell %in% names(sc_avail)) {
    cell
  } else {
    pick_first(names(sc_avail), c(paste0("^", cell, "$"), cell))
  }
  
  sc_val_col <- if ("val" %in% names(sc_avail)) {
    "val"
  } else {
    pick_first(names(sc_avail), c("^val$", "val"))
  }
  
  sc_tcr_col <- if ("tcr" %in% names(sc_avail)) {
    "tcr"
  } else {
    pick_first(names(sc_avail), c("^tcr$", "tcr"))
  }
  
  map_tcr_col <- if ("tcr" %in% names(tcr_map)) {
    "tcr"
  } else {
    pick_first(names(tcr_map), c("^tcr$", "tcr"))
  }
  
  map_guide_col <- if ("guide_id" %in% names(tcr_map)) {
    "guide_id"
  } else {
    pick_first(names(tcr_map), c("^guide_id$", "guide", "guideid"))
  }
  
  if (any(is.na(c(sc_cell_col, sc_val_col, sc_tcr_col)))) {
    stop(
      "sc_data_availibility.xlsx 缺少必要列（", toupper(cell), " / VAL / TCR）。现有列名：\n",
      paste(names(sc_avail), collapse = ", ")
    )
  }
  
  if (any(is.na(c(map_tcr_col, map_guide_col)))) {
    stop(
      "tcr_to_guide_id_matches.xlsx 缺少必要列（TCR / guide_id）。现有列名：\n",
      paste(names(tcr_map), collapse = ", ")
    )
  }
  
  sc_avail <- sc_avail %>%
    mutate(
      cell_flag = to_logical(.data[[sc_cell_col]]),
      val_flag  = to_logical(.data[[sc_val_col]]),
      tcr       = as.character(.data[[sc_tcr_col]])
    )
  
  ## keep your old special rule
  if (cell == "cd4") {
    sc_avail <- sc_avail %>% mutate(val_flag = ifelse(tcr == "18", TRUE, val_flag))
  } else if (cell == "cd8") {
    sc_avail <- sc_avail %>% mutate(val_flag = ifelse(tcr == "18", FALSE, val_flag))
  }
  
  tcr_map <- tcr_map %>%
    transmute(
      tcr      = as.character(.data[[map_tcr_col]]),
      guide_id = as.character(.data[[map_guide_col]])
    ) %>%
    filter(!is.na(tcr), !is.na(guide_id), tcr != "", guide_id != "") %>%
    distinct(tcr, guide_id)
  
  sc_avail %>%
    filter(cell_flag) %>%
    distinct(tcr, val_flag) %>%
    left_join(tcr_map, by = "tcr") %>%
    filter(!is.na(guide_id), guide_id != "") %>%
    distinct(guide_id, .keep_all = TRUE) %>%
    mutate(
      vax_group = ifelse(isTRUE(val_flag), "VAL", "NONVAL")
    )
}

###############################################
## READ NONVAX HITS
###############################################
read_nonvax_hits <- function(path, cell = c("cd4", "cd8")) {
  cell <- match.arg(cell)
  
  df <- read_excel(path) %>% as.data.frame()
  colnames(df) <- clean_names_unique(colnames(df))
  
  guide_col <- if ("guide_id" %in% names(df)) {
    "guide_id"
  } else {
    pick_first(names(df), c("^guide_id$", "guide", "guideid"))
  }
  
  nv_col <- if ("nv_tcr" %in% names(df)) {
    "nv_tcr"
  } else {
    pick_first(names(df), c("^nv_tcr$", "nv_tcr", "^tcr$"))
  }
  
  cell_col <- if (cell %in% names(df)) {
    cell
  } else {
    pick_first(names(df), c(paste0("^", cell, "$"), cell))
  }
  
  df %>%
    transmute(
      guide_id = as.character(.data[[guide_col]]),
      nv_tcr   = as.character(.data[[nv_col]]),
      cell_ok  = to_logical(.data[[cell_col]])
    ) %>%
    filter(cell_ok, !is.na(guide_id), guide_id != "", !is.na(nv_tcr), nv_tcr != "") %>%
    distinct(guide_id, .keep_all = TRUE) %>%
    mutate(label = paste0("NV_TCR", nv_tcr))
}

###############################################
## MAIN PIPELINE
###############################################
run_nonvax_pipeline <- function(cell = c("cd4", "cd8"),
                                count_xlsx,
                                nonvax_hits_path,
                                sc_avail_path,
                                tcr_map_path,
                                out_root = base_dir) {
  
  cell <- match.arg(cell)
  CELL <- toupper(cell)
  
  message("==== Running NONVAX: ", CELL, " ====")
  
  ###############################################
  ## 1. counts and DESeq2
  ###############################################
  raw_df <- read_count_xlsx(count_xlsx)
  x <- make_count_and_meta(raw_df, label_for_error = CELL)
  count_data <- x$count_data
  meta <- x$meta
  
  message(CELL, " dimensions: ", nrow(count_data), " guides x ", ncol(count_data), " samples")
  message(CELL, " condition table:")
  print(table(meta$condition))
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = meta,
    design    = ~ condition
  )
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("condition", "A", "B"),
    alpha = 0.05,
    independentFiltering = FALSE
  )
  
  res_tb <- res %>%
    data.frame() %>%
    rownames_to_column("guide_id") %>%
    as_tibble() %>%
    mutate(
      guide_id = as.character(guide_id),
      padj0.05 = case_when(
        is.na(padj)  ~ "NA",
        padj <= 0.05 ~ "YES",
        TRUE         ~ "NO"
      ),
      logfc1 = case_when(
        is.na(log2FoldChange)    ~ "NA",
        abs(log2FoldChange) >= 1 ~ "YES",
        TRUE                     ~ "NO"
      ),
      logfc0.5 = case_when(
        is.na(log2FoldChange)      ~ "NA",
        abs(log2FoldChange) >= 0.5 ~ "YES",
        TRUE                       ~ "NO"
      )
    )
  
  write_excel_csv(
    res_tb,
    file.path(out_root, paste0(CELL, "_NONVAX_DESeq2_results_A_vs_B_annotated.csv"))
  )
  
  hits_un <- res_tb %>%
    filter(!is.na(padj), padj < 0.05)
  
  ###############################################
  ## 2. annotation tables
  ###############################################
  nonvax_hits <- read_nonvax_hits(nonvax_hits_path, cell = cell)
  old_vax_tbl <- read_old_vax_guides(
    cell = cell,
    sc_avail_path = sc_avail_path,
    tcr_map_path  = tcr_map_path
  )
  
  message(CELL, " nonvax hits used: ", nrow(nonvax_hits))
  print(nonvax_hits)
  
  hit_stats <- nonvax_hits %>%
    left_join(res_tb, by = "guide_id") %>%
    mutate(
      hit_side = case_when(
        !is.na(log2FoldChange) & log2FoldChange > 0 ~ "A",
        !is.na(log2FoldChange) & log2FoldChange < 0 ~ "B",
        TRUE ~ NA_character_
      ),
      hit_color = case_when(
        hit_side == "A" ~ A_GREEN,
        hit_side == "B" ~ B_RED,
        TRUE ~ "black"
      )
    )
  
  old_vax_stats <- old_vax_tbl %>%
    left_join(res_tb, by = "guide_id")
  
  ###############################################
  ## 3. plot-ready data
  ###############################################
  res_ma <- res_tb %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      x_ma = log10(baseMean + 1)
    )
  
  hits_un_ma <- hits_un %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      x_ma = log10(baseMean + 1)
    )
  
  hit_stats_ma <- hit_stats %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      x_ma = log10(baseMean + 1)
    )
  
  old_vax_ma <- old_vax_stats %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      x_ma = log10(baseMean + 1)
    )
  
  volc_df <- res_tb %>%
    filter(!is.na(padj), is.finite(padj), padj > 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      neglog10padj = -log10(padj)
    )
  
  hits_un_v <- hits_un %>%
    filter(!is.na(padj), is.finite(padj), padj > 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      neglog10padj = -log10(padj)
    )
  
  hit_stats_v <- hit_stats %>%
    filter(!is.na(padj), is.finite(padj), padj > 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      neglog10padj = -log10(padj)
    )
  
  old_vax_v <- old_vax_stats %>%
    filter(!is.na(padj), is.finite(padj), padj > 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      neglog10padj = -log10(padj)
    )
  
  ###############################################
  ## 4. limits
  ###############################################
  ###############################################
  ## 4. limits
  ###############################################
  ma_y_lim <- nice_limit_sym(res_ma$log2FoldChange, min_limit = 1)
  y_lim_ma <- c(-ma_y_lim, ma_y_lim)
  x_mid_ma <- median(res_ma$x_ma, na.rm = TRUE)
  
  ma_x_min <- min(res_ma$x_ma, na.rm = TRUE)
  ma_x_max <- max(res_ma$x_ma, na.rm = TRUE)
  ma_x_pad <- max((ma_x_max - ma_x_min) * 0.03, 0.05)
  x_lim_ma <- c(ma_x_min - ma_x_pad, ma_x_max + ma_x_pad)
  
  volcano_x_lim <- nice_limit_sym(volc_df$log2FoldChange, min_limit = 1)
  volcano_y_lim <- nice_limit_pos(volc_df$neglog10padj, min_limit = 2)
  x_lim_vol <- c(-volcano_x_lim, volcano_x_lim)
  y_lim_vol <- c(0, volcano_y_lim)
  
  ###############################################
  ## 5. MA plot
  ###############################################
  MA <- res_ma %>%
    ggplot(aes(
      x = x_ma,
      y = log2FoldChange,
      group = guide_id,
      text = paste0(
        "guide_id: ", guide_id,
        "<br>baseMean: ", signif(baseMean, 4),
        "<br>log2FC: ", signif(log2FoldChange, 4),
        "<br>padj: ", ifelse(is.na(padj), "NA", signif(padj, 4))
      )
    )) +
    ma_bg_layers(x_mid = x_mid_ma, y_lim = y_lim_ma) +
    geom_point(colour = scales::alpha(BASE_GREY, 0.50), size = 2.3) +
    geom_point(data = hits_un_ma, colour = SIG_GOLD, size = 4.0) +
    geom_point(
      data = old_vax_ma %>% filter(vax_group == "NONVAL"),
      shape = 21,
      fill = scales::alpha(VAX_FILL_GREY, 0.90),
      colour = VAX_NONVAL_BORDER,
      stroke = 0.28,
      size = 2.2
    ) +
    geom_point(
      data = old_vax_ma %>% filter(vax_group == "VAL"),
      shape = 21,
      fill = scales::alpha(VAX_FILL_GREY, 0.90),
      colour = VAX_VAL_BORDER,
      stroke = 0.28,
      size = 2.2
    ) +
    geom_point(
      data = hit_stats_ma,
      aes(color = hit_color),
      size = 2.8
    ) +
    geom_text_repel(
      data = hit_stats_ma,
      aes(label = label, color = hit_color),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.size = 0.35,
      segment.alpha = 0.9,
      direction = "y"
    ) +
    scale_color_identity() +
    theme_light(base_size = 11) +
    coord_cartesian(
      xlim = x_lim_ma,
      ylim = y_lim_ma
    ) +
    ggtitle(paste0("NONVAX MA ", CELL, " A=green, B=red")) +
    xlab("log10(baseMean + 1)") +
    ylab("log2FoldChange (A vs B)")
  ###############################################
  ## 6. Volcano plot
  ###############################################
  padj_line_y <- -log10(0.05)
  
  volcano <- volc_df %>%
    ggplot(aes(
      x = log2FoldChange,
      y = neglog10padj,
      text = paste0(
        "guide_id: ", guide_id,
        "<br>log2FC: ", signif(log2FoldChange, 4),
        "<br>padj: ", signif(padj, 4),
        "<br>-log10(padj): ", signif(neglog10padj, 4)
      )
    )) +
    volcano_bg_layers(x_lim = x_lim_vol, y_lim = y_lim_vol) +
    geom_hline(yintercept = padj_line_y, linetype = "dashed") +
    geom_point(colour = scales::alpha(BASE_GREY, 0.50), size = 0.6) +
    geom_point(data = hits_un_v, colour = SIG_GOLD, size = 2.8) +
    geom_point(
      data = old_vax_v %>% filter(vax_group == "NONVAL"),
      shape = 21,
      fill = scales::alpha(VAX_FILL_GREY, 0.90),
      colour = VAX_NONVAL_BORDER,
      stroke = 0.28,
      size = 1.8
    ) +
    geom_point(
      data = old_vax_v %>% filter(vax_group == "VAL"),
      shape = 21,
      fill = scales::alpha(VAX_FILL_GREY, 0.90),
      colour = VAX_VAL_BORDER,
      stroke = 0.28,
      size = 1.8
    ) +
    geom_point(
      data = hit_stats_v,
      aes(color = hit_color),
      size = 2.4
    ) +
    geom_text_repel(
      data = hit_stats_v,
      aes(label = label, color = hit_color),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.size = 0.35,
      segment.alpha = 0.9,
      direction = "y"
    ) +
    scale_color_identity() +
    theme_light(base_size = 11) +
    coord_cartesian(
      xlim = x_lim_vol,
      ylim = y_lim_vol
    ) +
    ggtitle(paste0("NONVAX Volcano ", CELL, " A=green, B=red")) +
    xlab("log2FoldChange (A vs B)") +
    ylab("-log10(padj)")
  
  ###############################################
  ## 7. save
  ###############################################
  outdir <- file.path(out_root, paste0("figures_nonvax_", cell))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  save_gg(MA, outdir, paste0("NONVAX_MA_", CELL, "_bgA_B"))
  save_gg(volcano, outdir, paste0("NONVAX_Volcano_", CELL, "_bgA_B"))
  
  save_plotly(MA, outdir, paste0("NONVAX_MA_", CELL, "_interactive"))
  save_plotly(volcano, outdir, paste0("NONVAX_Volcano_", CELL, "_interactive"))
  
  write_excel_csv(
    hit_stats %>%
      select(label, nv_tcr, guide_id, log2FoldChange, baseMean, pvalue, padj, hit_side),
    file.path(out_root, paste0(CELL, "_NONVAX_hits_joined.csv"))
  )
  
  write_excel_csv(
    old_vax_stats %>%
      select(tcr, guide_id, vax_group, log2FoldChange, baseMean, pvalue, padj),
    file.path(out_root, paste0(CELL, "_OLDVAX_guides_joined.csv"))
  )
  
  message("==== Done: ", CELL, " (saved to ", outdir, ") ====")
  
  invisible(list(
    results = res_tb,
    MA = MA,
    volcano = volcano,
    nonvax_hits = hit_stats,
    old_vax_stats = old_vax_stats
  ))
}

###############################################
## RUN BOTH
###############################################
out_cd4 <- run_nonvax_pipeline(
  cell = "cd4",
  count_xlsx = cd4_xlsx,
  nonvax_hits_path = nonvax_hits_path,
  sc_avail_path = sc_avail_path,
  tcr_map_path = tcr_map_path,
  out_root = base_dir
)

out_cd8 <- run_nonvax_pipeline(
  cell = "cd8",
  count_xlsx = cd8_xlsx,
  nonvax_hits_path = nonvax_hits_path,
  sc_avail_path = sc_avail_path,
  tcr_map_path = tcr_map_path,
  out_root = base_dir
)