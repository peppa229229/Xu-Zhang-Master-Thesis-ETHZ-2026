###############################################
## FULL PIPELINE (CD4 + CD8)
## UNIFORM EXPORT SIZE FOR ALL FIGURES
## + FIXED TCR18 SPECIAL RULE BUG
## + INTERACTIVE ggplot (hover shows guide_id)
###############################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(readxl)
library(scales)

library(ragg)   # agg_png
library(grid)   # pheatmap drawing

# interactive
library(plotly)
library(htmlwidgets)

setwd("C:/Users/peppa/Downloads/run3/")

###############################################
## GLOBAL FIGURE SIZE (ALL FIGURES SAME)
###############################################
FIG_W_MM <- 190
FIG_H_MM <- 140
FIG_DPI  <- 600

# for pheatmap (keep consistent physical size)
MM_PER_INCH <- 25.4
PH_DPI <- 300  # pheatmap PNG dpi

###############################################
## Helper functions
###############################################

clean_names_simple <- function(nms) {
  nms <- gsub("[^A-Za-z0-9]+", "_", nms)
  nms <- gsub("^_|_$", "", nms)
  tolower(nms)
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

## --- colors ---
VAL_GREEN <- "#2E7D32"  # A-side
VAL_RED   <- "#C62828"  # B-side

ma_bg_layers <- function(x_mid, y_lim = c(-2, 2),
                         fill_pos = "#DFF0D8", fill_neg = "#F8D7DA",
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
                              fill_pos = "#DFF0D8", fill_neg = "#F8D7DA",
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

## ---- unified ggplot saver (ALL SAME SIZE) ----
save_gg <- function(p, outdir, fname, width_mm = FIG_W_MM, height_mm = FIG_H_MM, dpi = FIG_DPI) {
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

## ---- unified pheatmap saver (same physical size as ggplots) ----
save_pheatmap_uniform <- function(ph_obj, outdir, fname_base,
                                  fig_w_mm = FIG_W_MM, fig_h_mm = FIG_H_MM,
                                  png_dpi = PH_DPI) {
  pdf_w_in <- fig_w_mm / MM_PER_INCH
  pdf_h_in <- fig_h_mm / MM_PER_INCH
  
  png_w_px <- round((fig_w_mm / MM_PER_INCH) * png_dpi)
  png_h_px <- round((fig_h_mm / MM_PER_INCH) * png_dpi)
  
  ragg::agg_png(
    filename = file.path(outdir, paste0(fname_base, ".png")),
    width = png_w_px, height = png_h_px, res = png_dpi, background = "white"
  )
  grid::grid.newpage(); grid::grid.draw(ph_obj$gtable)
  dev.off()
  
  pdf(file.path(outdir, paste0(fname_base, ".pdf")), width = pdf_w_in, height = pdf_h_in)
  grid::grid.newpage(); grid::grid.draw(ph_obj$gtable)
  dev.off()
}

## ---- interactive html saver ----
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
## Main pipeline
###############################################

run_pipeline <- function(cell = c("cd4", "cd8"),
                         count_xlsx_dir = "C:/Users/peppa/Downloads/run3/",
                         sc_avail_path  = "C:/Users/peppa/Downloads/sc_data_availibility.xlsx",
                         tcr_map_path   = "C:/Users/peppa/Downloads/tcr_to_guide_id_matches.xlsx") {
  
  cell <- match.arg(cell)
  CELL <- toupper(cell)
  
  message("==== Running: ", CELL, " ====")
  
  ###############################################
  ## 1. Prepare count matrix and sample metadata
  ###############################################
  
  counts_file <- file.path(count_xlsx_dir, paste0(cell, ".xlsx"))
  COUNTS <- read_excel(counts_file)
  
  count_data <- COUNTS %>%
    as.data.frame() %>%
    column_to_rownames("guide_id")
  
  count_data <- round(as.matrix(count_data))
  
  condition <- ifelse(grepl("_A_", colnames(count_data)), "B", "A") # inverted
  
  meta <- data.frame(
    row.names = colnames(count_data),
    condition = factor(condition)
  )
  
  stopifnot(all(colnames(count_data) %in% rownames(meta)))
  stopifnot(all(colnames(count_data) == rownames(meta)))
  
  ###############################################
  ## 2. estimateSizeFactors first, export normalized counts
  ###############################################
  
  dds_sf <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = meta,
    design    = ~ condition
  )
  
  dds_sf <- estimateSizeFactors(dds_sf)
  
  normalized_counts <- counts(dds_sf, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("guide_id")
  
  ###############################################
  ## 3. QC: PCA, sample correlation, boxplots
  ###############################################
  
  rld <- rlog(dds_sf, blind = TRUE)
  
  pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p_pca <- ggplot(
    pcaData,
    aes(x = PC1, y = PC2, color = condition, label = name)
  ) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_light(base_size = 11) +
    ggtitle(paste0("PCA ", CELL, " A vs B"))
  
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  counts_unnormalized_long <- count_data %>%
    as.data.frame() %>%
    rownames_to_column("guide_id") %>%
    pivot_longer(cols = -guide_id, names_to = "sample", values_to = "counts")
  
  normalized_counts_long <- normalized_counts %>%
    pivot_longer(cols = -guide_id, names_to = "sample", values_to = "norm_counts")
  
  cul <- counts_unnormalized_long %>%
    ggplot(aes(x = sample, y = log10(counts + 1))) +
    geom_boxplot() +
    theme_light(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90),
      plot.title  = element_text(size = 8, face = "bold")
    ) +
    xlab("") +
    ylab("Counts (log10)") +
    ggtitle(paste0("Un-normalized data ", CELL, " A vs B"))
  
  ncl <- normalized_counts_long %>%
    ggplot(aes(x = sample, y = log10(norm_counts + 1))) +
    geom_boxplot() +
    theme_light(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90),
      plot.title  = element_text(size = 8, face = "bold")
    ) +
    xlab("") +
    ylab("Counts (log10)") +
    ggtitle(paste0("Normalized data ", CELL, " A vs B"))
  
  p_hist <- normalized_counts_long %>%
    ggplot(aes(x = log10(norm_counts + 1))) +
    geom_histogram(aes(y = ..density..), binwidth = 0.25,
                   colour = "black", fill = "white") +
    geom_density(alpha = 0.2, fill = "#FF6666") +
    facet_wrap(sample ~ .) +
    theme_light(base_size = 11) +
    ylab("Density") +
    xlab("Normalized reads per guide (log10)") +
    ggtitle(paste0("Normalized counts distribution: ", CELL))
  
  ###############################################
  ## 4. Run DESeq (A vs B)
  ###############################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = meta,
    design    = ~ condition
  )
  
  dds <- DESeq(dds)
  
  contrast_condition <- c("condition", "A", "B")
  res_table_unshrunken <- results(
    dds,
    contrast = contrast_condition,
    alpha = 0.05
  )
  
  res_table_unshrunken_tb <- res_table_unshrunken %>%
    data.frame() %>%
    rownames_to_column(var = "guide_id") %>%
    as_tibble()
  
  hits_un <- res_table_unshrunken_tb %>%
    filter(!is.na(padj), padj < 0.05)
  
  ###############################################
  ## 5. sc_data_availibility + TCR->guide_id mapping
  ###############################################
  
  sc_avail <- read_excel(sc_avail_path) %>%
    rename_with(clean_names_simple)
  
  tcr_map <- read_excel(tcr_map_path) %>%
    rename_with(clean_names_simple)
  
  sc_cell_col <- if (cell %in% names(sc_avail)) cell else pick_first(names(sc_avail), c(paste0("^", cell, "$"), cell))
  sc_val_col  <- if ("val" %in% names(sc_avail)) "val" else pick_first(names(sc_avail), c("^val$", "val"))
  sc_tcr_col  <- if ("tcr" %in% names(sc_avail)) "tcr" else pick_first(names(sc_avail), c("^tcr$", "tcr"))
  
  map_tcr_col   <- if ("tcr" %in% names(tcr_map)) "tcr" else pick_first(names(tcr_map), c("^tcr$", "tcr"))
  map_guide_col <- if ("guide_id" %in% names(tcr_map)) "guide_id" else pick_first(names(tcr_map), c("^guide_id$", "guide", "guideid"))
  
  if (any(is.na(c(sc_cell_col, sc_val_col, sc_tcr_col)))) {
    stop("sc_data_availibility.xlsx 缺少必要列（", CELL, " / VAL / TCR）。请检查表头。")
  }
  if (any(is.na(c(map_tcr_col, map_guide_col)))) {
    stop("tcr_to_guide_id_matches.xlsx 缺少必要列（TCR / guide_id）。请检查表头。")
  }
  
  sc_avail <- sc_avail %>%
    mutate(
      cell_flag = to_logical(.data[[sc_cell_col]]),
      val       = to_logical(.data[[sc_val_col]]),
      tcr       = as.character(.data[[sc_tcr_col]])
    )
  
  ## ---- SPECIAL RULE (FIXED): TCR18 VAL in CD4 only ----
  if (cell == "cd4") {
    sc_avail <- sc_avail %>% mutate(val = ifelse(tcr == "18", TRUE, val))
  } else if (cell == "cd8") {
    sc_avail <- sc_avail %>% mutate(val = ifelse(tcr == "18", FALSE, val))
  }
  
  tcr_map <- tcr_map %>%
    transmute(
      tcr      = as.character(.data[[map_tcr_col]]),
      guide_id = as.character(.data[[map_guide_col]])
    ) %>%
    filter(!is.na(tcr), !is.na(guide_id)) %>%
    distinct(tcr, guide_id)
  
  # CELL==TRUE -> guide_id list
  cell_guides <- sc_avail %>%
    filter(cell_flag) %>%
    distinct(tcr) %>%
    left_join(tcr_map, by = "tcr") %>%
    filter(!is.na(guide_id)) %>%
    distinct(guide_id)
  
  cell_guides_stats <- cell_guides %>%
    inner_join(res_table_unshrunken_tb, by = "guide_id")
  
  # CELL==TRUE & VAL==TRUE -> guide_id + label "TCR{编号}"
  val_guides <- sc_avail %>%
    filter(cell_flag, val) %>%
    distinct(tcr) %>%
    left_join(tcr_map, by = "tcr") %>%
    filter(!is.na(guide_id)) %>%
    distinct(tcr, guide_id) %>%
    mutate(label = paste0("TCR", tcr))
  
  val_guides_stats <- val_guides %>%
    inner_join(res_table_unshrunken_tb, by = "guide_id") %>%
    distinct(guide_id, tcr, .keep_all = TRUE) %>%
    mutate(
      val_side  = ifelse(log2FoldChange > 0, "A", "B"),
      val_color = ifelse(log2FoldChange > 0, VAL_GREEN, VAL_RED)
    )
  
  # CELL-only (exclude VAL TRUE): black + X
  cell_only_stats <- cell_guides_stats %>%
    anti_join(val_guides_stats %>% distinct(guide_id), by = "guide_id")
  
  # ---- Filter only for plotting (avoid NA / Inf) ----
  res_ma <- res_table_unshrunken_tb %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  hits_un_ma <- hits_un %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  cell_guides_stats_ma <- cell_guides_stats %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  cell_only_stats_ma <- cell_only_stats %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  val_guides_stats_ma <- val_guides_stats %>%
    filter(!is.na(baseMean), is.finite(baseMean), baseMean >= 0,
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      val_side  = ifelse(log2FoldChange > 0, "A", "B"),
      val_color = ifelse(log2FoldChange > 0, VAL_GREEN, VAL_RED)
    )
  
  ###############################################
  ## 6. MA plots (INTERACTIVE: hover shows guide_id)
  ###############################################
  
  y_lim_ma <- c(-2, 2)
  x_mid_ma <- median(log10(res_ma$baseMean + 1), na.rm = TRUE)
  
  MA_1 <- res_ma %>%
    ggplot(aes(
      x = log10(baseMean + 1),
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
    geom_point(colour = scales::alpha("grey80", 0.5), size = 2.3) +
    geom_point(data = hits_un_ma, colour = scales::alpha("gold", 0.35), size = 4) +
    geom_point(data = cell_guides_stats_ma, colour = "red", size = 2.6) +
    theme_light(base_size = 11) +
    coord_cartesian(ylim = y_lim_ma) +
    ggtitle(paste0("MA_1 ", CELL, " guides in red"))
  
  MA_2 <- res_ma %>%
    ggplot(aes(
      x = log10(baseMean + 1),
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
    geom_point(colour = scales::alpha("grey80", 0.5), size = 2.3) +
    geom_point(data = hits_un_ma, colour = scales::alpha("gold", 0.35), size = 4) +
    geom_point(data = cell_only_stats_ma, colour = "black", size = 2.2) +
    geom_point(data = cell_only_stats_ma, shape = 4, colour = "black", size = 2.2, stroke = 0.9) +
    geom_point(
      data = val_guides_stats_ma,
      aes(color = val_color),
      size = 2.6
    ) +
    geom_text_repel(
      data = val_guides_stats_ma,
      aes(label = label, color = val_color),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.color = val_guides_stats_ma$val_color,
      segment.size = 0.35,
      segment.alpha = 0.9
    ) +
    scale_color_identity() +
    theme_light(base_size = 11) +
    coord_cartesian(ylim = y_lim_ma) +
    ggtitle(paste0("MA_2 ", CELL, " VAL: A=green, B=red"))
  
  ###############################################
  ## 7. Volcano plots (INTERACTIVE: hover shows guide_id)
  ###############################################
  
  volcano_base <- res_table_unshrunken_tb %>%
    filter(!is.na(padj), is.finite(padj),
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  padj_line_y <- -log10(0.05)
  
  hits_un_v <- hits_un %>%
    filter(!is.na(padj), is.finite(padj),
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  cell_guides_stats_v <- cell_guides_stats %>%
    filter(!is.na(padj), is.finite(padj),
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  cell_only_stats_v <- cell_only_stats %>%
    filter(!is.na(padj), is.finite(padj),
           !is.na(log2FoldChange), is.finite(log2FoldChange))
  
  val_guides_stats_v <- val_guides_stats %>%
    filter(!is.na(padj), is.finite(padj),
           !is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
    mutate(
      val_side  = ifelse(log2FoldChange > 0, "A", "B"),
      val_color = ifelse(log2FoldChange > 0, VAL_GREEN, VAL_RED)
    )
  
  x_lim_vol <- c(-3, 3)
  y_lim_vol <- if (cell == "cd4") c(0, 15) else c(0, 60)
  
  volcano_1 <- volcano_base %>%
    ggplot(aes(
      x = log2FoldChange,
      y = -log10(padj),
      text = paste0(
        "guide_id: ", guide_id,
        "<br>log2FC: ", signif(log2FoldChange, 4),
        "<br>padj: ", signif(padj, 4),
        "<br>-log10(padj): ", signif(-log10(padj), 4)
      )
    )) +
    volcano_bg_layers(x_lim = x_lim_vol, y_lim = y_lim_vol) +
    geom_hline(yintercept = padj_line_y, linetype = "dashed") +
    geom_point(colour = scales::alpha("grey80", 0.5), size = 0.6) +
    geom_point(data = hits_un_v, colour = scales::alpha("gold", 0.35), size = 2.8) +
    geom_point(data = cell_guides_stats_v, colour = "red", size = 2.0) +
    theme_light(base_size = 11) +
    xlim(x_lim_vol[1], x_lim_vol[2]) +
    coord_cartesian(ylim = y_lim_vol) +
    ggtitle(paste0("Volcano_1 ", CELL, " guides in red"))
  
  volcano_2 <- volcano_base %>%
    ggplot(aes(
      x = log2FoldChange,
      y = -log10(padj),
      text = paste0(
        "guide_id: ", guide_id,
        "<br>log2FC: ", signif(log2FoldChange, 4),
        "<br>padj: ", signif(padj, 4),
        "<br>-log10(padj): ", signif(-log10(padj), 4)
      )
    )) +
    volcano_bg_layers(x_lim = x_lim_vol, y_lim = y_lim_vol) +
    geom_hline(yintercept = padj_line_y, linetype = "dashed") +
    geom_point(colour = scales::alpha("grey80", 0.5), size = 0.6) +
    geom_point(data = hits_un_v, colour = scales::alpha("gold", 0.35), size = 2.8) +
    geom_point(data = cell_only_stats_v, colour = "black", size = 1.8) +
    geom_point(data = cell_only_stats_v, shape = 4, colour = "black", size = 2.0, stroke = 0.9) +
    geom_point(
      data = val_guides_stats_v,
      aes(color = val_color),
      size = 2.2
    ) +
    geom_text_repel(
      data = val_guides_stats_v,
      aes(label = label, color = val_color),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.15,
      min.segment.length = 0,
      segment.color = val_guides_stats_v$val_color,
      segment.size = 0.35,
      segment.alpha = 0.9
    ) +
    scale_color_identity() +
    theme_light(base_size = 11) +
    xlim(x_lim_vol[1], x_lim_vol[2]) +
    coord_cartesian(ylim = y_lim_vol) +
    ggtitle(paste0("Volcano_2 ", CELL, " VAL: A=green, B=red"))
  
  ###############################################
  ## 8. Export results
  ###############################################
  
  res_table_annotated <- res_table_unshrunken_tb %>%
    mutate(
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
    res_table_annotated,
    paste0(CELL, "_DESeq2_results_A_vs_B_annotated.csv")
  )
  
  ###############################################
  ## 9. Save figures (ALL SAME SIZE) + interactive HTML
  ###############################################
  
  outdir <- paste0("figures_", cell)
  dir.create(outdir, showWarnings = FALSE)
  
  save_gg(p_pca, outdir, paste0("PCA_", CELL, "_A_vs_B"))
  save_gg(cul, outdir, "Boxplot_raw")
  save_gg(ncl, outdir, "Boxplot_normalized")
  save_gg(cul + ncl, outdir, "Boxplot_raw_vs_normalized")
  save_gg(p_hist, outdir, "NormCount_hist_density_by_sample")
  
  save_gg(MA_1, outdir, paste0("MA_1_", CELL, "_bgA_B"))
  save_gg(MA_2, outdir, paste0("MA_2_", CELL, "_bgA_B_VAL_signColor"))
  save_gg(volcano_1, outdir, paste0("Volcano_1_", CELL, "_bgA_B"))
  save_gg(volcano_2, outdir, paste0("Volcano_2_", CELL, "_bgA_B_VAL_signColor"))
  
  # interactive html (hover shows guide_id)
  save_plotly(MA_1, outdir, paste0("MA_1_", CELL, "_interactive"))
  save_plotly(MA_2, outdir, paste0("MA_2_", CELL, "_interactive"))
  save_plotly(volcano_1, outdir, paste0("Volcano_1_", CELL, "_interactive"))
  save_plotly(volcano_2, outdir, paste0("Volcano_2_", CELL, "_interactive"))
  
  ph_cor <- pheatmap(rld_cor, silent = TRUE)
  save_pheatmap_uniform(ph_cor, outdir, "Sample_correlation_heatmap")
  
  message("==== Done: ", CELL, " (saved to ", outdir, ") ====")
  
  invisible(list(
    p_pca = p_pca,
    MA_1 = MA_1, MA_2 = MA_2,
    volcano_1 = volcano_1, volcano_2 = volcano_2
  ))
}

###############################################
## Run BOTH CD4 and CD8
###############################################
out_cd4 <- run_pipeline("cd4")
out_cd8 <- run_pipeline("cd8")