suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(scales)
  library(ggrepel)
})

## =========================================================
## 0) Input / output
## =========================================================
bulk_file <- "C:/Users/peppa/Downloads/p108_betas_merged_typed_for_Connor.csv"
sc_file   <- "C:/Users/peppa/Downloads/sc_data_availibility.xlsx"
nv_file   <- "C:/Users/peppa/Downloads/nonvax_hits_summary.xlsx"

out_dir <- "C:/Users/peppa/Downloads"

out_png_bg   <- file.path(out_dir, "p108_bulk_background_only.png")
out_pdf_bg   <- file.path(out_dir, "p108_bulk_background_only.pdf")

out_png_hits <- file.path(out_dir, "p108_bulk_sc_nv_only.png")
out_pdf_hits <- file.path(out_dir, "p108_bulk_sc_nv_only.pdf")

out_csv_used <- file.path(out_dir, "p108_bulk_sc_nv_used_matches.csv")

LOD <- 1e-4

## =========================================================
## 1) Helper functions
## =========================================================
clean_names_simple <- function(x) {
  x %>%
    gsub("\\s+", "_", .) %>%
    gsub("\\.+", "_", .) %>%
    gsub("[^A-Za-z0-9_]", "", .)
}

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

parse_logical_loose <- function(x) {
  y <- tolower(trimws(as.character(x)))
  case_when(
    y %in% c("true", "t", "1", "yes", "y")  ~ TRUE,
    y %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

## =========================================================
## 2) Read bulk data
## =========================================================
bulk_raw <- read.csv(bulk_file, stringsAsFactors = FALSE)

df_long <- bulk_raw %>%
  select(Beta_clonotype, Type, p108_pretreatment, p108_prevax, p108_postvax, p108_w32) %>%
  mutate(Beta_clonotype = trimws(as.character(Beta_clonotype))) %>%
  pivot_longer(
    cols = c(p108_pretreatment, p108_prevax, p108_postvax, p108_w32),
    names_to = "Timepoint",
    values_to = "Frequency"
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("p108_pretreatment", "p108_prevax", "p108_postvax", "p108_w32")
    ),
    Frequency = suppressWarnings(as.numeric(Frequency)),
    Frequency = ifelse(Frequency <= 0, NA_real_, Frequency)
  )

## =========================================================
## 3) Read SC / NV Excel
## =========================================================
sc_raw <- read_excel(sc_file)
nv_raw <- read_excel(nv_file)

colnames(sc_raw) <- clean_names_simple(colnames(sc_raw))
colnames(nv_raw) <- clean_names_simple(colnames(nv_raw))

message("SC columns: ", paste(colnames(sc_raw), collapse = ", "))
message("NV columns: ", paste(colnames(nv_raw), collapse = ", "))

sc_tcr_col  <- find_col(sc_raw, c("TCR", "tcr", "TCR_ID", "TCRid"))
sc_beta_col <- find_col(sc_raw, c("beta_clonotype", "Beta_clonotype", "beta", "Beta"))
sc_val_col  <- find_col(sc_raw, c("val", "VAL", "Val"))

nv_tcr_col  <- find_col(nv_raw, c("NV_TCR", "nv_tcr", "NVTCR", "NVTCR"))
nv_beta_col <- find_col(nv_raw, c("beta_clonotype", "Beta_clonotype", "beta", "Beta"))

## =========================================================
## 4) Build highlight map
##    rules:
##    - SC val=TRUE  : red thick line + label
##    - SC val=FALSE : black thin line, no label
##    - NV           : blue line + label
## =========================================================
highlight_sc <- sc_raw %>%
  transmute(
    Beta_clonotype = trimws(as.character(.data[[sc_beta_col]])),
    Label          = paste0("TCR", as.character(.data[[sc_tcr_col]])),
    Source         = "SC",
    val            = parse_logical_loose(.data[[sc_val_col]])
  ) %>%
  filter(!is.na(Beta_clonotype), Beta_clonotype != "") %>%
  mutate(
    HighlightGroup = case_when(
      val %in% TRUE  ~ "SC_VAL_TRUE",
      val %in% FALSE ~ "SC_VAL_FALSE",
      TRUE           ~ "SC_VAL_FALSE"
    )
  ) %>%
  distinct(Beta_clonotype, Label, Source, HighlightGroup)

highlight_nv <- nv_raw %>%
  transmute(
    Beta_clonotype = trimws(as.character(.data[[nv_beta_col]])),
    Label          = paste0("NV_TCR", as.character(.data[[nv_tcr_col]])),
    Source         = "NV",
    HighlightGroup = "NV"
  ) %>%
  filter(!is.na(Beta_clonotype), Beta_clonotype != "") %>%
  distinct(Beta_clonotype, Label, Source, HighlightGroup)

highlight_map <- bind_rows(highlight_sc, highlight_nv) %>%
  distinct()

## only keep entries that actually exist in bulk data
highlight_map_used <- highlight_map %>%
  semi_join(df_long %>% distinct(Beta_clonotype), by = "Beta_clonotype")

write.csv(highlight_map_used, out_csv_used, row.names = FALSE)

message("Used SC matches: ", sum(highlight_map_used$Source == "SC"))
message("Used NV matches: ", sum(highlight_map_used$Source == "NV"))

## =========================================================
## 5) Data for plot 1: background only
##    no thick line, no labels
## =========================================================
df_bg <- df_long

## =========================================================
## 6) Data for plot 2: only SC/NV hits
##    no background lines
## =========================================================
df_hits <- df_long %>%
  inner_join(highlight_map_used, by = "Beta_clonotype")

## labels:
## - SC val=TRUE : label
## - NV          : label
## - SC val=FALSE: no label
df_label <- df_hits %>%
  filter(HighlightGroup %in% c("SC_VAL_TRUE", "NV")) %>%
  mutate(tp_idx = as.integer(Timepoint)) %>%
  filter(!is.na(Frequency)) %>%
  group_by(Beta_clonotype, Label, Source, HighlightGroup) %>%
  slice_max(tp_idx, n = 1, with_ties = FALSE) %>%
  ungroup()

## =========================================================
## 7) Shared theme pieces
## =========================================================
x_scale_shared <- scale_x_discrete(
  labels = c(
    p108_pretreatment = "Pre-treatment",
    p108_prevax       = "Post-Nivolumab",
    p108_postvax      = "Post-Vaccine",
    p108_w32          = "W32"
  ),
  expand = expansion(mult = c(0.02, 0.35))
)

y_scale_shared <- scale_y_log10(
  breaks = c(1e-3, 1e-2, 1e-1, 1),
  labels = c("0.001", "0.01", "0.1", "1")
)

theme_shared <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 80, 5.5, 5.5)
  )

## =========================================================
## 8) Plot 1: background only
##    no highlighted lines, no labels
## =========================================================
p_bg <- ggplot() +
  geom_line(
    data = df_bg %>% filter(Type == "Existing"),
    aes(Timepoint, Frequency, group = Beta_clonotype),
    color = "grey80", linewidth = 0.35, alpha = 0.7, na.rm = TRUE
  ) +
  geom_point(
    data = df_bg %>% filter(Type == "Existing"),
    aes(Timepoint, Frequency, group = Beta_clonotype),
    color = "grey80", size = 0.7, alpha = 0.7, na.rm = TRUE
  ) +
  geom_line(
    data = df_bg %>% filter(Type == "Post-Nivolumab"),
    aes(Timepoint, Frequency, group = Beta_clonotype),
    color = "#D9A441", linewidth = 0.35, alpha = 0.75, na.rm = TRUE
  ) +
  geom_point(
    data = df_bg %>% filter(Type == "Post-Nivolumab"),
    aes(Timepoint, Frequency, group = Beta_clonotype),
    color = "#D9A441", size = 0.7, alpha = 0.75, na.rm = TRUE
  ) +
  geom_line(
    data = df_bg %>% filter(Type == "Post-Vaccine"),
    aes(Timepoint, Frequency, group = Beta_clonotype),
    color = "#2C7FB8", linewidth = 0.35, alpha = 0.8, na.rm = TRUE
  ) +
  geom_point(
    data = df_bg %>% filter(Type == "Post-Vaccine"),
    aes(Timepoint, Frequency, group = Beta_clonotype),
    color = "#2C7FB8", size = 0.7, alpha = 0.8, na.rm = TRUE
  ) +
  geom_hline(
    yintercept = LOD,
    color = "red",
    linewidth = 0.7,
    linetype = "longdash"
  ) +
  annotate(
    "text",
    x = 1,
    y = LOD,
    label = "LOD",
    color = "red",
    vjust = -0.6,
    hjust = 0,
    size = 3.8
  ) +
  x_scale_shared +
  y_scale_shared +
  coord_cartesian(ylim = c(1e-4, 1.5), clip = "off") +
  labs(x = NULL, y = "Frequency") +
  theme_shared

## =========================================================
## 9) Plot 2: only SC/NV highlighted hits
##    no background lines
##    SC val=FALSE: black thin line, no label
## =========================================================
p_hits <- ggplot() +
  ## SC val = TRUE : red thick
  geom_line(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_TRUE"),
    aes(Timepoint, Frequency, group = interaction(Beta_clonotype, Label)),
    color = "red", linewidth = 1.6, alpha = 0.98, na.rm = TRUE
  ) +
  geom_point(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_TRUE"),
    aes(Timepoint, Frequency, group = interaction(Beta_clonotype, Label)),
    color = "red", size = 1.4, alpha = 0.98, na.rm = TRUE
  ) +
  
  ## SC val = FALSE : black thin, no label
  geom_line(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_FALSE"),
    aes(Timepoint, Frequency, group = interaction(Beta_clonotype, Label)),
    color = "black", linewidth = 0.55, alpha = 0.95, na.rm = TRUE
  ) +
  geom_point(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_FALSE"),
    aes(Timepoint, Frequency, group = interaction(Beta_clonotype, Label)),
    color = "black", size = 0.8, alpha = 0.95, na.rm = TRUE
  ) +
  
  ## NV : blue thick
  geom_line(
    data = df_hits %>% filter(HighlightGroup == "NV"),
    aes(Timepoint, Frequency, group = interaction(Beta_clonotype, Label)),
    color = "blue", linewidth = 1.6, alpha = 0.98, na.rm = TRUE
  ) +
  geom_point(
    data = df_hits %>% filter(HighlightGroup == "NV"),
    aes(Timepoint, Frequency, group = interaction(Beta_clonotype, Label)),
    color = "blue", size = 1.4, alpha = 0.98, na.rm = TRUE
  ) +
  
  ## labels only for SC val TRUE and NV
  geom_text_repel(
    data = df_label %>% filter(HighlightGroup == "SC_VAL_TRUE"),
    aes(Timepoint, Frequency, label = Label),
    color = "red",
    size = 3.3,
    direction = "y",
    hjust = 0,
    nudge_x = 0.20,
    segment.color = "red",
    segment.size = 0.35,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    seed = 123
  ) +
  geom_text_repel(
    data = df_label %>% filter(HighlightGroup == "NV"),
    aes(Timepoint, Frequency, label = Label),
    color = "blue",
    size = 3.3,
    direction = "y",
    hjust = 0,
    nudge_x = 0.20,
    segment.color = "blue",
    segment.size = 0.35,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    seed = 123
  ) +
  geom_hline(
    yintercept = LOD,
    color = "red",
    linewidth = 0.7,
    linetype = "longdash"
  ) +
  annotate(
    "text",
    x = 1,
    y = LOD,
    label = "LOD",
    color = "red",
    vjust = -0.6,
    hjust = 0,
    size = 3.8
  ) +
  x_scale_shared +
  y_scale_shared +
  coord_cartesian(ylim = c(1e-4, 1.5), clip = "off") +
  labs(x = NULL, y = "Frequency") +
  theme_shared

## =========================================================
## 10) Print and save
## =========================================================


ggsave(out_png_bg, plot = p_bg, width = 8.5, height = 5.8, units = "in", dpi = 600)
ggsave(out_pdf_bg, plot = p_bg, width = 8.5, height = 5.8, units = "in", dpi = 600, device = cairo_pdf)

ggsave(out_png_hits, plot = p_hits, width = 8.5, height = 5.8, units = "in", dpi = 600)
ggsave(out_pdf_hits, plot = p_hits, width = 8.5, height = 5.8, units = "in", dpi = 600, device = cairo_pdf)
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(scales)
})

## =========================================================
## 0) Input / output
## =========================================================
bulk_file <- "C:/Users/peppa/Downloads/p108_betas_merged_typed_for_Connor.csv"
sc_file   <- "C:/Users/peppa/Downloads/sc_data_availibility.xlsx"
nv_file   <- "C:/Users/peppa/Downloads/nonvax_hits_summary.xlsx"

out_dir <- "C:/Users/peppa/Downloads"

out_png_hits <- file.path(out_dir, "p108_bulk_sc_nv_only_no_background_spacedlabels_W32sorted.png")
out_pdf_hits <- file.path(out_dir, "p108_bulk_sc_nv_only_no_background_spacedlabels_W32sorted.pdf")
out_csv_used <- file.path(out_dir, "p108_bulk_sc_nv_used_matches.csv")

LOD  <- 1e-4
YMIN <- 1e-4
YMAX <- 1.5

## =========================================================
## 1) Helper functions
## =========================================================
clean_names_simple <- function(x) {
  x %>%
    gsub("\\s+", "_", .) %>%
    gsub("\\.+", "_", .) %>%
    gsub("[^A-Za-z0-9_]", "", .)
}

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

parse_logical_loose <- function(x) {
  y <- tolower(trimws(as.character(x)))
  case_when(
    y %in% c("true", "t", "1", "yes", "y")  ~ TRUE,
    y %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

## 在 log10 尺度上生成均匀分开的标签位置
spread_y_log <- function(n, ymin = 1e-4, ymax = 1.5, top_pad = 0.10, bottom_pad = 0.14) {
  if (n == 0) return(numeric(0))
  if (n == 1) return(10^((log10(ymin) + log10(ymax)) / 2))
  
  log_min <- log10(ymin)
  log_max <- log10(ymax)
  rng <- log_max - log_min
  
  usable_min <- log_min + rng * bottom_pad
  usable_max <- log_max - rng * top_pad
  
  10^seq(usable_max, usable_min, length.out = n)
}

## =========================================================
## 2) Read bulk data
## =========================================================
bulk_raw <- read.csv(bulk_file, stringsAsFactors = FALSE)

df_long <- bulk_raw %>%
  select(Beta_clonotype, Type, p108_pretreatment, p108_prevax, p108_postvax, p108_w32) %>%
  mutate(Beta_clonotype = trimws(as.character(Beta_clonotype))) %>%
  pivot_longer(
    cols = c(p108_pretreatment, p108_prevax, p108_postvax, p108_w32),
    names_to = "Timepoint",
    values_to = "Frequency"
  ) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("p108_pretreatment", "p108_prevax", "p108_postvax", "p108_w32")
    ),
    Timepoint_num = case_when(
      Timepoint == "p108_pretreatment" ~ 1,
      Timepoint == "p108_prevax"       ~ 2,
      Timepoint == "p108_postvax"      ~ 3,
      Timepoint == "p108_w32"          ~ 4
    ),
    Frequency = suppressWarnings(as.numeric(Frequency)),
    Frequency = ifelse(Frequency <= 0, NA_real_, Frequency)
  )
df_w32 <- df_long %>%
  filter(Timepoint == "p108_w32") %>%
  select(Beta_clonotype, W32_Frequency = Frequency)

## =========================================================
## 3) Read SC / NV Excel
## =========================================================
sc_raw <- read_excel(sc_file)
nv_raw <- read_excel(nv_file)

colnames(sc_raw) <- clean_names_simple(colnames(sc_raw))
colnames(nv_raw) <- clean_names_simple(colnames(nv_raw))

message("SC columns: ", paste(colnames(sc_raw), collapse = ", "))
message("NV columns: ", paste(colnames(nv_raw), collapse = ", "))

sc_tcr_col  <- find_col(sc_raw, c("TCR", "tcr", "TCR_ID", "TCRid"))
sc_beta_col <- find_col(sc_raw, c("beta_clonotype", "Beta_clonotype", "beta", "Beta"))
sc_val_col  <- find_col(sc_raw, c("val", "VAL", "Val"))

nv_tcr_col  <- find_col(nv_raw, c("NV_TCR", "nv_tcr", "NVTCR", "NVTCR"))
nv_beta_col <- find_col(nv_raw, c("beta_clonotype", "Beta_clonotype", "beta", "Beta"))



## =========================================================
## 4) Build highlight map
## =========================================================
highlight_sc <- sc_raw %>%
  transmute(
    Beta_clonotype = trimws(as.character(.data[[sc_beta_col]])),
    Label          = paste0("TCR", as.character(.data[[sc_tcr_col]])),
    Source         = "SC",
    val            = parse_logical_loose(.data[[sc_val_col]])
  ) %>%
  filter(!is.na(Beta_clonotype), Beta_clonotype != "") %>%
  mutate(
    HighlightGroup = case_when(
      val %in% TRUE  ~ "SC_VAL_TRUE",
      val %in% FALSE ~ "SC_VAL_FALSE",
      TRUE           ~ "SC_VAL_FALSE"
    )
  ) %>%
  distinct(Beta_clonotype, Label, Source, HighlightGroup)

highlight_nv <- nv_raw %>%
  transmute(
    Beta_clonotype = trimws(as.character(.data[[nv_beta_col]])),
    Label          = paste0("NV_TCR", as.character(.data[[nv_tcr_col]])),
    Source         = "NV",
    HighlightGroup = "NV"
  ) %>%
  filter(!is.na(Beta_clonotype), Beta_clonotype != "") %>%
  distinct(Beta_clonotype, Label, Source, HighlightGroup)

highlight_map <- bind_rows(highlight_sc, highlight_nv) %>%
  distinct()

highlight_map_used <- highlight_map %>%
  semi_join(df_long %>% distinct(Beta_clonotype), by = "Beta_clonotype")

write.csv(highlight_map_used, out_csv_used, row.names = FALSE)

message("Used SC matches: ", sum(highlight_map_used$Source == "SC"))
message("Used NV matches: ", sum(highlight_map_used$Source == "NV"))

## =========================================================
## 5) Data for hit-only plot
## =========================================================
df_hits <- df_long %>%
  inner_join(highlight_map_used, by = "Beta_clonotype")

## =========================================================
## 6) Build label table manually

## =========================================================
df_label <- df_hits %>%
  filter(HighlightGroup %in% c("SC_VAL_TRUE", "NV")) %>%
  filter(!is.na(Frequency)) %>%
  group_by(Beta_clonotype, Label, Source, HighlightGroup) %>%
  slice_max(Timepoint_num, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(df_w32, by = "Beta_clonotype") %>%
  mutate(
    W32_Frequency = ifelse(is.na(W32_Frequency), -Inf, W32_Frequency)
  ) %>%
  arrange(desc(W32_Frequency), desc(Frequency), Label) %>%
  mutate(
    x_point_num = Timepoint_num,
    x_label_num = 4.58,
    label_y     = spread_y_log(n(), ymin = YMIN, ymax = YMAX, top_pad = 0.10, bottom_pad = 0.14),
    text_color  = case_when(
      HighlightGroup == "SC_VAL_TRUE" ~ "red",
      HighlightGroup == "NV"          ~ "blue",
      TRUE                            ~ "black"
    )
  )

df_label_red  <- df_label %>% filter(HighlightGroup == "SC_VAL_TRUE")
df_label_blue <- df_label %>% filter(HighlightGroup == "NV")

## =========================================================
## 7) Plot
## =========================================================
p_hits <- ggplot() +
  ## SC val = TRUE : red thick
  geom_line(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_TRUE"),
    aes(x = Timepoint_num, y = Frequency, group = interaction(Beta_clonotype, Label)),
    color = "red", linewidth = 1.6, alpha = 0.98, na.rm = TRUE
  ) +
  geom_point(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_TRUE"),
    aes(x = Timepoint_num, y = Frequency, group = interaction(Beta_clonotype, Label)),
    color = "red", size = 1.4, alpha = 0.98, na.rm = TRUE
  ) +
  
  ## SC val = FALSE : black thin, no label
  geom_line(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_FALSE"),
    aes(x = Timepoint_num, y = Frequency, group = interaction(Beta_clonotype, Label)),
    color = "black", linewidth = 0.55, alpha = 0.95, na.rm = TRUE
  ) +
  geom_point(
    data = df_hits %>% filter(HighlightGroup == "SC_VAL_FALSE"),
    aes(x = Timepoint_num, y = Frequency, group = interaction(Beta_clonotype, Label)),
    color = "black", size = 0.8, alpha = 0.95, na.rm = TRUE
  ) +
  
  ## NV : blue thick
  geom_line(
    data = df_hits %>% filter(HighlightGroup == "NV"),
    aes(x = Timepoint_num, y = Frequency, group = interaction(Beta_clonotype, Label)),
    color = "blue", linewidth = 1.6, alpha = 0.98, na.rm = TRUE
  ) +
  geom_point(
    data = df_hits %>% filter(HighlightGroup == "NV"),
    aes(x = Timepoint_num, y = Frequency, group = interaction(Beta_clonotype, Label)),
    color = "blue", size = 1.4, alpha = 0.98, na.rm = TRUE
  ) +
  
  ## connector lines
  geom_segment(
    data = df_label_red,
    aes(x = x_point_num, xend = x_label_num - 0.03, y = Frequency, yend = label_y),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 0.40,
    alpha = 0.9
  ) +
  geom_segment(
    data = df_label_blue,
    aes(x = x_point_num, xend = x_label_num - 0.03, y = Frequency, yend = label_y),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.40,
    alpha = 0.9
  ) +
  
  ## label text
  geom_text(
    data = df_label_red,
    aes(x = x_label_num, y = label_y, label = Label),
    inherit.aes = FALSE,
    color = "red",
    hjust = 0,
    size = 3.4
  ) +
  geom_text(
    data = df_label_blue,
    aes(x = x_label_num, y = label_y, label = Label),
    inherit.aes = FALSE,
    color = "blue",
    hjust = 0,
    size = 3.4
  ) +
  
  geom_hline(
    yintercept = LOD,
    color = "red",
    linewidth = 0.7,
    linetype = "longdash"
  ) +
  annotate(
    "text",
    x = 1,
    y = LOD,
    label = "LOD",
    color = "red",
    vjust = -0.6,
    hjust = 0,
    size = 3.8
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("Pre-treatment", "Post-Nivolumab", "Post-Vaccine", "W32"),
    limits = c(1, 4.95),
    expand = c(0.02, 0.02)
  ) +
  scale_y_log10(
    breaks = c(1e-3, 1e-2, 1e-1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  coord_cartesian(ylim = c(YMIN, YMAX), clip = "off") +
  labs(x = NULL, y = "Frequency") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 120, 5.5, 5.5)
  )

## =========================================================
## 8) Print and save
## =========================================================
print(p_hits)

ggsave(out_png_hits, plot = p_hits, width = 9.0, height = 6.6, units = "in", dpi = 600)
ggsave(out_pdf_hits, plot = p_hits, width = 9.0, height = 6.6, units = "in", dpi = 600, device = cairo_pdf)