# ========================================================
# Arm-level CNV analysis and heatmap plotting
# ========================================================

# ----------------------------------------
# Set working directory and load libraries
# ----------------------------------------
data_dir <- "/home//mjehangir/telomere-sv-analysis/data/"
setwd(data_dir)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(patchwork)
library(tibble)
library(grid)

# ----------------------------------------
# Read CNV segmentation file and chromosome arms info
#    Required columns: ID, chrom, loc.start, loc.end, seg.mean
# ----------------------------------------
df_cnv <- read.delim("Mean_normalized_v3_tab.seg", header = TRUE)

df_arms <- read.delim("chm13_p_q_arm_length.tsv", header = TRUE)

# Add "chr" prefix to chromosome names
df_cnv$chrom <- paste0("chr", df_cnv$chrom)

# ----------------------------------------
# Assign Gain, Loss, Neutral based on seg.mean
# ----------------------------------------
df_cnv$Gain_Loss <- case_when(
  df_cnv$seg.mean < -0.3 ~ "Loss",   # Negative CNV -> Loss
  df_cnv$seg.mean > 0.3  ~ "Gain",   # Positive CNV -> Gain
  TRUE ~ "Neutral"                    # Everything else -> Neutral
)

# ----------------------------------------
# Merge CNVs with arms to assign p/q and arm length
# ----------------------------------------
df_cnv_with_arms <- df_cnv %>%
  rowwise() %>%
  mutate(
    Arm = df_arms %>%
      filter(Chromosome == chrom,
             loc.start >= Start,
             loc.end <= End) %>%
      pull(Arm) %>% first(),

    Arm_Length = df_arms %>%
      filter(Chromosome == chrom,
             loc.start >= Start,
             loc.end <= End) %>%
      pull(Length) %>% first()
  ) %>%
  ungroup()

# ----------------------------------------
# Summarize CNVs per arm
# ----------------------------------------
segment_summary_by_arm <- df_cnv_with_arms %>%
  filter(Gain_Loss %in% c("Gain", "Loss")) %>%
  mutate(cnv_length = loc.end - loc.start + 1) %>%
  group_by(ID, chrom, Arm, Gain_Loss, Arm_Length) %>%
  summarise(
    total_cnv_length = sum(cnv_length),
    segment_count = n(),
    .groups = 'drop'
  )

write.table(segment_summary_by_arm, "arm_level_CNV_segments.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------
# Calculate CNV rates and proportions per Mb
# ----------------------------------------
cnv_summary_final <- segment_summary_by_arm %>%
  mutate(
    arm_length_mb = Arm_Length / 1e6,                    # Arm length in Mb
    total_cnv_length = total_cnv_length / 1e6,          # CNV length in Mb
    cnv_rate_per_mb = segment_count / arm_length_mb,    # CNV segments per Mb
    cnv_proportion_per_mb = total_cnv_length / arm_length_mb # Fraction of arm affected
  ) %>%
  mutate(arm = paste0(chrom, Arm))  # Create combined chromosome-arm label

write.table(cnv_summary_final, "cnvs_summary_rate_prop_pq.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------------------
# Split by CNV type
# ----------------------------------------
cnv_split <- split(cnv_summary_final, cnv_summary_final$Gain_Loss)

# ----------------------------------------
# Build matrices per CNV type for rate and proportion
# ----------------------------------------
# Helper to clean matrices
clean_matrix <- function(mat) {
  if (is.null(mat) || !is.matrix(mat)) {
    warning("Input is not a valid matrix.")
    return(NULL)
  }
  mat[is.na(mat) | is.infinite(mat)] <- 0
  return(mat)
}

# Function to generate matrix from CNV data
create_matrix <- function(df, value_col) {
  df_clean <- df %>%
    filter(!is.na(arm), !is.na(ID), !is.na({{value_col}})) %>%
    distinct(arm, ID, .keep_all = TRUE)
  
  df_wide <- df_clean %>%
    select(arm, ID, {{value_col}}) %>%
    pivot_wider(names_from = ID, values_from = {{value_col}}, values_fill = 0)
  
  if (!"arm" %in% colnames(df_wide)) return(NULL)
  
  df_mat <- as.data.frame(df_wide)
  rownames(df_mat) <- df_mat$arm
  df_mat <- df_mat[, -1, drop = FALSE]
  as.matrix(df_mat)
}

rate_matrices_pq <- lapply(cnv_split, create_matrix, value_col = cnv_rate_per_mb)
prop_matrices_pq <- lapply(cnv_split, create_matrix, value_col = cnv_proportion_per_mb)

rate_matrices_pq$Gain <- clean_matrix(rate_matrices_pq$Gain)
rate_matrices_pq$Loss <- clean_matrix(rate_matrices_pq$Loss)
prop_matrices_pq$Gain <- clean_matrix(prop_matrices_pq$Gain)
prop_matrices_pq$Loss <- clean_matrix(prop_matrices_pq$Loss)

# ----------------------------------------
# Heatmap plotting helper
# ----------------------------------------
safe_col_fun <- function(x, colors) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    delta <- if (rng[1] == 0) 1 else abs(rng[1]) * 0.5
    rng <- rng + c(-delta, delta)
  }
  colorRamp2(rng, colors)
}

# Set plot size if using Jupyter/RStudio
options(repr.plot.width = 10, repr.plot.height = 15, repr.plot.res = 200)

# Define color scales
col_fun_gain_rate <- safe_col_fun(rate_matrices_pq$Gain, c("white", "red"))
col_fun_loss_rate <- safe_col_fun(rate_matrices_pq$Loss, c("white", "blue"))
col_fun_gain_prop <- safe_col_fun(log(prop_matrices_pq$Gain + 1), c("white", "red"))
col_fun_loss_prop <- safe_col_fun(log(prop_matrices_pq$Loss + 1), c("white", "blue"))

# Create 2x2 CNV heatmaps
cnv_heatmaps <- list(
  Heatmap(rate_matrices_pq$Gain, name = "CNV_Gain_Rate", col = col_fun_gain_rate,
          column_title = "CNV Gain Rate", cluster_rows = TRUE, cluster_columns = TRUE,
          row_km = 2, na_col = "gray", rect_gp = gpar(col = "gray", lwd = 0.5)),
  Heatmap(log(prop_matrices_pq$Gain + 1), name = "CNV_Gain_Prop", col = col_fun_gain_prop,
          column_title = "CNV Gain Proportion (log)", cluster_rows = TRUE, cluster_columns = TRUE,
          row_km = 2, na_col = "gray", rect_gp = gpar(col = "gray", lwd = 0.5)),
  Heatmap(rate_matrices_pq$Loss, name = "CNV_Loss_Rate", col = col_fun_loss_rate,
          column_title = "CNV Loss Rate", cluster_rows = TRUE, cluster_columns = TRUE,
          row_km = 2, na_col = "gray", rect_gp = gpar(col = "gray", lwd = 0.5)),
  Heatmap(log(prop_matrices_pq$Loss + 1), name = "CNV_Loss_Prop", col = col_fun_loss_prop,
          column_title = "CNV Loss Proportion (log)", cluster_rows = TRUE, cluster_columns = TRUE,
          row_km = 2, na_col = "gray", rect_gp = gpar(col = "gray", lwd = 0.5))
)

# Draw in 2x2 grid layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
for (i in seq_along(cnv_heatmaps)) {
  row <- ceiling(i / 2)
  col <- ifelse(i %% 2 == 1, 1, 2)
  vp <- viewport(layout.pos.row = row, layout.pos.col = col)
  pushViewport(vp)
  draw(cnv_heatmaps[[i]], newpage = FALSE, heatmap_legend_side = "bottom")
  upViewport()
}

# ----------------------------------------
# Save heatmaps as PDF
# ----------------------------------------
output_dir <- "/home/mjehangir/telomere-sv-analysis/plots/fig2/"
heatmap_titles <- c("gain_rate_heatmap.pdf", "gain_prop_heatmap.pdf", 
                    "loss_rate_heatmap.pdf", "loss_prop_heatmap.pdf")

for (i in seq_along(cnv_heatmaps)) {
  pdf(file.path(output_dir, heatmap_titles[i]), width = 8, height = 8)
  draw(cnv_heatmaps[[i]], newpage = FALSE, heatmap_legend_side = "bottom")
  dev.off()
}
