

#### Nephelometry Data Analysis Script ####

## Etienne Bremand ###

# R version 4.3.1 #



#### Package Management Script ####
# Description: Checks, installs (if needed), and loads required R packages for nephelometry analysis.

```{r}
packages <- c(
  "dplyr",       # v1.1.4
  "ggplot2",     # v4.0.0
  "grid",        # v4.3.1
  "gridExtra",   # v2.3
  "pracma",      # v2.4.6
  "stringi",     # v1.8.7
  "stringr",     # v1.5.2
  "tidyverse",   # v2.0.0
  "tidyr",       # v1.3.1
  "writexl"      # v1.5.3
)

# Install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("→ Installing package:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all packages quietly
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

cat("\nAll packages successfully loaded.\n")

```




#### Nephelometry Data Import Script ####
# Description: Sets working directory, imports nephelometry data, and filters for a specific strain.

```{r}
################################# TO MODIFY ##########################################

working_directory <- ""


data_nephelo <- ""

######################################################################################

# Set working directory
setwd(working_directory)

# Import data and filter for the chosen strain

data <- read.table(data_nephelo,
                   sep = "\t",
                   dec = ",",
                   header = TRUE,
                   fill = TRUE,
                   quote = "",
                   check.names = FALSE)

# 2. Vérification immédiate
cat("Nombre de colonnes chargées :", ncol(data), "\n")
head(colnames(data), 10) # Affiche les 10 premiers noms de colonnes
#data <- data[data$Souche == souche,]

# Display first rows for verification
#data <- data[data$Appareil == "Nephelo",]

head(data)

```


#### Nephelometry Parameters Setup ####
# Description: Define columns for measurements, conditions, replicates, and timing between measurements.

```{r}
################################# TO MODIFY ##########################################

beginning <-    # Column number of the first measurement
condition <-   # Column number containing condition names to group
rep <-          # Column number for replicates

timing <- 10     # Time interval between measurements in minutes

######################################################################################
```




#### Data Cleaning and Time Setup ####
# Description: Keep non-empty measurement columns, rename them in minutes, and calculate time vectors per condition.

```{r}
col.condition <- condition
col.measurements <- beginning:ncol(data)

# Keep only columns with at least one non-NA value
data_meas <- data[, col.measurements]
non_na_cols <- which(colSums(!is.na(data_meas)) > 0)
col.measurements <- beginning - 1 + non_na_cols
data <- data[, c(1:(beginning-1), col.measurements)]

# Rename measurement columns in minutes
n_meta <- beginning - 1
col_meas <- (n_meta + 1):ncol(data)
minutes_vect <- seq(0, by = timing, length.out = length(col_meas))
colnames(data)[col_meas] <- paste0("min_", minutes_vect)

# Check if all measurement columns are already numeric
initially_numeric <- all(sapply(data[, col_meas], is.numeric))

if (initially_numeric) {
  cat("Check Numeric: OK - All measurement columns are numeric.\n")
} else {
  # Attempt conversion
  data[, col_meas] <- lapply(data[, col_meas], function(x) as.numeric(as.character(x)))
  
  # Verify again after conversion
  still_not_numeric <- any(!sapply(data[, col_meas], is.numeric))
  
  if (still_not_numeric) {
    cat("Check Numeric: Error - Even after conversion, some values are not numeric.\n")
  } else {
    cat("Check Numeric: Fixed - Values were not numeric but have been correctly changed.\n")
  }
}

# Extract all and unique conditions
conditions.toutes <- data[[col.condition]]
conditions.uniques <- unique(conditions.toutes)

# Time vector in hours
temps <- seq(0, (length(col_meas) - 1) * timing / 60, by = timing / 60)

# ---- Max time per row ----
max.time.par.ligne <- apply(data[, col_meas], 1, function(x) {
  # last.valid is the index of the last non-NA value
  last.valid <- if(all(is.na(x))) NA else max(which(!is.na(x)))
  if (is.finite(last.valid)) temps[last.valid] else NA
})

min.time.max <- min(max.time.par.ligne, na.rm = TRUE)
max.time.max <- max(max.time.par.ligne, na.rm = TRUE)

# Summary info (kept at the end)
cat("\n--- Condition Info ---\n")
cat("Total rows:", nrow(data), "\n")
cat("Distinct conditions:", length(conditions.uniques), "\n")
cat("Example conditions:", paste(head(conditions.uniques, 5), collapse = ", "), "\n\n")

cat("--- Growth Duration Analysis ---\n")
cat("Minimum duration among all rows:", round(min.time.max, 2), "hours.\n")
cat("Maximum duration among all rows:", round(max.time.max, 2), "hours.\n")

# ---- Output Directory Setup ----
output_dir <- "output_result"

if (!file.exists(output_dir)) {
  dir.create(output_dir)
  cat("\nOutput Directory: Created successfully at '", output_dir, "'.\n", sep = "")
} else {
  cat("\nOutput Directory: OK - '", output_dir, "' already exists.\n", sep = "")
}

```










#### Curve Smoothing and Plateau Fixing ####
# Description: Corrects initial fog spike, sets baseline to 0, optionally smooths curves (moving average or LOESS), 
#              and blocks values after maximum. Generates line plots before/after plateau fixation.

```{r}
################################# TO MODIFY ##########################################

smooth.curves  <- TRUE        # Smooth curves if TRUE
smooth.window  <- 10           # Window size in HOURS (e.g., 5 means smoothing over a 5h neighborhood)

export_graphs  <- TRUE        # TRUE = export plots in "graph_after_before_cleaning" folder

######################################################################################

cols_mesures <- beginning:ncol(data)
mesures <- data[, cols_mesures]

# ---- 0. Sauvegarde des données strictement brutes ----
data_brut_total <- data 

# ---- 1 & 2. Correction synchronisée (Fog Spike + Baseline) ----
mesures_corrected <- t(apply(mesures, 1, function(x){
  x <- as.numeric(x)
  if(all(is.na(x))) return(x)
  search_limit <- floor(length(x) / 2)
  min_val <- min(x[1:search_limit], na.rm=TRUE)
  min_idx <- which.min(x[1:search_limit])
  x[1:min_idx] <- min_val
  x_final <- x - min_val
  x_final[x_final < 0] <- 0
  return(x_final)
}))

# ---- 3. Préparation des datasets pour comparaison ----
data2 <- cbind(data[,1:(beginning-1)], mesures_corrected)
colnames(data2) <- colnames(data)

# ---- 4. Lissage LOESS avec fenêtre en heures (si TRUE) ----
mesures_temp <- data2[, beginning:ncol(data2)]
if(smooth.curves){
  # Calcul de la durée totale de l'expérience en heures
  duree_totale_h <- max(temps, na.rm = TRUE)
  
  # Calcul du span (proportion) basé sur la fenêtre en heures
  # On s'assure que le span est compris entre 0.05 et 1
  loess_span <- max(0.05, min(smooth.window / duree_totale_h, 1))
  
  mesures_temp <- t(apply(mesures_temp, 1, function(x){
    if(all(is.na(x))) return(x)
    df <- data.frame(y=x, x=seq_along(x))
    # Application du LOESS
    mod <- try(loess(y~x, data=df, span=loess_span, degree=2, na.action=na.exclude), silent=TRUE)
    if(inherits(mod,"try-error")) x else predict(mod, newdata=df$x)
  }))
}

# data3 contient les données corrigées et lissées
data3 <- cbind(data2[,1:(beginning-1)], mesures_temp)
colnames(data3) <- colnames(data)

# ---- 8.5 Calcul de l'échelle Y globale (kNTU) ----
unit_factor <- 1000
val_max_globale <- max(c(as.matrix(data_brut_total[, beginning:ncol(data_brut_total)]), 
                         as.matrix(data3[, beginning:ncol(data3)])), 
                       na.rm = TRUE) / unit_factor

y_limit_fixe <- val_max_globale * 1.05

# ---- Setup Export Directory ----
if (export_graphs) {
  setwd(working_directory)
  graph_dir <- file.path("output_result", "graph_after_before_cleaning")
  if (!file.exists(graph_dir)) {
    dir.create(graph_dir, recursive = TRUE)
    cat("Graph Directory: Created successfully at '", graph_dir, "'.\n", sep = "")
  }
}

# ---- 9. Plot: Total Raw vs Processed (Unit: kNTU) ----
right_title <- paste0("Processed (Fog+Base", ifelse(smooth.curves, "+LOESS", ""), ")")

for(cond in conditions.uniques){
  idx <- which(data[[col.condition]] == cond)
  
  make_df_with_labels <- function(dataset) {
    df_list <- lapply(seq_along(idx), function(i){
      li <- idx[i]
      temp <- dataset[li, beginning:ncol(dataset)]
      valid <- !is.na(temp)
      if(!any(valid)) return(NULL)
      
      val_rep <- as.character(dataset[li, rep])
      
      data.frame(
        LineID = as.factor(li), 
        Rep    = as.factor(val_rep),
        Time   = temps[valid], 
        Value  = as.numeric(temp[valid]) / unit_factor
      )
    })
    
    df <- do.call(rbind, df_list)
    if(is.null(df)) return(NULL)
    
    labels_df <- df %>% group_by(LineID) %>% slice_tail(n = 1)
    return(list(data = df, labels = labels_df))
  }

  res_raw  <- make_df_with_labels(data_brut_total)
  res_proc <- make_df_with_labels(data3)

  if(is.null(res_raw)) next

  # Graphique de gauche (Brut)
  p1 <- ggplot(res_raw$data, aes(x = Time, y = Value, group = LineID, color = Rep)) +
    geom_line(alpha = 0.6, linewidth = 0.8) +
    geom_text(data = res_raw$labels, aes(label = LineID), hjust = 0, nudge_x = 0.5, size = 3, color = "black") +
    ggtitle(paste0(cond, " - Raw Data")) + 
    xlab("Time (h)") + ylab("kNTU") +
    ylim(0, y_limit_fixe) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + 
    scale_color_brewer(palette = "Set1") +
    theme_minimal() + theme(legend.position = "bottom")

  # Graphique de droite (Traité)
  p2 <- ggplot(res_proc$data, aes(x = Time, y = Value, group = LineID, color = Rep)) +
    geom_line(alpha = 0.6, linewidth = 0.8) +
    geom_text(data = res_proc$labels, aes(label = LineID), hjust = 0, nudge_x = 0.5, size = 3, color = "black") +
    ggtitle(paste0(cond, " - ", right_title)) + 
    xlab("Time (h)") + ylab("kNTU") +
    ylim(0, y_limit_fixe) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + 
    scale_color_brewer(palette = "Set1") +
    theme_minimal() + theme(legend.position = "bottom")

  # ---- Export or Display Control ----
  if (export_graphs) {
    clean_cond_name <- gsub("[^A-Za-z0-9_-]", "_", cond)
    file_name <- file.path(graph_dir, paste0("", clean_cond_name, ".png"))
    
    # 1. Open png device and plot to save file
    png(file_name, width = 1000, height = 500, res = 100)
    grid.arrange(p1, p2, ncol = 2, widths = c(1, 1))
    dev.off() 
    
    # 2. Re-plot immediately to show it inside the RMarkdown script
    grid.arrange(p1, p2, ncol = 2, widths = c(1, 1))
    
  } else {
    # If FALSE, just print it in RStudio / HTML output
    grid.arrange(p1, p2, ncol = 2, widths = c(1, 1))
  }
}
```






#### Line removing and Graph Export ####
# Description: Removes selected lines, re-fixes the plateau for specific lines, 
#              and generates line plots by condition. Plots are displayed and optionally exported.

```{r}
################################# TO MODIFY ##########################################

lines_to_remove   <- c() # Lines to remove

export_graphs     <- TRUE                        # TRUE = export plots
export_data_clean <- TRUE                        # TRUE = export data_ntu_clean.xlsx

######################################################################################

setwd(working_directory)

# ---- 1. On repart de data3 (la sortie de ton bloc précédent) ----
data_source <- data3
cols_mesures <- beginning:ncol(data_source)

# On garde une trace des IDs originaux avant de supprimer des lignes
data_source$Original_LineID <- 1:nrow(data_source)

# Correction des noms de colonnes
col.condition.name <- colnames(data_source)[condition]
col.rep.name       <- colnames(data_source)[rep]

# ---- 2. Création de data4 (Le dataset "Clean") ----
lines_to_keep <- setdiff(1:nrow(data_source), lines_to_remove)
data4 <- data_source[lines_to_keep, ]

cat("Cleaning completed:\n")
cat("- Removed lines:", length(lines_to_remove), "\n")
cat("- Remaining lines (data4):", nrow(data4), "\n")

# ---- Setup Final Graph Directory ----
if (export_graphs) {
  # Build path: output_result/final_graph
  final_graph_dir <- file.path(output_dir, "final_growth_curve")
  
  if (!file.exists(final_graph_dir)) {
    dir.create(final_graph_dir, recursive = TRUE)
    cat("Graph Directory: Created successfully at '", final_graph_dir, "'.\n", sep = "")
  }
}

# ---- 3. Plot des courbes nettoyées (Format Large & kNTU) ----

unit_factor <- 1000
# On réutilise la limite Y fixe calculée au bloc précédent pour la cohérence
y_limit_clean <- max(data4[, cols_mesures], na.rm=TRUE) / unit_factor * 1.05

for(cond in unique(data4[[col.condition.name]])){
  
  idx <- which(data4[[col.condition.name]] == cond)
  
  df_plot <- do.call(rbind, lapply(idx, function(li){
    temp <- as.numeric(data4[li, cols_mesures])
    valid <- !is.na(temp)
    data.frame(
      LineID = data4$Original_LineID[li],
      Rep    = as.factor(data4[li, col.rep.name]),
      Time   = temps[valid],
      Value  = temp[valid] / unit_factor # Passage en kNTU
    )
  }))
  
  if(is.null(df_plot) || nrow(df_plot) == 0) next
  
  df_last <- df_plot %>% group_by(LineID) %>% slice_tail(n = 1)
  
  p <- ggplot(df_plot, aes(x = Time, y = Value, group = LineID, color = Rep)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_text(data = df_last, aes(label = LineID), 
              hjust = 0, nudge_x = 0.5, size = 3, color = "black") +
    ggtitle(paste0(cond, "")) + 
    xlab("Time (h)") + ylab("kNTU") +
    ylim(0, y_limit_fixe) + # On utilise la même échelle que le bloc de comparaison
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_color_brewer(palette = "Set1") + 
    theme_minimal() + 
    theme(legend.position = "bottom")

  # Affichage
  print(p)
  
  if(export_graphs){
    setwd(working_directory)
    # Target the new final_graph directory
    output_file <- file.path(final_graph_dir, paste0("clean_graph_", gsub("[^A-Za-z0-9]", "_", cond), ".png"))
    ggsave(filename = output_file, plot = p, width = 10, height = 6, dpi = 300)
  }
}

# ---- 4. Export Excel  ----
if(export_data_clean){
  setwd(working_directory)
  library(writexl)
  # Build the clean output path: output_result/data_ntu_clean.xlsx
  excel_output_path <- file.path(output_dir, "data_ntu_clean.xlsx")
  write_xlsx(data4, path = excel_output_path)
  cat("File exported successfully:", excel_output_path, "\n")
}

```





#### Definition of thresholds for lag phase and AUC calculation ####
# Description: Set the NTU threshold to determine the start of the lag phase
#              and the time threshold up to which the AUC will be calculated.
#              Optionally, export the mean curves with thresholds as a PNG file.

```{r}

################################# TO MODIFY ##########################################

max_time     <- 100 # Maximum time threshold (hours) for all subsequent calculations
export_graph <- TRUE # TRUE = export plots, FALSE = just display

#######################################################################################

# ---- Create summary table with mean and standard error per condition ----
cols_mesures <- beginning:ncol(data4)
unit_factor <- 1000 # Factor to convert NTU to kNTU

df_long <- do.call(rbind, lapply(conditions.uniques, function(cond) {

  idx <- which(data4[[col.condition]] == cond)
  
  if(length(idx) == 0) return(NULL)

  temp_mat <- as.matrix(data4[idx, cols_mesures])
  valid_cols <- colSums(!is.na(temp_mat)) > 0
  
  if(!any(valid_cols)) return(NULL)

  temp_mat <- temp_mat[, valid_cols, drop = FALSE]
  temps_valid <- temps[valid_cols]

  # Mean and SE (divided by unit_factor to switch to kNTU)
  mean_vals <- colMeans(temp_mat, na.rm = TRUE) / unit_factor
  se_vals   <- apply(temp_mat, 2, function(x) {
    s <- sd(x, na.rm = TRUE)
    n <- sum(!is.na(x))
    if(is.na(s) | n == 0) return(0) else return(s/sqrt(n))
  }) / unit_factor

  if(length(temps_valid) == length(mean_vals)) {
    data.frame(
      Condition = cond,
      Time = temps_valid,
      Mean = mean_vals,
      SE   = se_vals
    )
  } else {
    NULL
  }
}))

# ---- Plot mean curves with SE and analysis time limit ----
p <- ggplot(df_long, aes(x = Time, y = Mean, color = Condition, fill = Condition)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE), alpha = 0.2, color = NA) +
  # Vertical line indicating the analysis limit
  geom_vline(xintercept = max_time, linetype = "dotted", color = "black") +
  xlab("Time (h)") + ylab("kNTU") + # Axe Y modifié en kNTU
  ggtitle(paste0("Mean curves per condition")) +
  theme_minimal() +
  theme(legend.position = "none")

# ---- Export plot if requested ----
if(export_graph){
  setwd(working_directory)
  output_file <- file.path(output_dir, "Mean_curves_analysis_limit.png")
  ggsave(filename = output_file, plot = p, width = 8, height = 5, dpi = 300)
  cat("File exported successfully:", output_file, "\n")
}

# ---- Display the plot ----
suppressWarnings(print(p))

```








#### Maximal Slope Calculation ####
# Description: Computes the maximal slope for each growth curve between the lag phase threshold 
#              and the plateau. Uses a fixed sliding window of specified size (window_size).
#              Generates plots per condition showing individual curves, lag phase threshold,
#              and the segment with maximal slope. Graphs can optionally be exported.

```{r}
################################# TO MODIFY ##########################################

window_size   <- 10           # Sliding window size (hours)
export_graphs <- TRUE        

######################################################################################

cols_mesures <- beginning:ncol(data4)
mat_data <- as.matrix(data4[, cols_mesures])
temps_vect <- temps[1:ncol(mat_data)]
pts_window <- round(window_size / (timing / 60))

# Fonction de pente ultra-rapide
fast_slope <- function(x, y) {
  n <- length(x)
  if(n < 3) return(NA)
  denominator <- (n * sum(x^2) - (sum(x))^2)
  if(is.na(denominator) || denominator == 0) return(NA)
  return((n * sum(x * y) - sum(x) * sum(y)) / denominator)
}

pentes_list <- list()

for(i in 1:nrow(mat_data)) {
  y_all <- mat_data[i, ]
  mask <- !is.na(y_all) & temps_vect <= max_time
  y <- y_all[mask]
  t <- temps_vect[mask]
  
  if(length(y) < pts_window || length(t) < pts_window) next
  
  max_val_courbe <- max(y, na.rm = TRUE)
  idx_max_absolu <- which(y == max_val_courbe)[1]
  t_max_absolu <- t[idx_max_absolu]
  
  n_windows <- idx_max_absolu - pts_window + 1
  if(n_windows < 1) n_windows <- length(y) - pts_window + 1
  
  slopes <- sapply(1:max(1, n_windows), function(j) {
    idx <- j:(j + pts_window - 1)
    if(max(idx) <= length(t)) fast_slope(t[idx], y[idx]) else NA
  })
  
  if(all(is.na(slopes))) next
  best_idx <- which.max(slopes)
  max_pente <- slopes[best_idx]
  
  t_win <- t[best_idx:(best_idx + pts_window - 1)]
  y_win <- y[best_idx:(best_idx + pts_window - 1)]
  intercept <- mean(y_win) - max_pente * mean(t_win)
  
  t_lag_geo <- -intercept / max_pente
  t_plateau_geo_projet <- (max_val_courbe - intercept) / max_pente
  
  if(is.na(t_plateau_geo_projet) || is.na(t_lag_geo)) {
    t_plateau_geo_final <- t_max_absolu
  } else if (t_plateau_geo_projet > max(t, na.rm=TRUE) * 1.5 || t_plateau_geo_projet < t_lag_geo) {
    t_plateau_geo_final <- t_max_absolu 
  } else {
    t_plateau_geo_final <- t_plateau_geo_projet 
  }
  
  pentes_list[[length(pentes_list) + 1]] <- data.frame(
    LineID = i,
    Original_LineID = data4$Original_LineID[i],
    slope = max_pente,
    intercept = intercept,
    t_lag_geo = max(0, t_lag_geo, na.rm = TRUE),
    t_plateau_geo = t_plateau_geo_final,
    max_val = max_val_courbe,
    t1_plot = t_win[1],
    t2_plot = t_win[length(t_win)],
    y1_plot = max_pente * t_win[1] + intercept,
    y2_plot = max_pente * t_win[length(t_win)] + intercept
  )
}

pentes_df <- do.call(rbind, pentes_list)

# ---- VISUALISATION (Échelles Fixes basées sur le max des courbes) ----

unit_factor <- 1000 # Factor to convert NTU to kNTU
global_max_y <- (max(mat_data, na.rm = TRUE) / unit_factor) * 1.05  # +5% margin in kNTU
global_max_x <- max(temps_vect, na.rm = TRUE)

# ---- Setup Verif_slope Directory inside output_result ----
output_dir_slope <- file.path(output_dir, "Verif_slope")
if(export_graphs && !dir.exists(output_dir_slope)) {
  dir.create(output_dir_slope, recursive = TRUE)
  cat("Slope Verification Directory: Created successfully at '", output_dir_slope, "'.\n", sep = "")
}

for(cond in unique(data4[[col.condition]])) {
  idx_rows <- which(data4[[col.condition]] == cond)
  
  df_plot <- do.call(rbind, lapply(idx_rows, function(r) {
    val <- as.numeric(data4[r, cols_mesures])
    mask <- !is.na(val)
    data.frame(LineID = r, Time = temps_vect[mask], Value = val[mask] / unit_factor) # Conversion to kNTU
  }))
  
  # Conversion of slopes and intercepts parameters for the kNTU scale visualization
  df_pentes_cond <- data.frame()
  if(!is.null(pentes_df)) {
    df_pentes_cond <- pentes_df[pentes_df$LineID %in% idx_rows, ]
    if(nrow(df_pentes_cond) > 0) {
      df_pentes_cond$slope     <- df_pentes_cond$slope / unit_factor
      df_pentes_cond$intercept <- df_pentes_cond$intercept / unit_factor
      df_pentes_cond$max_val   <- df_pentes_cond$max_val / unit_factor
      df_pentes_cond$y1_plot   <- df_pentes_cond$y1_plot / unit_factor
      df_pentes_cond$y2_plot   <- df_pentes_cond$y2_plot / unit_factor
    }
  }

  p <- ggplot(df_plot, aes(x = Time, y = Value, group = LineID)) +
    geom_line(alpha = 0.2, color = "black") + 
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") + 
    geom_hline(data = df_pentes_cond, aes(yintercept = max_val), linetype = "dotted", color = "blue") +
    geom_abline(data = df_pentes_cond, aes(slope = slope, intercept = intercept), 
                linetype = "dashed", color = "orange", alpha = 0.4) +
    geom_segment(data = df_pentes_cond, 
                 aes(x = t1_plot, xend = t2_plot, y = y1_plot, yend = y2_plot), 
                 linewidth = 1.5, color = "orange") +
    geom_point(data = df_pentes_cond, aes(x = t_lag_geo, y = 0), color = "red", size = 2.5) +
    geom_point(data = df_pentes_cond, aes(x = t_plateau_geo, y = max_val), color = "blue", size = 2.5) +
    coord_cartesian(xlim = c(0, global_max_x), ylim = c(-0.02 * global_max_y, global_max_y)) +
    ggtitle(paste0(cond)) + xlab("Time (h)") + ylab("kNTU") + # Labels adjusted
    theme_minimal()
  
  print(p)
  
  if(export_graphs){
    setwd(working_directory)
    clean_name <- gsub("[^A-Za-z0-9_-]", "_", cond)
    output_file <- file.path(output_dir_slope, paste0("slope_check_", clean_name, ".png"))
    ggsave(filename = output_file, plot = p, width = 10, height = 6, dpi = 300)
  }
}
```


















#### Kinetic Parameter Calculation ####
# Description: Computes key kinetic parameters for each line including:
#              - AUC: Area Under the Curve until the max_time threshold
#              - Max_slope: Maximum growth rate (maximal sliding window slope)
#              - Maximum_capacity: Highest turbidity value reached (absolute plateau)
#              - Lag_time: Geometric germination time (tangent intersection with y=0)
#              - Plateau_time: Geometric plateau entry time (tangent intersection with y=max)
#              - Exponential_duration: Calculated duration of active growth phase
# The results can optionally be exported to an Excel file.

```{r}

################################# TO MODIFY ##########################################

export_results <- TRUE   # TRUE = export results to Excel, FALSE = just display head()

######################################################################################

# ---- Measurement columns ----
cols_mesures <- beginning:ncol(data4)

# ---- Main calculation function (Geometric & Final) ----
calc_params_final <- function(i, time, data_matrix, pentes_df, max_time_limit) {
  
  values <- as.numeric(data_matrix[i, ])
  valid <- !is.na(values)
  t_vec <- time[valid]
  v_vec <- values[valid]
  
  # 1. Récupération des données de la tangente pour cette ligne
  orig_id <- data4$Original_LineID[i]
  pente_info <- pentes_df[pentes_df$Original_LineID == orig_id, ]
  
  if (nrow(pente_info) == 0) return(rep(NA, 6)) # 6 paramètres restants
  
  # Paramètres géométriques
  t_lag_geo     <- pente_info$t_lag_geo
  t_plateau_geo <- pente_info$t_plateau_geo
  max_pente     <- pente_info$slope
  max_ntu       <- pente_info$max_val
  
  # --- 2. AUC Totale (jusqu'à max_time) ---
  idx_auc_tot <- which(t_vec <= max_time_limit)
  auc_total <- if(length(idx_auc_tot) > 1) pracma::trapz(t_vec[idx_auc_tot], v_vec[idx_auc_tot]) else NA
  
  # --- 3. Durée de la phase exponentielle ---
  exp_duration <- t_plateau_geo - t_lag_geo
  
  return(c(
    auc_total,      # AUC totale
    max_pente,      # Pente max
    max_ntu,        # Max NTU
    t_lag_geo,      # Temps de germination
    t_plateau_geo,  # Temps pour arriver au plateau
    exp_duration    # Durée phase expo
  ))
}

# ---- Apply function to all rows ----
resultats_raw <- vapply(
  1:nrow(data4),
  function(i) calc_params_final(i, temps_vect, mat_data, pentes_df, max_time),
  numeric(6)
)
resultats_final <- t(resultats_raw)

# ---- Assign explicit column names ----
colnames(resultats_final) <- c(
  "AUC",                 # AUC totale jusqu'à max_time
  "Max_slope",           # Pente maximale
  "Maximum_capacity",    # Valeur NTU max atteinte
  "Lag_time",            # Temps de germination (Intersection y=0)
  "Plateau_time",        # Temps de plateau (Intersection y=max)
  "Exponential_duration" # Durée de la phase expo
)

# ---- Combine with metadata ----
meta_cols <- setdiff(names(data4), names(data4)[cols_mesures])
resultats_cinetiques <- cbind(data4[, meta_cols, drop = FALSE], resultats_final)

# ---- Export results if requested ----
setwd(working_directory)

if(export_results){
  # Target the global output_result directory
  output_file <- file.path(output_dir, "kinetic_parameters_final.xlsx")
  writexl::write_xlsx(resultats_cinetiques, path = output_file)
  cat("Export of parameters for each sample successfully saved at:", output_file, "\n")
}

head(resultats_cinetiques)
```


















#### Summary and Graphical Representation of Kinetic Parameters ####
# Description: Generates summary tables of kinetic parameters per condition and per replicate,
#              optionally exports them to Excel, and produces barplots or boxplots
#              for each kinetic parameter. User can select error type (SD or SE) and graph type.

```{r}

################################# TO MODIFY ##########################################

erreur_type         <- "se"       # "sd" or "se"
graph_type          <- "barplot"  # "barplot" or "boxplot"

export_by_rep       <- TRUE    # TRUE = export tables by replicate
export_by_condition <- TRUE    # TRUE = export tables by condition
export_graphs       <- TRUE    # TRUE = export graphs to folder

######################################################################################

# 1. Identification des colonnes
col_condition <- colnames(resultats_cinetiques)[condition]

if(!is.na(rep) && rep <= ncol(resultats_cinetiques)) {
  col_replicat <- colnames(resultats_cinetiques)[rep]
} else {
  col_replicat <- NA
}

# On identifie les paramètres cinétiques
parametres_existants <- c("AUC", "Max_slope", "Maximum_capacity", "Lag_time", "Plateau_time", "Exponential_duration")

# ---- Summary by condition ----
resume_condition <- resultats_cinetiques %>%
  group_by(across(all_of(col_condition))) %>%
  summarise(
    n = n(),
    across(
      all_of(parametres_existants),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE),
        se   = ~sd(.x, na.rm = TRUE)/sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Filtrage des colonnes selon erreur_type pour plus de clarté
suffixe_err <- paste0("_", erreur_type)
resume_condition <- resume_condition %>%
  select(all_of(col_condition), n, ends_with("_mean"), ends_with(suffixe_err))

# ---- Summary by replicate (si applicable) ----
if(!is.na(col_replicat)) {
  resultats_replicats <- resultats_cinetiques %>%
    group_by(across(all_of(c(col_condition, col_replicat)))) %>%
    summarise(
      n = n(),
      across(
        all_of(parametres_existants),
        list(
          mean = ~mean(.x, na.rm=TRUE),
          sd   = ~sd(.x, na.rm=TRUE),
          se   = ~sd(.x, na.rm=TRUE)/sqrt(sum(!is.na(.x)))
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) %>%
    select(all_of(c(col_condition, col_replicat)), n, ends_with("_mean"), ends_with(suffixe_err))
} else {
  resultats_replicats <- resultats_cinetiques
}

setwd(working_directory)
# ---- Prepare output folder for graphs inside output_result ----
if(export_graphs) {
  output_dir_param <- file.path(output_dir, "Graphs_by_parameter")
  if(!dir.exists(output_dir_param)) {
    dir.create(output_dir_param, recursive = TRUE)
    cat("Graph Parameter Directory: Created successfully at '", output_dir_param, "'.\n", sep = "")
  }
}

# ---- Graphs for each parameter ----
for(param in parametres_existants) {
  y_mean <- paste0(param, "_mean")
  y_err  <- paste0(param, "_", erreur_type)
  
  if(graph_type == "barplot") {
    p <- ggplot(resume_condition, aes(x = .data[[col_condition]], y = .data[[y_mean]])) +
      geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
      geom_errorbar(aes(
        ymin = .data[[y_mean]] - .data[[y_err]],
        ymax = .data[[y_mean]] + .data[[y_err]]
      ), width = 0.2) +
      xlab("Condition") + ylab(param) +
      ggtitle(paste("Mean", param, "by condition (", erreur_type, ")")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else if(graph_type == "boxplot") {
    p <- ggplot(resultats_cinetiques, aes(x = .data[[col_condition]], y = .data[[param]])) +
      geom_boxplot(fill = "lightblue", color = "black") +
      xlab("Condition") + ylab(param) +
      ggtitle(paste("Boxplot of", param, "by condition")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  print(p)
  
  if(export_graphs){
    file_name <- file.path(output_dir_param, paste0(gsub("[^A-Za-z0-9]", "_", param), ".png"))
    ggsave(filename = file_name, plot=p, width=8, height=5, dpi=300)
  }
}

setwd(working_directory)

# ---- Export Excel to output_result ----
if(export_by_rep) {
  excel_rep_path <- file.path(output_dir, "parameters_by_replicate.xlsx")
  write_xlsx(resultats_replicats, path = excel_rep_path)
  cat("File exported successfully:", excel_rep_path, "\n")
}

if(export_by_condition) {
  excel_cond_path <- file.path(output_dir, "parameters_by_condition.xlsx")
  write_xlsx(resume_condition, path = excel_cond_path)
  cat("File exported successfully:", excel_cond_path, "\n")
}

# Affichage final
resume_condition

```



