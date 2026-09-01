

# Script for generating figures for the publication:

# Comparative Evaluation of Growth Assessment Methods in Filamentous Fungi Reveals the Potential of Nephelometry for Nutrient Profiling


# Load dataset
```{r}

data <- read.delim2("")

```


# Load required packages
```{r}

if(!require(tidyverse)) install.packages("tidyverse")
if(!require(GGally)) install.packages("GGally") 
if(!require(Hmisc)) install.packages("Hmisc")
if(!require(corrplot)) install.packages("corrplot")
if(!require(patchwork)) install.packages("patchwork")
if(!require(cowplot)) install.packages("cowplot")
if(!require(multcompView)) install.packages("multcompView")

library(tidyverse)
library(GGally)
library(Hmisc)
library(cowplot)
library(patchwork)
library(multcompView)

```


# Load colors assigned to each nutrient
```{r}

custom_colors <- c(
  "Ammonium"      = "#C3223B",
  "Asparagine"    = "#EF6D6D",
  "Control"       = "#000000",
  "GABA"          = "#E5E87B",
  "Glutamate"     = "#0BA279",
  "Glycine"       = "#69955C",
  "Nitrate"       = "#51FF00",
  "Nitrite"       = "#73A6F0",
  "Phenylalanine" = "#3B8ED2",
  "Proline"       = "#5200FF",
  "Serine"        = "#9B8BDE",
  "Taurine"       = "#CEB551",
  "Tryptophan"    = "#808D85"
)

```








# Data cleaning and normality testing

```{r}

df_clean_wide <- data %>%
  # 1. Clean numeric values
  mutate(Valeur = as.numeric(str_replace_all(Valeur, ",", "."))) %>%
  # Exclude Arginine / Serine if necessary
  filter(!Condition %in% c("Arginine", "Serine")) %>%
  
  # 2. Construct target variable names
  mutate(
    Suffixe_Time = if_else(is.na(Timing) | Timing == "NA", "", paste0("_", Timing)),
    Var_Name = if_else(
      Appareil == "Boite",
      Parameter,
      paste0(Parameter, Suffixe_Time)
    ),
    # Clean special characters in variable names
    Var_Name = str_replace_all(Var_Name, "[- ]", "_")
  ) %>%
  
  # 3. Handle duplicates / technical replicates before pivoting
  group_by(Condition, Rep, Var_Name) %>%
  summarise(Valeur = mean(Valeur, na.rm = TRUE), .groups = "drop") %>%
  
  # 4. Convert to wide format
  pivot_wider(names_from = Var_Name, values_from = Valeur) %>%
  rename(Nitrogen_source = Condition, Biological_replicate = Rep) %>%
  mutate(Nitrogen_source = as.factor(Nitrogen_source)) %>%
  
  # 5. Divide nephelometric MCC values by 1000 without changing column names
  mutate(across(contains("MCC") & contains("Nephelo"), ~ . / 1000))

# Aggregate at the biological level (mean for each nitrogen source)
df_agreg_azote <- df_clean_wide %>%
  select(-Biological_replicate) %>% 
  group_by(Nitrogen_source) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = "drop")


# Quick normality assessment for key variables
vars_num <- df_agreg_azote %>% select(where(is.numeric)) %>% names()

normality_results <- map_dfr(vars_num, function(v) {
  x <- na.omit(df_agreg_azote[[v]])
  if(length(x) >= 3) {
    sh <- shapiro.test(x)
    tibble(Variable = v, W = round(sh$statistic, 4), P_Value = sh$p.value)
  }
})

cat("\n=== NORMALITY TEST RESULTS (SHAPIRO-WILK) ===\n")
print(normality_results)

```



# Figure 2
# Correlation matrix at 33 h for nephelometry

```{r}
library(tidyverse)
library(Hmisc)

# 1. Data preparation
df_corr <- data %>%
  filter(
    (Appareil == "Nephelo" & Timing == "33") |
      (Appareil == "Boite" & Parameter == "Dry_weight")
  ) %>%
  mutate(
    Valeur = as.numeric(Valeur),
    Parameter = case_when(
      Parameter == "Nephelo-AUC" ~ "AUC",
      Parameter == "Nephelo-MCC" ~ "MCC",
      Parameter == "Nephelo-Growth_rate" ~ "Growth Rate",
      Parameter == "Nephelo-Lag_time" ~ "Lag Time",
      Parameter == "Nephelo-Plateau_time" ~ "MCC Time",
      Parameter == "Nephelo-Exponential_duration" ~ "Exponential duration",
      Parameter == "Dry_weight" ~ "Dry weight"
    )
  ) %>%
  group_by(Condition, Parameter) %>%
  summarise(Valeur = mean(Valeur, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Parameter, values_from = Valeur) %>%
  select(any_of(c(
    "Dry weight", "AUC", "MCC", "Growth Rate",
    "Lag Time", "MCC Time", "Exponential duration"
  )))

# 2. Spearman correlations
ct <- rcorr(as.matrix(df_corr), type = "spearman")

# 3. Heatmap data
ordre <- c(
  "Dry weight", "AUC", "MCC", "Growth Rate",
  "Lag Time", "MCC Time", "Exponential duration"
)

matrice_data <- as.data.frame(as.table(ct$r)) %>%
  rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
  mutate(p = as.vector(ct$P)) %>%
  filter(match(Var1, ordre) > match(Var2, ordre)) %>%
  mutate(
    p_symbol = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    Var2 = factor(Var2, levels = ordre),
    Var1 = factor(Var1, levels = rev(ordre)),
    label_txt = if_else(
      p_symbol == "ns",
      sprintf("%.2f", r),
      paste0(sprintf("%.2f", r), p_symbol)
    )
  )

# 4. Heatmap
p_matrice_nephelo_33 <- ggplot(
  matrice_data,
  aes(x = Var2, y = Var1, fill = r)
) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = label_txt), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_matrice_nephelo_33)


```



# Figure 3
# Correlation plots comparing growth assessment methods

```{r}

# 1. List of the three X variables to compare with Dry_weight as Y
vars_x <- c("Diameter", "Spectro_AUC_33", "Nephelo_AUC_33")

# 2. Generate the three plots with Dry_weight on the y-axis
liste_graphiques_spearman <- map(vars_x, function(var_x) {
  
  # Check whether the variable is present in the dataframe
  if(!var_x %in% names(df_agreg_azote)) {
    warning(paste("La variable", var_x, "n'est pas présente dans df_agreg_azote."))
    return(NULL)
  }
  
  # Calculate the Spearman correlation test
  test_cor <- cor.test(df_agreg_azote[[var_x]], df_agreg_azote$Dry_weight, method = "spearman", exact = FALSE)
  r_val <- round(test_cor$estimate, 2)
  p_val <- test_cor$p.value
  p_text <- if_else(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
  
  # Construct the scatter plot (X = var_x, Y = Dry_weight)
  ggplot(df_agreg_azote, aes(x = .data[[var_x]], y = Dry_weight)) +
    geom_point(aes(color = Nitrogen_source), size = 3.5, alpha = 0.9) +
    scale_color_manual(values = custom_colors) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey50", linetype = "dashed", linewidth = 0.6) +
    labs(
      title = paste0(var_x),
      subtitle = paste0("Spearman r = ", r_val, " (", p_text, ")"),
      x = var_x,
      y = "Dry weight (mg)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5, face = "italic"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
})

# Remove any NULL elements
liste_graphiques_spearman <- compact(liste_graphiques_spearman)

# 3. Extract the common legend
p_legende <- ggplot(df_agreg_azote, aes(x = Diameter, y = Dry_weight, color = Nitrogen_source)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = custom_colors) +
  labs(color = "Nitrogen source") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8.5)
  )

legende_commune <- cowplot::get_legend(p_legende)

# 4. Assemble the final figure (three panels + legend)
grille_3_plots <- cowplot::plot_grid(plotlist = liste_graphiques_spearman, nrow = 1)
combinaison_finale_spearman <- cowplot::plot_grid(grille_3_plots, legende_commune, rel_widths = c(0.85, 0.15), ncol = 2)

# Display the figure
print(combinaison_finale_spearman)


```


# Figure 4
# Boxplots of standardized values for six nutrient conditions

```{r}

# 1. Order of nitrogen sources
ordre_azote <- c("Asparagine", "Ammonium", "Nitrate", "Glycine", 
                 "Tryptophan", "Nitrite", "Control")

# 2. Desired order of variables on the x-axis: dry weight, diameter, spectrophotometry, nephelometry
ordre_variables <- c("Dry Weight", "Diameter", "Spectro_AUC", "Nephelo_AUC")

# 3. Select and filter variables of interest
df_boxplot <- df_clean_wide %>%
  filter(Nitrogen_source %in% ordre_azote) %>%
  select(
    Nitrogen_source, 
    Dry_weight, 
    Diameter, 
    Spectro_AUC_33,
    Nephelo_AUC_33
  )

# 4. Standardization (Z-score)
df_standardise <- df_boxplot %>%
  mutate(across(where(is.numeric), ~ as.vector(scale(.)), .names = "std_{.col}")) %>%
  select(Nitrogen_source, starts_with("std_"))

# 5. Convert to long format, reorder factors, and remove NA values
df_long_auc <- df_standardise %>%
  pivot_longer(cols = starts_with("std_"), names_to = "Mesure", values_to = "Valeur") %>%
  filter(!is.na(Valeur) & !is.nan(Valeur)) %>%
  mutate(Mesure = case_when(
    Mesure == "std_Dry_weight"     ~ "Dry Weight",
    Mesure == "std_Diameter"       ~ "Diameter",
    Mesure == "std_Spectro_AUC_33" ~ "Spectro_AUC",
    Mesure == "std_Nephelo_AUC_33" ~ "Nephelo_AUC"
  )) %>%
  mutate(
    Nitrogen_source = factor(Nitrogen_source, levels = ordre_azote),
    Mesure          = factor(Mesure, levels = ordre_variables)
  )

# 6. Generate the plot
p_boxplot_auc <- ggplot(df_long_auc, aes(x = Mesure, y = Valeur, fill = Mesure)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5, color = "black", shape = 21, stroke = 0.3) +
  facet_wrap(~ Nitrogen_source, nrow = 1, drop = TRUE) + 
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Standardized Value (Z-score)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8.5), 
    strip.background = element_rect(fill = "grey93", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9.5),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2")

# Display the plot
print(p_boxplot_auc)

```


# Figure 5
# Boxplots of coefficients of variation

```{r}

# 1. Calculate the overall CV (%) by condition/nitrogen source
df_cv_global_cond <- df_clean_wide %>%
  group_by(Nitrogen_source) %>%
  summarise(
    `Dry Weight`  = (sd(Dry_weight, na.rm = TRUE) / mean(Dry_weight, na.rm = TRUE)) * 100,
    `Diameter`    = (sd(Diameter, na.rm = TRUE) / mean(Dry_weight, na.rm = TRUE)) * 100,
    `Spectro_AUC` = (sd(Spectro_AUC_33, na.rm = TRUE) / mean(Spectro_AUC_33, na.rm = TRUE)) * 100,
    `Nephelo_AUC` = (sd(Nephelo_AUC_33, na.rm = TRUE) / mean(Nephelo_AUC_33, na.rm = TRUE)) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(`Dry Weight`, `Diameter`, `Spectro_AUC`, `Nephelo_AUC`), 
    names_to = "Method", 
    values_to = "CV_percent"
  ) %>%
  filter(!is.na(CV_percent) & !is.nan(CV_percent) & !is.infinite(CV_percent)) %>%
  mutate(Method = factor(Method, levels = c("Dry Weight", "Diameter", "Spectro_AUC", "Nephelo_AUC")))

# 2. Tukey HSD post-hoc test for statistical significance letters
cv_model2   <- aov(CV_percent ~ Method, data = df_cv_global_cond)
cv_tukey2   <- TukeyHSD(cv_model2)
cv_letters2 <- multcompLetters4(cv_model2, cv_tukey2)$Method$Letters

df_letters2 <- data.frame(
  Method = names(cv_letters2),
  Letter = cv_letters2
) %>%
  left_join(df_cv_global_cond %>% group_by(Method) %>% summarise(Y_pos = max(CV_percent, na.rm = TRUE) + 3), by = "Method")

# 3. Plot
p_cv_global <- ggplot(df_cv_global_cond, aes(x = Method, y = CV_percent)) +
  geom_boxplot(fill = "grey95", color = "black", alpha = 0.6, outlier.shape = NA) + 
  geom_jitter(aes(color = Nitrogen_source), width = 0.2, size = 3, alpha = 0.9) +
  scale_color_manual(values = custom_colors) + 
  geom_text(data = df_letters2, aes(x = Method, y = Y_pos, label = Letter),
            vjust = 0, size = 5, fontface = "bold", inherit.aes = FALSE) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "",
    y = "Coefficient of Variation (CV %)",
    color = "Nitrogen Source"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, face = "bold", angle = 35, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_cv_global)


# Calculate mean CV (%) by method (with standard deviation)
df_cv_summary <- df_cv_global_cond %>%
  group_by(Method) %>%
  summarise(
    CV_moyen = mean(CV_percent, na.rm = TRUE),
    CV_sd    = sd(CV_percent, na.rm = TRUE),
    CV_min   = min(CV_percent, na.rm = TRUE),
    CV_max   = max(CV_percent, na.rm = TRUE),
    .groups  = "drop"
  )

print(as.data.frame(df_cv_summary))

```




# Figure 6 (data used for Figure 6)
# Spearman correlation matrices after 100 h of growth, with and without nitrite

```{r}

library(tidyverse)
library(Hmisc)

# Display order of variables
ordre <- c(
  "Dry weight", "AUC", "MCC", "Growth Rate",
  "Lag Time", "MCC Time", "Exponential duration"
)

# Function to generate the heatmap
creer_heatmap_100h <- function(df_input, titre_legende = "Spearman r") {
  
  # 1. Data preparation
  df_corr <- df_input %>%
    filter(
      (Appareil == "Nephelo" & Timing == "100") |
        (Appareil == "Boite" & Parameter == "Dry_weight")
    ) %>%
    mutate(
      Valeur = as.numeric(Valeur),
      Parameter = case_when(
        Parameter == "Nephelo-AUC" ~ "AUC",
        Parameter == "Nephelo-MCC" ~ "MCC",
        Parameter == "Nephelo-Growth_rate" ~ "Growth Rate",
        Parameter == "Nephelo-Lag_time" ~ "Lag Time",
        Parameter == "Nephelo-Plateau_time" ~ "MCC Time",
        Parameter == "Nephelo-Exponential_duration" ~ "Exponential duration",
        Parameter == "Dry_weight" ~ "Dry weight"
      )
    ) %>%
    group_by(Condition, Parameter) %>%
    summarise(Valeur = mean(Valeur, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Parameter, values_from = Valeur) %>%
    select(any_of(ordre))
  
  # 2. Spearman correlations
  ct <- rcorr(as.matrix(df_corr), type = "spearman")
  
  # 3. Format data for the heatmap
  matrice_data <- as.data.frame(as.table(ct$r)) %>%
    rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
    mutate(p = as.vector(ct$P)) %>%
    filter(match(Var1, ordre) > match(Var2, ordre)) %>%
    mutate(
      p_symbol = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      Var2 = factor(Var2, levels = ordre),
      Var1 = factor(Var1, levels = rev(ordre)),
      label_txt = if_else(
        p_symbol == "ns",
        sprintf("%.2f", r),
        paste0(sprintf("%.2f", r), p_symbol)
      )
    )
  
  # 4. Generate ggplot heatmap
  ggplot(matrice_data, aes(x = Var2, y = Var1, fill = r)) +
    geom_tile(color = "black", linewidth = 0.5) +
    geom_text(aes(label = label_txt), size = 4, fontface = "bold") +
    scale_fill_gradient2(
      low = "#4575b4",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      limits = c(-1, 1),
      name = titre_legende
    ) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, face = "bold", color = "black"),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
}



# 1. Correlation matrix after 100 h (all conditions)
p_matrice_nephelo_100 <- creer_heatmap_100h(data)

# 2. Correlation matrix after 100 h (excluding nitrite)
p_matrice_nephelo_100_sans_nitrite <- creer_heatmap_100h(
  data %>% filter(Condition != "Nitrite")
)

# Display the two figures
print(p_matrice_nephelo_100)
print(p_matrice_nephelo_100_sans_nitrite)

```

# Figure 7
# Glucose and asparagine concentration matrix

```{r}

donnees_matrice <- read.csv2("")

ordre_glucose    <- c("Glucose20", "Glucose10", "Glucose5", "Glucose1")
ordre_asparagine <- c("Asparagine5", "Asparagine4", "Asparagine3", "Asparagine2", "Asparagine1")

df_samples <- donnees_matrice %>%
  filter(Condition_C != "Glucose0") %>%
  pivot_longer(cols = matches("^[0-9]+ h"), names_to = "Time_Raw", values_to = "Value_Raw") %>%
  mutate(Value = as.numeric(str_replace_all(Value_Raw, ",", "."))) %>%
  mutate(
    heures = as.numeric(str_extract(Time_Raw, "^[0-9]+")),
    minutes = str_extract(Time_Raw, "[0-9]+(?=\\s*min)"),
    minutes = if_else(is.na(minutes), 0, as.numeric(minutes)),
    Time_Hours = heures + (minutes / 60)
  ) %>%
  mutate(
    Condition_C = factor(Condition_C, levels = ordre_glucose),
    Condition_N = factor(Condition_N, levels = ordre_asparagine)
  ) %>%
  filter(Time_Hours >= 0 & Time_Hours <= 100)

df_stats_facet <- df_samples %>%
  group_by(Condition_C, Condition_N, Conca, Time_Hours) %>%
  summarise(
    Moyenne_Instantanee = mean(Value, na.rm = TRUE),
    SD_Instantane = sd(Value, na.rm = TRUE),
    CV_Instantane = if_else(Moyenne_Instantanee > 10, (SD_Instantane / Moyenne_Instantanee) * 100, 0),
    .groups = "drop_last"
  ) %>%
  summarise(
    Moyenne_CV_Temporel = mean(CV_Instantane, na.rm = TRUE),
    N_Points = n(),
    Erreur_Type = sd(CV_Instantane, na.rm = TRUE) / sqrt(N_Points),
    Marge_IC = qt(0.975, df = N_Points - 1) * Erreur_Type,
    .groups = "drop"
  ) %>%
  mutate(
    Label_CV = paste0("Avg CV = ", round(Moyenne_CV_Temporel, 2), "% ± ", round(Marge_IC, 2), "%")
  )

min_global <- min(df_samples$Value, na.rm = TRUE)
df_stats_facet <- df_stats_facet %>% mutate(X_pos = 5, Y_pos = min_global * 1.05)

df_summary_matrix <- df_samples %>%
  group_by(Condition_C, Condition_N, Conca, Time_Hours) %>%
  summarise(
    Moyenne = mean(Value, na.rm = TRUE),
    Ecart_Type = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(
  ggplot(df_summary_matrix, aes(x = Time_Hours, y = Moyenne, color = Conca, fill = Conca)) +
    geom_ribbon(aes(ymin = Moyenne - Ecart_Type, ymax = Moyenne + Ecart_Type), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_text(data = df_stats_facet, aes(x = X_pos, y = Y_pos, label = Label_CV),
              hjust = 0, vjust = 0, size = 2.8, color = "black", fontface = "bold", inherit.aes = FALSE) +
    coord_cartesian(xlim = c(0, 100)) +
    facet_grid(Condition_C ~ Condition_N) +
    labs(
      title = "Kinetic Profile Variability Standardized by Growth Peak",
      subtitle = "Temporal average CV (%) computed point-by-point (mean ± 95% confidence interval)",
      x = "Time (hours)",
      y = "Growth Signal Output (Nephelometry)"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey92"),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
)

```


# Figure 8
# Spearman correlation matrix between nephelometric growth parameters after nutrient optimization

```{r}

library(tidyverse)
library(Hmisc)

# 1. Data preparation
df_corr <- data %>%
  filter(
    (Appareil == "Nephelo" & Timing == "Opti") |
      (Appareil == "Boite" & Parameter == "Dry_weight")
  ) %>%
  filter(!Condition %in% c("Arginine", "Serine")) %>%
  mutate(
    Valeur = as.numeric(str_replace_all(Valeur, ",", ".")),
    Parameter = case_when(
      Parameter == "Nephelo-AUC" ~ "AUC",
      Parameter == "Nephelo-MCC" ~ "MCC",
      Parameter == "Nephelo-Growth_rate" ~ "Growth Rate",
      Parameter == "Nephelo-Lag_time" ~ "Lag Time",
      Parameter == "Nephelo-Plateau_time" ~ "MCC Time",
      Parameter == "Nephelo-Exponential_duration" ~ "Exponential duration",
      Parameter == "Dry_weight" ~ "Dry weight",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Parameter)) %>%
  
  # Average technical replicates
  group_by(Condition, Rep, Parameter) %>%
  summarise(Valeur = mean(Valeur, na.rm = TRUE), .groups = "drop") %>%
  
  # Convert to wide format
  pivot_wider(names_from = Parameter, values_from = Valeur) %>%
  
  # Average biological replicates
  group_by(Condition) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  
  # Select variables in the desired order
  select(any_of(c(
    "Dry weight", "AUC", "MCC", "Growth Rate",
    "Lag Time", "MCC Time", "Exponential duration"
  )))

# 2. Spearman correlations
ct <- rcorr(as.matrix(df_corr), type = "spearman")

# 3. Heatmap data
ordre <- c(
  "Dry weight", "AUC", "MCC", "Growth Rate",
  "Lag Time", "MCC Time", "Exponential duration"
)

matrice_data <- as.data.frame(as.table(ct$r)) %>%
  rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
  mutate(p = as.vector(ct$P)) %>%
  filter(match(Var1, ordre) > match(Var2, ordre)) %>%
  mutate(
    p_symbol = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    Var2 = factor(Var2, levels = ordre),
    Var1 = factor(Var1, levels = rev(ordre)),
    label_txt = if_else(
      p_symbol == "ns",
      sprintf("%.2f", r),
      paste0(sprintf("%.2f", r), p_symbol)
    )
  )

# 4. Heatmap
p_matrice_nephelo_opti <- ggplot(
  matrice_data,
  aes(x = Var2, y = Var1, fill = r)
) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = label_txt), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_matrice_nephelo_opti)

```

# Supplementary Figures 2A, 2B, and 4
# Growth curves

```{r}

library(readxl)
library(tidyverse)

nephelo <- read_excel("Comparaison nephelo spectro/output_result_nephelo_33h/data_ntu_clean.xlsx")
spectro <- read_excel("Comparaison nephelo spectro/output_result_spectro_33h/data_ntu_clean.xlsx")
nephelo_opti <- read_excel("Optimisation nephelo/output_result/data_ntu_clean.xlsx")

custom_colors <- c(
  "Ammonium"="#C3223B", "Asparagine"="#EF6D6D", "Control"="#000000", "GABA"="#E5E87B",
  "Glutamate"="#0BA279", "Glycine"="#69955C", "Nitrate"="#51FF00", "Nitrite"="#73A6F0",
  "Phenylalanine"="#3B8ED2", "Proline"="#5200FF", "Serine"="#9B8BDE", "Taurine"="#CEB551",
  "Tryptophan"="#808D85"
)

plot_growth <- function(df, titre, unit="kNTU", div=1000, line33=FALSE) {
  
  # 1. Convert to long format and clean data
  d_long <- df %>%
    mutate(Condition = recode(Condition, "Tryptophane" = "Tryptophan")) %>%
    pivot_longer(starts_with("min_"), names_to = "Time", values_to = "Signal") %>%
    mutate(
      Time_h = as.numeric(str_remove(Time, "min_")) / 60,
      Signal = Signal / div
    ) %>%
    filter(Time_h <= 100)
  
  # 2. Average technical replicates within each biological replicate
  d_rep <- d_long %>%
    group_by(Condition, Rep, Time_h) %>%
    summarise(Signal_rep = mean(Signal, na.rm = TRUE), .groups = "drop")
  
  # 3. Order facets according to maximum signal
  ord <- d_rep %>%
    group_by(Condition) %>%
    summarise(M = max(Signal_rep, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(M)) %>%
    pull(Condition)
  
  # 4. Calculate overall mean and SE across biological replicates
  d_final <- d_rep %>%
    mutate(Condition = factor(Condition, levels = ord)) %>%
    group_by(Condition, Time_h) %>%
    summarise(
      Mean = mean(Signal_rep, na.rm = TRUE),
      SE   = sd(Signal_rep, na.rm = TRUE) / sqrt(sum(!is.na(Signal_rep))),
      .groups = "drop"
    )
  
  # 5. Generate the plot
  p <- ggplot(d_final, aes(Time_h, Mean)) +
    geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Condition), alpha = 0.14, color = NA) +
    geom_line(aes(color = Condition), linewidth = 1.15)
  
  if (line33) {
    p <- p + geom_vline(xintercept = 33, linetype = "dashed", linewidth = 0.6, color = "grey40")
  }
  
  p <- p + facet_wrap(~Condition, nrow = 3, ncol = 4) +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max(d_final$Mean + d_final$SE, na.rm = TRUE)), expand = expansion(mult = c(0, 0.03))) +
    labs(title = titre, x = "Time (h)", y = unit) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      strip.background = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 8.5, color = "black"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0.9, "lines")
    )
  
  return(p)
}


p_nephelo <- plot_growth(nephelo, "Nephelometry", unit = "kNTU", div = 1000, line33 = TRUE)
p_spectro <- plot_growth(spectro, "Spectrophotometry", unit = "OD", div = 1, line33 = TRUE)
p_opti    <- plot_growth(nephelo_opti, "Optimized nephelometry", unit = "kNTU", div = 1000, line33 = FALSE)

print(p_nephelo)
print(p_spectro)
print(p_opti)

```

# Supplementary Figure 3
# Spearman correlation matrix between spectrophotometric parameters after 33 h of growth and dry weight

```{r}

library(tidyverse)
library(Hmisc)

# 1. Data preparation
df_corr <- data %>%
  filter(
    (Appareil == "Spectro" & Timing == "33") |
      (Appareil == "Boite" & Parameter %in% c("Dry_weight", "Diameter"))
  ) %>%
  filter(!Condition %in% c("Arginine", "Serine")) %>%
  mutate(
    Valeur = as.numeric(str_replace_all(Valeur, ",", ".")),
    Parameter = case_when(
      Parameter == "Spectro-AUC" ~ "AUC",
      Parameter == "Spectro-MCC" ~ "MCC",
      Parameter == "Spectro-Growth_rate" ~ "Growth Rate",
      Parameter == "Spectro-Lag_time" ~ "Lag Time",
      Parameter == "Spectro-Plateau_time" ~ "MCC Time",
      Parameter == "Spectro-Exponential_duration" ~ "Exponential duration",
      Parameter == "Dry_weight" ~ "Dry weight",
      Parameter == "Diameter" ~ "Diameter",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Parameter)) %>%
  
  # Average technical replicates
  group_by(Condition, Rep, Parameter) %>%
  summarise(Valeur = mean(Valeur, na.rm = TRUE), .groups = "drop") %>%
  
  # Convert to wide format
  pivot_wider(names_from = Parameter, values_from = Valeur) %>%
  
  # Average biological replicates
  group_by(Condition) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  
  # Select variables in the desired order
  select(any_of(c(
    "Dry weight", "Diameter", "AUC", "MCC", "Growth Rate",
    "Lag Time", "MCC Time", "Exponential duration"
  )))

# 2. Spearman correlations
ct <- rcorr(as.matrix(df_corr), type = "spearman")

# 3. Heatmap data
ordre <- c(
  "Dry weight", "Diameter", "AUC", "MCC", "Growth Rate",
  "Lag Time", "MCC Time", "Exponential duration"
)

matrice_data <- as.data.frame(as.table(ct$r)) %>%
  rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
  mutate(p = as.vector(ct$P)) %>%
  filter(match(Var1, ordre) > match(Var2, ordre)) %>%
  mutate(
    p_symbol = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    Var2 = factor(Var2, levels = ordre),
    Var1 = factor(Var1, levels = rev(ordre)),
    label_txt = if_else(
      p_symbol == "ns",
      sprintf("%.2f", r),
      paste0(sprintf("%.2f", r), p_symbol)
    )
  )

# 4. Heatmap
p_matrice_spectro_33 <- ggplot(
  matrice_data,
  aes(x = Var2, y = Var1, fill = r)
) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = label_txt), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman r"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_matrice_spectro_33)

```

# Supplementary Figure 5
# Correlation matrix with scatterplots for AUC, MCC, growth rate, and lag time after 100 h of growth under optimized nutrient conditions


```{r}

library(tidyverse)
library(GGally)

# Dynamically select optimized variables
col_target <- intersect(
  c("Dry_weight",
    "Nephelo_MCC_Opti",
    "Nephelo_AUC_Opti",
    "Nephelo_Growth_rate_Opti",
    "Nephelo_Lag_time_Opti"),
  names(df_agreg_azote)
)

df_matrix_plot <- df_agreg_azote %>% select(Nitrogen_source, all_of(col_target))

# Function for scatterplots with regression line and Spearman correlation
points_et_cor_s <- function(data, mapping, ...) {
  x_var <- rlang::as_name(mapping$x)
  y_var <- rlang::as_name(mapping$y)
  
  test_cor <- cor.test(data[[x_var]], data[[y_var]], method = "spearman", exact = FALSE)
  r_val <- round(test_cor$estimate, 2)
  p_val <- test_cor$p.value
  p_text <- if_else(p_val < 0.001, "p < 0.001", paste0("p = ", round(p_val, 3)))
  label_txt <- paste0("r = ", r_val, "\n", p_text)
  
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 3, alpha = 0.8) +
    # Add global regression line (not segmented by condition)
    geom_smooth(
      method = "lm", 
      formula = y ~ x, 
      se = FALSE, 
      color = "grey40", 
      linetype = "dashed", 
      linewidth = 0.6,
      inherit.aes = FALSE,
      aes(x = .data[[x_var]], y = .data[[y_var]])
    ) +
    annotate("label", x = -Inf, y = Inf, label = label_txt, 
             hjust = -0.1, vjust = 1.1, size = 3, 
             label.size = 0, fill = ggplot2::alpha("white", 0.6)) +
    theme_bw()
}

# Generate the ggpairs correlation matrix
p_matrix_s <- ggpairs(
  df_matrix_plot, 
  columns = 2:ncol(df_matrix_plot),
  mapping = aes(color = Nitrogen_source),
  lower = list(continuous = points_et_cor_s),
  upper = list(continuous = "blank"), 
  diag = list(continuous = wrap("densityDiag", alpha = 0.8)),
  legend = c(1, 1) 
) + 
  theme_bw() +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(), 
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  labs(color = "Nitrogen Source")

# Clean up the upper triangle
for(i in 1:p_matrix_s$nrow) {
  for(j in 1:p_matrix_s$ncol) {
    if(i < j) { 
      p_matrix_s[i, j] <- p_matrix_s[i, j] + 
        theme(panel.border = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
    }
  }
}

print(p_matrix_s)

```



