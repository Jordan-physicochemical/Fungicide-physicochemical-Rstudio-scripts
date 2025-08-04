Figure 2 (3 plots: plot1 uptake, plot2 mobility, plot3 direction of translocation)
# === LOAD LIBRARIES === plot3
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(stringr)
library(car)        # For Levene's test
library(dunn.test)  # For Dunn's post-hoc test

# === PARAMETERS ===
fungicides_to_plot <- c("AZX", "CYP", "FLU", "MFN")
included_tissues   <- c("T", "B")  # Only T and B
treatment_mapping <- list(
  FLU = c("Mix", "Impact", "FLU"),
  CYP = c("Mix", "Amistar", "CYP"),
  AZX = c("Mix", "Amistar", "AZX"),
  MFN = c("Mix", "Belanty", "MFN")
)
fungicide_order <- c("FLU", "CYP", "AZX", "MFN")
tissues    <- c("T", "B")
fungicides <- fungicides_to_plot

# === LOAD DATA ===
data_raw <- read_csv("FORM5 processed data.csv", skip=1)
colnames(data_raw) <- str_trim(colnames(data_raw))
data <- data_raw %>% filter(DAT == 7)

# === Generate both unmodified and modified versions ===

# Pivot to long format
long_data_orig <- data %>%
  filter(Tissue %in% included_tissues) %>%
  pivot_longer(
    cols = all_of(fungicides_to_plot),
    names_to = "Fungicide",
    values_to = "Value"
  ) %>%
  mutate(Value = as.numeric(Value))

# Make a *copy* for modification and add a flag
long_data_mod <- long_data_orig %>% mutate(Replaced_B_zero = FALSE)

for (fungicide in fungicides_to_plot) {
  fungicide_data <- long_data_mod %>%
    filter(Fungicide == fungicide)
  wide_data <- fungicide_data %>%
    select(Repetition, Treatment, Tissue, Value) %>%
    pivot_wider(names_from = Tissue, values_from = Value, names_prefix = "Tissue_")
  to_replace <- wide_data %>%
    filter(Tissue_B == 0 & Tissue_T > 0)
  if (nrow(to_replace) > 0) {
    for (i in 1:nrow(to_replace)) {
      idx <- which(
        long_data_mod$Fungicide == fungicide &
          long_data_mod$Repetition == to_replace$Repetition[i] &
          long_data_mod$Treatment == to_replace$Treatment[i] &
          long_data_mod$Tissue == "B" &
          long_data_mod$Value == 0
      )
      if(length(idx) == 1) {
        long_data_mod$Value[idx] <- 125
        long_data_mod$Replaced_B_zero[idx] <- TRUE
      }
    }
  }
}

# === Shared function to calculate proportions ===
append_props <- function(df) {
  df %>%
    group_by(Repetition, Fungicide, Treatment) %>%
    mutate(Total_TB = sum(Value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Proportion = (Value / Total_TB) * 100)
}
long_data_orig <- append_props(long_data_orig)
long_data_mod  <- append_props(long_data_mod)

# === Shared filter and factor setup ===
apply_filters <- function(df) {
  df %>%
    filter(
      (Fungicide == "FLU" & Treatment %in% treatment_mapping$FLU) |
        (Fungicide == "CYP" & Treatment %in% treatment_mapping$CYP) |
        (Fungicide == "AZX" & Treatment %in% treatment_mapping$AZX) |
        (Fungicide == "MFN" & Treatment %in% treatment_mapping$MFN)
    ) %>%
    mutate(
      Treatment = case_when(
        Fungicide == "FLU" ~ factor(Treatment, levels = treatment_mapping$FLU),
        Fungicide == "CYP" ~ factor(Treatment, levels = treatment_mapping$CYP),
        Fungicide == "AZX" ~ factor(Treatment, levels = treatment_mapping$AZX),
        Fungicide == "MFN" ~ factor(Treatment, levels = treatment_mapping$MFN),
        TRUE ~ as.factor(Treatment)
      ),
      Tissue    = factor(Tissue, levels = included_tissues),
      Fungicide = factor(Fungicide, levels = fungicide_order)
    )
}
long_data_orig <- apply_filters(long_data_orig)
long_data_mod  <- apply_filters(long_data_mod)

# === XY label data (use modified for plotting) ===
xy_label_data <- long_data_mod %>%
  group_by(Fungicide, Treatment, Tissue) %>%
  summarise(
    n_total = sum(!is.na(Value)),               # Total reps (including 0s)
    n_nonzero = sum(Value > 0, na.rm = TRUE),   # Reps with value > 0
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("(", n_nonzero, "/", n_total, ")"),
    y_pos = 1
  ) %>%
  mutate(
    Treatment = case_when(
      Fungicide == "FLU" ~ factor(Treatment, levels = treatment_mapping$FLU),
      Fungicide == "CYP" ~ factor(Treatment, levels = treatment_mapping$CYP),
      Fungicide == "AZX" ~ factor(Treatment, levels = treatment_mapping$AZX),
      Fungicide == "MFN" ~ factor(Treatment, levels = treatment_mapping$MFN),
      TRUE ~ as.factor(Treatment)
    ),
    Fungicide = factor(Fungicide, levels = fungicide_order),
    Tissue = factor(Tissue, levels = included_tissues)
  )

# === Plotting (use long_data_mod for highlighting replaced points) ===
plot3 <- ggplot(long_data_mod %>% filter(Value > 0), 
                aes(x = Tissue, y = Proportion, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  # non-replaced points
  geom_jitter(
    data = long_data_mod %>% filter(Value > 0, !Replaced_B_zero),
    aes(color = Treatment),
    size = 1.2, alpha = 0.6,
    position = position_dodge(width = 0.8)
  ) +
  # replaced points in bold red
  geom_jitter(
    data = long_data_mod %>% filter(Value > 0, Replaced_B_zero),
    color = "red", shape = 16, size = 4, alpha = 0.85,
    position = position_dodge(width = 0.8)
  ) +
  geom_text(data = xy_label_data, 
            aes(x = Tissue, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.2, fontface = "plain",
            position = position_dodge(width = 0.8)) +
  facet_wrap(~ Fungicide, nrow = 4, ncol = 1) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = "Tissue",
    y = "Proportion (%) of Uptake in Tissue (log scale)",
    title = "Tissue Distribution of Treatments for Each Fungicide (7 DAT, Log Scale)\nRed points: imputed B=0 set to 125"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
plot3

# === ANALYSIS FUNCTIONS RUN ON BOTH DATASETS ===

stats_suite <- function(df, label_for_output = "UNMODIFIED") {
  cat("\n\n=====================", label_for_output, "DATA =====================\n")
  for(fungicide in fungicides) {
    for(tissue in tissues) {
      data_subset <- df %>%
        filter(Fungicide == fungicide, Tissue == tissue, !is.na(Proportion))
      if (length(unique(data_subset$Treatment)) < 2) {
        cat(sprintf("%s %s: Not enough treatment groups for analysis\n", fungicide, tissue))
        next
      }
      cat(sprintf("\n--- %s Fungicide, %s Tissue ---\n", fungicide, tissue))
      
      aov_mod <- aov(Proportion ~ Treatment, data = data_subset)
      shapiro_test <- shapiro.test(residuals(aov_mod))
      cat(sprintf("Shapiro-Wilk Test: W = %.4f, p = %.4f, n = %d\n",
                  shapiro_test$statistic, shapiro_test$p.value, length(residuals(aov_mod))))
      levene_test <- car::leveneTest(Proportion ~ Treatment, data = data_subset)
      cat(sprintf("Levene's Test: F(%d,%d) = %.4f, p = %.4f\n",
                  levene_test$Df[1], levene_test$Df[2], levene_test$`F value`[1], levene_test$`Pr(>F)`[1]))
      
      # Kruskal-Wallis and Dunn's
      kw_test <- kruskal.test(Proportion ~ Treatment, data = data_subset)
      cat(sprintf("Kruskal-Wallis: H = %.3f, df = %d, p = %.4f\n",
                  kw_test$statistic, kw_test$parameter, kw_test$p.value))
      if (kw_test$p.value < 0.05) {
        dunn_result <- dunn.test(data_subset$Proportion, data_subset$Treatment, 
                                 method = "bonferroni", alpha = 0.05, list = FALSE)
        sig_comparisons <- dunn_result$comparisons[dunn_result$P.adjusted < 0.05]
        if(length(sig_comparisons) == 0) {
          cat("(No significant Dunn's pairwise comparisons)\n")
        } else {
          cat("Significant Dunn's pairwise differences:\n")
          print(data.frame(
            Comparison = dunn_result$comparisons,
            z = dunn_result$Z,
            p_unadj = dunn_result$P,
            p_adj = dunn_result$P.adjusted,
            sig = dunn_result$P.adjusted < 0.05
          )[dunn_result$P.adjusted < 0.05, ])
        }
      }
    }
  }
}

# --- RUN ANALYSIS ---
stats_suite(long_data_orig, "UNMODIFIED (zeroes)")
stats_suite(long_data_mod,  "MODIFIED (zeroes set to 125)")

# === DESCRIPTIVE STATISTICS ===
cat("\n=== DESCRIPTIVE STATISTICS (MODIFIED DATA) ===\n")
descriptive_stats <- long_data_mod %>%
  group_by(Fungicide, Tissue, Treatment) %>%
  summarise(
    mean_proportion = mean(Proportion, na.rm = TRUE),
    sd_proportion = sd(Proportion, na.rm = TRUE),
    median_proportion = median(Proportion, na.rm = TRUE),
    n = sum(!is.na(Proportion)),
    .groups = "drop"
  ) %>%
  arrange(Fungicide, Tissue, Treatment)
print(descriptive_stats)



####plot1
library(ggplot2)
library(dplyr)
library(car)        # for Levene's test
library(dunn.test)  # for Dunn's test
library(multcompView)  # for compact letter display

# Define your custom treatment levels for each fungicide
custom_treatment_mapping <- list(
  FLU = c("Mix", "Impact", "FLU"),
  CYP = c("Mix", "Amistar", "CYP"),
  AZX = c("Mix", "Amistar", "AZX"),
  MFN = c("Mix", "Belanty", "MFN")
)

# Calculate total uptake over tissues T, M, B for each fungicide-treatment-repetition
# Filter out non-numeric values and missing data
sum_uptake <- long_data %>%
  filter(Tissue %in% c("T", "M", "B")) %>%
  # Convert Value to numeric, which will turn non-numeric values like "removed" to NA
  mutate(Value = as.numeric(as.character(Value))) %>%
  # Group and calculate total, but only for complete cases
  group_by(Fungicide, Treatment, Repetition) %>%
  # Only keep groups where all tissue values are numeric (not NA)
  filter(all(!is.na(Value))) %>%
  summarise(Total_Uptake = sum(Value, na.rm=TRUE), .groups = "drop") %>%
  # Additional filter to remove any remaining non-numeric or negative values
  filter(is.finite(Total_Uptake) & Total_Uptake >= 0) %>%
  mutate(
    Treatment = case_when(
      Fungicide == "FLU" ~ factor(Treatment, levels = custom_treatment_mapping$FLU),
      Fungicide == "CYP" ~ factor(Treatment, levels = custom_treatment_mapping$CYP),
      Fungicide == "AZX" ~ factor(Treatment, levels = custom_treatment_mapping$AZX),
      Fungicide == "MFN" ~ factor(Treatment, levels = custom_treatment_mapping$MFN),
      TRUE ~ as.factor(Treatment)
    ),
    Fungicide = factor(Fungicide, levels = c("FLU", "CYP", "AZX", "MFN"))
  )

# Report data exclusions
cat("\n================== DATA EXCLUSION REPORT ==================\n")
excluded_data <- long_data %>%
  filter(Tissue %in% c("T", "M", "B")) %>%
  mutate(Value_numeric = as.numeric(as.character(Value))) %>%
  group_by(Fungicide, Treatment, Repetition) %>%
  summarise(
    has_non_numeric = any(is.na(Value_numeric) & !is.na(Value)),
    has_missing = any(is.na(Value)),
    .groups = "drop"
  ) %>%
  filter(has_non_numeric | has_missing)

if (nrow(excluded_data) > 0) {
  cat("Excluded repetitions due to non-numeric or missing values:\n")
  print(excluded_data)
} else {
  cat("No repetitions excluded - all data is numeric and complete.\n")
}

# Report final sample sizes
sample_sizes <- sum_uptake %>%
  group_by(Fungicide, Treatment) %>%
  summarise(n_reps = n(), .groups = "drop")

cat("\nFinal sample sizes after exclusions:\n")
print(sample_sizes)

# === REPORT MEAN UPTAKE VALUES ===
cat("\n================== MEAN UPTAKE VALUES ==================\n")

mean_uptake_values <- sum_uptake %>%
  group_by(Fungicide, Treatment) %>%
  summarise(
    mean_uptake = mean(Total_Uptake, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  arrange(Fungicide, Treatment)

for(fungicide in c("FLU", "CYP", "AZX", "MFN")) {
  cat(sprintf("\n%s:\n", fungicide))
  fungicide_data <- mean_uptake_values %>% filter(Fungicide == fungicide)
  for(i in 1:nrow(fungicide_data)) {
    cat(sprintf("  %s: %.2f (n=%d)\n", 
                fungicide_data$Treatment[i], 
                fungicide_data$mean_uptake[i],
                fungicide_data$n_reps[i]))
  }
}

# Function to create proper p-value matrix for multcompLetters
create_pvalue_matrix <- function(treatments, comparisons, p_values) {
  n <- length(treatments)
  pmat <- matrix(1, nrow = n, ncol = n)
  rownames(pmat) <- treatments
  colnames(pmat) <- treatments
  
  # Fill in p-values from comparisons
  for (i in seq_along(comparisons)) {
    pair <- strsplit(comparisons[i], " - ")[[1]]
    if (length(pair) == 2 && all(pair %in% treatments)) {
      pmat[pair[1], pair[2]] <- p_values[i]
      pmat[pair[2], pair[1]] <- p_values[i]
    }
  }
  
  # Set diagonal to 1 (no difference between same groups)
  diag(pmat) <- 1
  
  return(pmat)
}

# Function to run tests and return CLD letters for plotting
analyze_fungicide <- function(df, fungicide_name) {
  cat("\n============================================\n")
  cat(sprintf("Analysis for Fungicide: %s\n", fungicide_name))
  cat("============================================\n")
  
  df_sub <- df %>% filter(Fungicide == fungicide_name)
  n_samples <- nrow(df_sub)
  
  if(length(unique(df_sub$Treatment)) < 2) {
    cat("Not enough treatment groups for statistical testing.\n")
    return(NULL)
  }
  
  # Check if any treatment has too few samples
  sample_check <- df_sub %>%
    group_by(Treatment) %>%
    summarise(n = n(), .groups = "drop")
  
  cat("Sample sizes per treatment:\n")
  print(sample_check)
  
  if(any(sample_check$n < 2)) {
    cat("Warning: Some treatments have fewer than 2 samples.\n")
  }
  
  # Shapiro-Wilk test on residuals of ANOVA
  aov_model <- aov(Total_Uptake ~ Treatment, data = df_sub)
  shapiro_res <- shapiro.test(residuals(aov_model))
  cat(sprintf("Shapiro-Wilk test for residuals: W = %.4f, p = %.4f, n = %d\n",
              shapiro_res$statistic, shapiro_res$p.value, n_samples))
  
  # Levene's test for homogeneity of variances
  levene_res <- car::leveneTest(Total_Uptake ~ Treatment, data = df_sub)
  cat(sprintf("Levene's test: F(%d, %d) = %.4f, p = %.4f\n",
              levene_res$Df[1], levene_res$Df[2], levene_res$`F value`[1], levene_res$`Pr(>F)`[1]))
  
  # Kruskal-Wallis test
  kw_res <- kruskal.test(Total_Uptake ~ Treatment, data = df_sub)
  cat(sprintf("Kruskal-Wallis: H = %.3f, df = %d, p = %.4f\n",
              kw_res$statistic, kw_res$parameter, kw_res$p.value))
  
  cld <- NULL
  treatments <- levels(df_sub$Treatment)
  
  if(kw_res$p.value < 0.05) {
    # Run Dunn's test
    dunn_res <- dunn.test(df_sub$Total_Uptake, df_sub$Treatment, 
                          method = "bonferroni", alpha = 0.05, list = FALSE)
    
    # Print significant comparisons
    sig_comparisons <- data.frame(
      Comparison = dunn_res$comparisons,
      Z_score = dunn_res$Z,
      p_unadjusted = dunn_res$P,
      p_adjusted = dunn_res$P.adjusted,
      Significant = dunn_res$P.adjusted < 0.05
    )
    
    sig_pairs <- sig_comparisons[sig_comparisons$Significant, ]
    if(nrow(sig_pairs) == 0) {
      cat("No significant pairwise differences found.\n")
      # Assign same letter to all groups
      cld <- setNames(rep("a", length(treatments)), treatments)
    } else {
      cat("Significant pairwise comparisons:\n")
      print(sig_pairs)
      
      # Create p-value matrix for multcompLetters
      pmat <- create_pvalue_matrix(treatments, dunn_res$comparisons, dunn_res$P.adjusted)
      
      # Generate compact letter display
      cld_result <- multcompView::multcompLetters(pmat, threshold = 0.05)
      cld <- cld_result$Letters
    }
  } else {
    cat("Kruskal-Wallis test not significant - all groups get same letter.\n")
    # Assign same letter to all groups
    cld <- setNames(rep("a", length(treatments)), treatments)
  }
  
  return(cld)
}

# Run analysis over all fungicides, collect letters
cld_list <- list()
for (fungicide in levels(sum_uptake$Fungicide)) {
  cld_list[[fungicide]] <- analyze_fungicide(sum_uptake, fungicide)
}

# Prepare data frame with letters for plotting
letter_labels <- do.call(rbind, lapply(names(cld_list), function(fg) {
  if(is.null(cld_list[[fg]])) return(NULL)
  data.frame(
    Fungicide = fg,
    Treatment = names(cld_list[[fg]]),
    Letters = cld_list[[fg]],
    stringsAsFactors = FALSE
  )
}))

# Ensure Treatment factor levels are correctly ordered for labels
if (!is.null(letter_labels) && nrow(letter_labels) > 0) {
  letter_labels <- letter_labels %>%
    rowwise() %>%
    mutate(Treatment = factor(Treatment, levels = custom_treatment_mapping[[Fungicide]])) %>%
    ungroup() %>%
    mutate(Fungicide = factor(Fungicide, levels = c("FLU", "CYP", "AZX", "MFN")))
  
  # Get max uptake per fungicide-treatment group to position letters just above boxes
  max_uptake <- sum_uptake %>%
    group_by(Fungicide, Treatment) %>%
    summarise(max_uptake = max(Total_Uptake, na.rm = TRUE), .groups = "drop")
  
  # Join max uptake data to letter labels
  letter_labels <- letter_labels %>%
    left_join(max_uptake, by = c("Fungicide", "Treatment")) %>%
    mutate(y_pos = max_uptake * 1.05)  # Position letters slightly above max
  
  cat("\nSignificance letters for plotting:\n")
  print(letter_labels)
} else {
  cat("\nNo significance letters generated.\n")
  letter_labels <- data.frame()
}

# Plot with significance letters
plot1 <- ggplot(sum_uptake, aes(x = Treatment, y = Total_Uptake, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6, aes(color = Treatment)) +
  facet_wrap(~ Fungicide, nrow = 4, ncol = 1, scales = "free") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = "Treatment",
    y = "Total Uptake (T + M + B)",
    title = "Total Fungicide Uptake per Treatment (7 DAT)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  # Add top margin for letters
  )

# Add significance letters if they exist
if (nrow(letter_labels) > 0) {
  plot1 <- plot1 + 
    geom_text(data = letter_labels, 
              aes(x = Treatment, y = y_pos, label = Letters), 
              inherit.aes = FALSE, 
              size = 4, 
              vjust = 0, 
              fontface = "bold") +
    coord_cartesian(clip = "off")  # So letters are visible above plot area
}

print(plot1)



# === LOAD LIBRARIES === plot2
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(car)           # For Levene's test
library(dunn.test)     # For Dunn's post-hoc test
library(multcompView)  # For compact letter display (significance letters)

# === LOAD DATA ===
data <- read_csv("FORM5 processed data.csv", skip = 1)
colnames(data) <- str_trim(colnames(data))

# === FILTER FOR 7 DAT ONLY ===
data <- data %>% filter(DAT == 7)

# === PARAMETERS ===
fungicides_to_plot <- c("AZX", "CYP", "FLU", "MFN")
included_tissues <- c("T", "M", "B")
treatment_mapping <- list(
  FLU = c("Mix", "Impact", "FLU"),
  CYP = c("Mix", "Amistar", "CYP"),
  AZX = c("Mix", "Amistar", "AZX"),
  MFN = c("Mix", "Belanty", "MFN")
)

# === FILTER AND PIVOT DATA ===
long_data <- data %>%
  filter(Tissue %in% included_tissues) %>%
  pivot_longer(cols = all_of(fungicides_to_plot),
               names_to = "Fungicide",
               values_to = "Value") %>%
  mutate(Value = as.numeric(Value))

# === CALCULATE Mobile and Immobile percentages ===
mobile_data <- long_data %>%
  group_by(Fungicide, Treatment, Repetition) %>%
  summarise(
    Mobile = sum(Value[Tissue %in% c("T", "B")], na.rm = TRUE),
    Immobile = sum(Value[Tissue == "M"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total = Mobile + Immobile,
    Mobile_pct = ifelse(Total > 0, (Mobile / Total) * 100, NA_real_),
    Immobile_pct = ifelse(Total > 0, (Immobile / Total) * 100, NA_real_)
  ) %>%
  filter(!is.na(Mobile_pct))

# === FILTER & ORDER TREATMENTS ===
mobile_data <- mobile_data %>%
  filter(
    (Fungicide == "FLU" & Treatment %in% treatment_mapping$FLU) |
      (Fungicide == "CYP" & Treatment %in% treatment_mapping$CYP) |
      (Fungicide == "AZX" & Treatment %in% treatment_mapping$AZX) |
      (Fungicide == "MFN" & Treatment %in% treatment_mapping$MFN)
  ) %>%
  mutate(
    Treatment = case_when(
      Fungicide == "FLU" ~ factor(Treatment, levels = treatment_mapping$FLU),
      Fungicide == "CYP" ~ factor(Treatment, levels = treatment_mapping$CYP),
      Fungicide == "AZX" ~ factor(Treatment, levels = treatment_mapping$AZX),
      Fungicide == "MFN" ~ factor(Treatment, levels = treatment_mapping$MFN),
      TRUE ~ as.factor(Treatment)
    ),
    Fungicide = factor(Fungicide, levels = c("FLU", "CYP", "AZX", "MFN"))
  )
# === DIAGNOSTIC TESTS FOR MOBILITY DATA ===
run_diagnostic_tests <- function(fungicide_name, data) {
  data_subset <- data %>%
    filter(Fungicide == fungicide_name)
  
  if(length(unique(data_subset$Treatment)) < 2) {
    cat(sprintf("%s: Not enough treatment groups\n", fungicide_name))
    return(NULL)
  }
  
  cat(sprintf("\n=== %s Fungicide (Diagnostic Tests) ===\n", fungicide_name))
  
  # Shapiro-Wilk test for normality of residuals
  temp_anova <- aov(Mobile_pct ~ Treatment, data = data_subset)
  shapiro_test <- shapiro.test(residuals(temp_anova))
  n_samples <- length(residuals(temp_anova))
  
  cat(sprintf("Shapiro-Wilk Test (Normality):\n  W = %.4f, p = %.4f, n = %d => %s\n",
              shapiro_test$statistic, shapiro_test$p.value, n_samples,
              ifelse(shapiro_test$p.value > 0.05, "Normal", "Non-normal")))
  
  # Levene's test for homoscedasticity
  levene_test <- car::leveneTest(Mobile_pct ~ Treatment, data = data_subset)
  df1 <- levene_test$Df[1]
  df2 <- levene_test$Df[2]
  
  cat(sprintf("Levene's Test (Homoscedasticity):\n  F(%d,%d) = %.4f, p = %.4f => %s\n\n",
              df1, df2, levene_test$`F value`[1], levene_test$`Pr(>F)`[1],
              ifelse(levene_test$`Pr(>F)`[1] > 0.05, "Homogeneous", "Heterogeneous")))
}

# Run diagnostic tests for each fungicide
cat("=== DIAGNOSTIC TESTS FOR MOBILITY ACROSS TREATMENTS ===\n")
for (fungicide in levels(mobile_data$Fungicide)) {
  run_diagnostic_tests(fungicide, mobile_data)
}

# === FUNCTION TO CREATE PROPER P-VALUE MATRIX FOR multcompLetters ===
create_pvalue_matrix <- function(treatments, comparisons, p_values) {
  n <- length(treatments)
  pmat <- matrix(1, nrow = n, ncol = n)
  rownames(pmat) <- treatments
  colnames(pmat) <- treatments
  
  # Fill in p-values from comparisons
  for (i in seq_along(comparisons)) {
    pair <- strsplit(comparisons[i], " - ")[[1]]
    if (length(pair) == 2 && all(pair %in% treatments)) {
      pmat[pair[1], pair[2]] <- p_values[i]
      pmat[pair[2], pair[1]] <- p_values[i]
    }
  }
  
  # Set diagonal to 1 (no difference between same groups)
  diag(pmat) <- 1
  
  return(pmat)
}

# === RUN KRUSKAL-WALLIS + DUNN'S TESTS AND GENERATE SIGNIFICANCE LETTERS ===
run_mobility_kw_dunn <- function(fungicide_name, data) {
  data_subset <- data %>%
    filter(Fungicide == fungicide_name)
  
  if(length(unique(data_subset$Treatment)) < 2) {
    cat(sprintf("%s: Not enough treatment groups\n", fungicide_name))
    return(NULL)
  }
  
  cat(sprintf("\n=== %s Fungicide (Mobility Analysis) ===\n", fungicide_name))
  
  # Kruskal-Wallis test
  kw_test <- kruskal.test(Mobile_pct ~ Treatment, data = data_subset)
  cat(sprintf("Kruskal-Wallis: H = %.3f, df = %d, p = %.4f\n",
              kw_test$statistic, kw_test$parameter, kw_test$p.value))
  
  treatments <- levels(data_subset$Treatment)
  letters <- NULL
  
  if(kw_test$p.value < 0.05) {
    # Run Dunn's test
    dunn_res <- dunn.test(data_subset$Mobile_pct, data_subset$Treatment,
                          method = "bonferroni", alpha = 0.05, list = FALSE)
    
    # Print significant comparisons
    sig_comparisons <- data.frame(
      Comparison = dunn_res$comparisons,
      Z_score = dunn_res$Z,
      p_unadjusted = dunn_res$P,
      p_adjusted = dunn_res$P.adjusted,
      Significant = dunn_res$P.adjusted < 0.05
    )
    
    sig_pairs <- sig_comparisons[sig_comparisons$Significant, ]
    if(nrow(sig_pairs) > 0) {
      cat("Significant pairwise comparisons:\n")
      print(sig_pairs)
      
      # Create p-value matrix for multcompLetters
      pmat <- create_pvalue_matrix(treatments, dunn_res$comparisons, dunn_res$P.adjusted)
      
      # Generate compact letter display
      cld_result <- multcompView::multcompLetters(pmat, threshold = 0.05)
      letters <- cld_result$Letters
    } else {
      cat("No significant pairwise differences found\n")
      # Assign same letter to all groups
      letters <- setNames(rep("a", length(treatments)), treatments)
    }
  } else {
    cat("Kruskal-Wallis test not significant - all groups get same letter\n")
    # Assign same letter to all groups
    letters <- setNames(rep("a", length(treatments)), treatments)
  }
  
  return(list(kw = kw_test, letters = letters))
}

# Run analysis for each fungicide
kw_results_mobility <- list()
for (fungicide in levels(mobile_data$Fungicide)) {
  kw_results_mobility[[fungicide]] <- run_mobility_kw_dunn(fungicide, mobile_data)
}

# === PREPARE LETTERS DATAFRAME FOR PLOTTING ===
letter_dfs <- do.call(rbind, lapply(names(kw_results_mobility), function(fung) {
  res <- kw_results_mobility[[fung]]
  if(is.null(res) || is.null(res$letters)) return(NULL)
  
  data.frame(
    Fungicide = fung,
    Treatment = names(res$letters),
    Letters = res$letters,
    stringsAsFactors = FALSE
  )
}))

# Ensure Treatment factor levels are correctly ordered
if (!is.null(letter_dfs) && nrow(letter_dfs) > 0) {
  letter_dfs <- letter_dfs %>%
    rowwise() %>%
    mutate(Treatment = factor(Treatment, levels = treatment_mapping[[Fungicide]])) %>%
    ungroup() %>%
    mutate(Fungicide = factor(Fungicide, levels = c("FLU", "CYP", "AZX", "MFN")))
  
  # Calculate maximum Mobile_pct per fungicide-treatment to position letters above bars
  max_mobile <- mobile_data %>%
    group_by(Fungicide, Treatment) %>%
    summarise(max_Mobile_pct = max(Mobile_pct, na.rm = TRUE), .groups = "drop")
  
  letter_dfs <- left_join(letter_dfs, max_mobile, by = c("Fungicide", "Treatment")) %>%
    mutate(y_pos = max_Mobile_pct * 1.05)  # 5% above max for spacing
  
  cat("\nSignificance letters for plotting:\n")
  print(letter_dfs)
} else {
  cat("\nNo significance letters generated.\n")
  letter_dfs <- data.frame()
}

# === PLOT WITH SIGNIFICANCE LETTERS ===
custom_colors <- c("Mobile" = "green3", "Immobile" = "darkorange")

bar_data <- mobile_data %>%
  select(Fungicide, Treatment, Repetition, Mobile_pct, Immobile_pct) %>%
  pivot_longer(cols = c(Immobile_pct, Mobile_pct),
               names_to = "Mobility",
               values_to = "Percentage") %>%
  mutate(
    Mobility = case_when(
      Mobility == "Mobile_pct" ~ "Mobile",
      Mobility == "Immobile_pct" ~ "Immobile",
      TRUE ~ Mobility
    )
  )

bar_summary <- bar_data %>%
  group_by(Fungicide, Treatment, Mobility) %>%
  summarise(Mean_Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

plot2 <- ggplot() +
  geom_bar(data = bar_summary,
           aes(x = Treatment, y = Mean_Percentage, fill = Mobility),
           stat = "identity") +
  geom_jitter(data = mobile_data,
              aes(x = Treatment, y = Mobile_pct),
              color = "black", width = 0.2, size = 1.5, alpha = 0.8) +
  facet_wrap(~ Fungicide, nrow = 4, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Partitioning of Fungicide Uptake into Mobile (T+B) and Immobile (M) Fractions",
    x = "Treatment",
    y = "Percentage of Total Uptake",
    fill = "Mobility"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  # Add top margin for letters
  )
plot2

# Add significance letters if they exist
if (nrow(letter_dfs) > 0) {
  plot2 <- plot2 + 
    geom_text(data = letter_dfs,
              aes(x = Treatment, y = y_pos, label = Letters),
              inherit.aes = FALSE,
              size = 4,
              vjust = 0,
              fontface = "bold") +
    coord_cartesian(clip = "off")   # Show letters outside plot panel
}

print(plot2)

library(patchwork)
# === COMBINE ALL PLOTS IN 1 ROW, 3 COLUMNS ===
final_plot <- plot1 + plot2 + plot3 + plot_layout(nrow = 1)
final_plot


library (svglite)
ggsave("Figure 2 compare uptake mobility distribution across treatments v5.svg", plot = final_plot, width = 12, height = 8, units = "in", dpi = 300)

