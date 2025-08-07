Figure 3 (3 plots: plot1 uptake, plot2 mobility, plot3 direction of translocation)
# === LOAD LIBRARIES ===plot3
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(stringr)
library(car)            # For Levene's test
library(multcompView)   # For significance letters
library(dunn.test)      # For Dunn's test

# === LOAD DATA ===
data <- read_csv("FORM5 processed data.csv", skip = 1)
colnames(data) <- str_trim(colnames(data))

# === PARAMETERS ===
fungicides_to_plot <- c("TRC", "MXL", "FLU", "CYP", "BSC", "FLX", 
                        "EPX", "TEB", "PCZ", "AZX", "MFN", 
                        "MND", "FPR", "PYR")
included_tissues <- c("T", "B")

# === FILTER FOR MIX ONLY and DAT 3, 7, 14 ===
long_data <- data %>%
  filter(Treatment == "Mix", Tissue %in% included_tissues, DAT %in% c(3, 7, 14)) %>%
  pivot_longer(
    cols = all_of(fungicides_to_plot),
    names_to = "Fungicide",
    values_to = "Value"
  ) %>%
  mutate(Value = as.numeric(Value))

# === CALCULATE T+B TOTAL PER GROUP ===
total_TB <- long_data %>%
  group_by(Repetition, Fungicide, DAT) %>%
  summarise(Total_TB = sum(Value, na.rm = TRUE), .groups = "drop")

# === MERGE AND COMPUTE PROPORTIONS ===
long_data <- long_data %>%
  left_join(total_TB, by = c("Repetition", "Fungicide", "DAT")) %>%
  mutate(Proportion = (Value / Total_TB) * 100) %>%
  filter(Total_TB > 0)

# === FACTOR LEVELS ===
long_data$Fungicide <- factor(long_data$Fungicide, levels = fungicides_to_plot)
long_data$DAT <- factor(long_data$DAT, levels = c(3, 7, 14))

# === SHAPIRO-WILK AND LEVENE'S TEST FOR B TISSUE (Pooled DATs) ===
cat("\n=== Shapiro-Wilk and Levene's Tests on B Tissue Proportions ===\n")

b_data <- long_data %>% filter(Tissue == "B")

# Shapiro-Wilk test (normality of residuals from one-way ANOVA model)
shapiro_model <- aov(Proportion ~ Fungicide, data = b_data)
shapiro_test <- shapiro.test(residuals(shapiro_model))
n_residuals <- length(residuals(shapiro_model))

cat(sprintf("Shapiro-Wilk Normality Test: W = %.4f, n = %d, p = %.4f (%s)\n",
            shapiro_test$statistic,
            n_residuals,
            shapiro_test$p.value,
            ifelse(shapiro_test$p.value > 0.05, "Normal", "Non-normal")))

# Levene's Test for homoscedasticity
levene_result <- car::leveneTest(Proportion ~ Fungicide, data = b_data)
cat(sprintf("Levene's Test: F(%d,%d) = %.4f, p = %.4f (%s)\n",
            levene_result$Df[1],
            levene_result$Df[2],
            levene_result$`F value`[1],
            levene_result$`Pr(>F)`[1],
            ifelse(levene_result$`Pr(>F)`[1] > 0.05, "Homogeneous", "Heterogeneous")))

cat(rep("=", 60), "\n")

# === KRUSKAL-WALLIS AND DUNN'S TESTS FOR B% PROPORTIONS ===
cat("\n=== Kruskal-Wallis Test for B% Proportions Across Fungicides ===\n")

kw_test <- kruskal.test(Proportion ~ Fungicide, data = b_data)
cat(sprintf("Kruskal-Wallis test: H = %.3f, df = %d, p = %.4f\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))

# Test for DAT effect
kw_dat <- kruskal.test(Proportion ~ DAT, data = b_data)
cat(sprintf("Kruskal-Wallis test (DAT): H = %.3f, df = %d, p = %.4f\n",
            kw_dat$statistic, kw_dat$parameter, kw_dat$p.value))

# Default: assign all "a"
letter_df <- data.frame(Fungicide = fungicides_to_plot, Letters = "a", stringsAsFactors = FALSE)

if (kw_test$p.value < 0.05) {
  cat("\n--- Dunn's Post Hoc Test (Bonferroni adjusted) ---\n")
  dunn_result <- dunn.test(x = b_data$Proportion, g = b_data$Fungicide, method = "bonferroni", list = FALSE)
  
  # Calculate number of groups and comparisons
  n_groups <- length(unique(b_data$Fungicide))
  n_comparisons <- n_groups * (n_groups - 1) / 2
  
  cat(sprintf("Number of groups: %d\n", n_groups))
  cat(sprintf("Total pairwise comparisons: %d\n", n_comparisons))
  cat(sprintf("Multiple comparison correction: Bonferroni\n"))
  
  # Extract significant results with all parameters
  sig_results <- data.frame(
    Comparison = dunn_result$comparisons,
    Z_statistic = dunn_result$Z,
    P_unadjusted = dunn_result$P,
    P_adjusted = dunn_result$P.adjusted
  ) %>% filter(P_adjusted < 0.05)
  
  cat(sprintf("\nSignificant pairwise differences (p < 0.05): %d out of %d comparisons\n", 
              nrow(sig_results), n_comparisons))
  
  if (nrow(sig_results) > 0) {
    cat("\nDetailed results for significant comparisons:\n")
    for(i in 1:nrow(sig_results)) {
      cat(sprintf("  %s: Z = %.3f, p_unadj = %.4f, p_adj = %.4f\n", 
                  sig_results$Comparison[i], 
                  sig_results$Z_statistic[i], 
                  sig_results$P_unadjusted[i],
                  sig_results$P_adjusted[i]))
    }
    
    # Create p-matrix for letter generation
    all_fungicides <- levels(factor(b_data$Fungicide, levels = fungicides_to_plot))
    p_matrix <- matrix(1, nrow = length(all_fungicides), ncol = length(all_fungicides),
                       dimnames = list(all_fungicides, all_fungicides))
    
    for (i in seq_along(dunn_result$comparisons)) {
      comparison <- strsplit(dunn_result$comparisons[i], " - ")[[1]]
      if (length(comparison) == 2) {
        f1 <- comparison[1]
        f2 <- comparison[2]
        p_val <- dunn_result$P.adjusted[i]
        p_matrix[f1, f2] <- p_val
        p_matrix[f2, f1] <- p_val
      }
    }
    
    letters_out <- multcompLetters(p_matrix < 0.05)$Letters
    letter_df <- data.frame(Fungicide = names(letters_out), Letters = letters_out, stringsAsFactors = FALSE)
    
    cat("\nCompact letter display generated for visualization\n")
    
  } else {
    cat("\nNo significant pairwise differences detected after Bonferroni correction.\n")
  }
  
} else {
  cat("\nKruskal-Wallis test not significant (p >= 0.05). No post-hoc testing performed.\n")
}

# === CALCULATE B% PROPORTIONS FOR BAR PLOT ===
b_bar_data <- b_data %>%
  group_by(Fungicide, DAT) %>%
  summarise(
    Mean_B_Proportion = mean(Proportion, na.rm = TRUE),
    SE = sd(Proportion, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(Mean_B_Proportion))

# Add significance letters
b_bar_data <- left_join(b_bar_data, letter_df, by = "Fungicide")

# === PLOT ===
log1p_trans <- scales::trans_new(
  name = "log1p",
  transform = log1p,
  inverse = expm1,
  breaks = function(x) c(0, 1, 3, 5, 10, 30, 50, 100),
  domain = c(0, Inf)
)

b_bar_data$DAT <- factor(b_bar_data$DAT, levels = c(3, 7, 14))
b_bar_data$Fungicide <- factor(b_bar_data$Fungicide, levels = fungicides_to_plot)

plot3 <- ggplot(b_bar_data, aes(x = DAT, y = Mean_B_Proportion, fill = DAT)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_jitter(
    data = long_data %>% filter(Tissue == "B"),
    aes(x = DAT, y = Proportion, color = DAT),
    width = 0.2, size = 1.5, alpha = 0.7, inherit.aes = FALSE
  ) +
  geom_text(
    data = b_bar_data %>% filter(DAT == 7 & !is.na(Letters)),
    aes(x = DAT, y = Mean_B_Proportion + 5, label = Letters),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ Fungicide, nrow = 1, ncol = 15, drop = FALSE) +
  scale_y_continuous(
    trans = log1p_trans,
    limits = c(0, 100)
  ) +
  scale_fill_brewer(palette = "Set2", name = "Days After Treatment") +
  scale_color_brewer(palette = "Set2", name = "Days After Treatment") +
  labs(
    y = "Mean B% Proportion (log(x + 1) scale)",
    title = "B% Proportion Across Time (Mix Treatment Only)",
    caption = "Bars = Mean; Points = Individual Reps; Letters = Significance"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, face = "bold", size = 11)
  )

print(plot3)



# === LOAD LIBRARIES ===plot1
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(dunn.test)  # For Dunn's test
library(car)        # For Levene's test

# === LOAD DATA ===
data <- read_csv("FORM5 processed data.csv", skip = 1)
colnames(data) <- str_trim(colnames(data))
formulation_data <- read_csv("FORM5 formulation quant.csv")
colnames(formulation_data) <- str_trim(colnames(formulation_data))

# === APPLIED AMOUNTS ===
amount_applied <- formulation_data[1, -1]
amount_applied <- as.numeric(amount_applied)
names(amount_applied) <- colnames(formulation_data)[-1]

# === PARAMETERS ===
included_tissues <- c("T", "M", "B")
fungicide_order <- c("TRC", "MXL", "FLU", "CYP", "BSC", "FLX", 
                     "EPX", "TEB", "PCZ", "AZX", "MFN", 
                     "MND", "FPR", "PYR")

# === FILTER & LONG FORMAT ===
long_data <- data %>%
  filter(Tissue %in% included_tissues,
         Treatment == "Mix",
         DAT %in% c(3, 7, 14)) %>%
  pivot_longer(cols = all_of(fungicide_order),
               names_to = "Fungicide",
               values_to = "Value") %>%
  mutate(Value = as.numeric(Value))

# === CALCULATE TOTAL UPTAKE PER SAMPLE ===
sum_uptake <- long_data %>%
  group_by(Fungicide, Repetition, DAT) %>%
  summarise(Total_Uptake = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Applied = amount_applied[Fungicide]) %>%
  mutate(Uptake_Percent = (Total_Uptake / Applied) * 100)

# === FACTOR LEVELS ===
sum_uptake$Fungicide <- factor(sum_uptake$Fungicide, levels = fungicide_order)
sum_uptake$DAT <- factor(sum_uptake$DAT, levels = c(3, 7, 14))

# === PARAMETRIC ASSUMPTIONS TESTING ===
cat("=== TESTING PARAMETRIC ASSUMPTIONS ===\n")
temp_anova <- aov(Uptake_Percent ~ Fungicide + DAT, data = sum_uptake)

# Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(residuals(temp_anova))
n_samples <- length(residuals(temp_anova))

cat(sprintf("Shapiro-Wilk: W = %.4f, p = %.4f, n = %d (%s)\n", 
            shapiro_test$statistic, shapiro_test$p.value, n_samples,
            ifelse(shapiro_test$p.value > 0.05, "Normal", "Non-normal")))

# Levene's test for homoscedasticity
levene_test <- car::leveneTest(Uptake_Percent ~ Fungicide, data = sum_uptake)
n_samples <- nrow(sum_uptake)
df1 <- levene_test$Df[1]  # numerator degrees of freedom
df2 <- levene_test$Df[2]  # denominator degrees of freedom

cat(sprintf("Levene's: F(%d,%d) = %.4f, p = %.4f, n = %d (%s)\n", 
            df1, df2, levene_test$`F value`[1], levene_test$`Pr(>F)`[1], n_samples,
            ifelse(levene_test$`Pr(>F)`[1] > 0.05, "Homogeneous", "Heterogeneous")))
cat("Using non-parametric tests due to assumption violations\n\n")

# === KRUSKAL-WALLIS TEST ===
cat("=== KRUSKAL-WALLIS TEST RESULTS ===\n")
kw_fungicide <- kruskal.test(Uptake_Percent ~ Fungicide, data = sum_uptake)
kw_dat <- kruskal.test(Uptake_Percent ~ DAT, data = sum_uptake)

cat(sprintf("Fungicide effect: H = %.3f, df = %d, p = %.4f\n", 
            kw_fungicide$statistic, kw_fungicide$parameter, kw_fungicide$p.value))
cat(sprintf("DAT effect: H = %.3f, df = %d, p = %.4f\n\n", 
            kw_dat$statistic, kw_dat$parameter, kw_dat$p.value))

# === DUNN'S POST-HOC TEST ===
if(kw_fungicide$p.value < 0.05) {
  dunn_result <- dunn.test(sum_uptake$Uptake_Percent, sum_uptake$Fungicide, 
                           method = "bonferroni", alpha = 0.05, list = FALSE)
  
  # Create significance matrix
  fungicide_levels <- levels(sum_uptake$Fungicide)
  n_fungi <- length(fungicide_levels)
  total_comparisons <- n_fungi * (n_fungi - 1) / 2
  sig_matrix <- matrix(FALSE, n_fungi, n_fungi)
  rownames(sig_matrix) <- colnames(sig_matrix) <- fungicide_levels
  
  # Create comparison results dataframe
  comparisons <- data.frame(
    Comparison = dunn_result$comparisons,
    Z_score = dunn_result$Z,
    P_unadj = dunn_result$P,
    P_adj = dunn_result$P.adjusted,
    Significant = dunn_result$P.adjusted < 0.05
  )
  
  # Print summary information
  cat(sprintf("Number of groups: %d\n", n_fungi))
  cat(sprintf("Total pairwise comparisons: %d\n", total_comparisons))
  cat("Multiple comparison correction: Bonferroni\n\n")
  
  # Count and report significant differences
  n_significant <- sum(comparisons$Significant)
  cat(sprintf("Significant pairwise differences (p < 0.05): %d out of %d comparisons\n\n", 
              n_significant, total_comparisons))
  
  # Print detailed results for significant comparisons
  if(n_significant > 0) {
    cat("Detailed results for significant comparisons:\n")
    sig_comparisons <- comparisons[comparisons$Significant, ]
    
    for(i in 1:nrow(sig_comparisons)) {
      cat(sprintf("  %s: Z = %.3f, p_unadj = %.4f, p_adj = %.4f\n", 
                  sig_comparisons$Comparison[i],
                  sig_comparisons$Z_score[i],
                  sig_comparisons$P_unadj[i],
                  sig_comparisons$P_adj[i]))
    }
    cat("\n")
  } else {
    cat("No significant pairwise differences found\n\n")
  }
  
  # Fill matrix with significant differences (for letter assignment)
  for(i in 1:nrow(comparisons)) {
    if(comparisons$Significant[i]) {
      fungi_pair <- strsplit(comparisons$Comparison[i], " - ")[[1]]
      if(length(fungi_pair) == 2 && all(fungi_pair %in% fungicide_levels)) {
        sig_matrix[fungi_pair[1], fungi_pair[2]] <- TRUE
        sig_matrix[fungi_pair[2], fungi_pair[1]] <- TRUE
      }
    }
  }
  
  # Generate proper letter groupings
  median_order <- sum_uptake %>%
    group_by(Fungicide) %>%
    summarise(Median = median(Uptake_Percent, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Median)) %>%
    pull(Fungicide)
  
  letters_assigned <- rep("", n_fungi)
  names(letters_assigned) <- fungicide_levels
  current_letter <- 1
  
  for(fungicide in median_order) {
    if(letters_assigned[fungicide] == "") {
      # Find fungi not significantly different from this one
      not_sig_diff <- names(which(!sig_matrix[fungicide, ]))
      not_sig_diff <- not_sig_diff[letters_assigned[not_sig_diff] == ""]
      
      # Assign same letter to all non-significantly different fungi
      letters_assigned[not_sig_diff] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  fungicide_letters <- data.frame(
    Fungicide = names(letters_assigned),
    Letters = letters_assigned
  )
  
} else {
  fungicide_letters <- data.frame(
    Fungicide = fungicide_order,
    Letters = "a"
  )
  cat("Kruskal-Wallis test not significant - no post-hoc testing performed\n\n")
}

# === SUMMARY ===
fungicide_summary <- sum_uptake %>%
  group_by(Fungicide) %>%
  summarise(
    Median_Uptake = round(median(Uptake_Percent, na.rm = TRUE), 2),
    IQR = round(IQR(Uptake_Percent, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  left_join(fungicide_letters, by = "Fungicide") %>%
  arrange(desc(Median_Uptake))

print(fungicide_summary)

# === PLOT ===
mean_uptake <- sum_uptake %>%
  group_by(Fungicide, DAT) %>%
  summarise(Median_Uptake = median(Uptake_Percent, na.rm = TRUE), .groups = "drop") %>%
  left_join(fungicide_letters, by = "Fungicide")

mean_uptake$Fungicide <- factor(mean_uptake$Fungicide, levels = fungicide_order)
mean_uptake$DAT <- factor(mean_uptake$DAT, levels = c(3, 7, 14))

log1p_trans <- scales::trans_new(
  name = "log1p",
  transform = log1p,
  inverse = expm1,
  breaks = function(x) c(0, 1, 3, 5, 10, 30, 50, 100),
  domain = c(0, Inf)
)

plot1 <- ggplot(mean_uptake, aes(x = DAT, y = Median_Uptake, fill = DAT)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_jitter(
    data = sum_uptake,
    aes(x = DAT, y = Uptake_Percent, color = DAT),
    width = 0.2, size = 1.5, alpha = 0.7,
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(x = 2, y = 90, label = Letters),
    size = 4, fontface = "bold", color = "red",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Fungicide, nrow = 1) +
  scale_y_continuous(trans = log1p_trans, limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set2", name = "Days After Treatment") +
  scale_color_brewer(palette = "Set2", name = "Days After Treatment") +
  labs(
    x = "Days After Treatment (DAT)",
    y = "Uptake (% of Applied Amount)\n[log(y + 1) scale]",
    title = "Fungicide Uptake Over Time",
    subtitle = paste0("Kruskal-Wallis: H = ", round(kw_fungicide$statistic, 2), 
                      ", p = ", round(kw_fungicide$p.value, 4)),
    caption = "Red letters = significance groups (same letter = not significantly different)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(plot1)

# === SAVE RESULTS ===
write.csv(fungicide_summary, "fungicide_summary_kruskal.csv", row.names = FALSE)
if(exists("comparisons")) {
  write.csv(comparisons, "fungicide_dunn_comparisons.csv", row.names = FALSE)
}

test_results <- data.frame(
  Test = c("Shapiro-Wilk", "Levene's", "Kruskal-Wallis (Fungicide)", "Kruskal-Wallis (DAT)"),
  Statistic = c(shapiro_test$statistic, levene_test$`F value`[1], 
                kw_fungicide$statistic, kw_dat$statistic),
  P_value = c(shapiro_test$p.value, levene_test$`Pr(>F)`[1], 
              kw_fungicide$p.value, kw_dat$p.value)
)

write.csv(test_results, "plot1statistical_tests.csv", row.names = FALSE)

cat("\nResults saved to CSV files\n")


# === LOAD LIBRARIES === plot2
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(car)        # For Levene's test
library(dunn.test)  # For Dunn's test

# === LOAD DATA ===
data <- read_csv("FORM5 processed data.csv", skip = 1)
colnames(data) <- str_trim(colnames(data))

# === PARAMETERS ===
fungicide_order <- c("TRC", "MXL", "FLU", "CYP", "BSC", "FLX", 
                     "EPX", "TEB", "PCZ", "AZX", "MFN", 
                     "MND", "FPR", "PYR")
included_tissues <- c("T", "M", "B")
selected_DATs <- c(3, 7, 14)

# === FILTER AND PIVOT DATA ===
long_data <- data %>%
  filter(Treatment == "Mix", DAT %in% selected_DATs, Tissue %in% included_tissues) %>%
  pivot_longer(cols = all_of(fungicide_order),
               names_to = "Fungicide",
               values_to = "Value") %>%
  mutate(Value = as.numeric(Value))

# === CALCULATE MOBILE/IMMOBILE % ===
mobile_data <- long_data %>%
  group_by(DAT, Fungicide, Repetition) %>%
  summarise(
    Mobile = sum(Value[Tissue %in% c("T", "B")], na.rm = TRUE),
    Immobile = sum(Value[Tissue == "M"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total = Mobile + Immobile,
    Mobile_pct = ifelse(Total > 0, (Mobile / Total) * 100, NA_real_),
    Immobile_pct = ifelse(Total > 0, (Immobile / Total) * 100, NA_real_),
    Fungicide = factor(Fungicide, levels = fungicide_order)
  ) %>%
  filter(!is.na(Mobile_pct))

# === DIAGNOSTIC TESTS ===
temp_anova <- aov(Mobile_pct ~ Fungicide, data = mobile_data)
shapiro_test <- shapiro.test(residuals(temp_anova))
levene_test <- car::leveneTest(Mobile_pct ~ Fungicide, data = mobile_data)
n_samples <- length(residuals(temp_anova))

# === REPORT TEST RESULTS ===
cat("==== DIAGNOSTIC TEST RESULTS ====\n")
cat(sprintf("Shapiro-Wilk Test (Normality of residuals):\n  W = %.4f, p = %.4f, n = %d => %s\n",
            shapiro_test$statistic, shapiro_test$p.value, n_samples,
            ifelse(shapiro_test$p.value > 0.05, "Normal", "Non-normal")))
cat(sprintf("Levene's Test (Homoscedasticity):\n  F(%d, %d) = %.4f, p = %.4f => %s\n",
            levene_test$Df[1], levene_test$Df[2], levene_test$`F value`[1],
            levene_test$`Pr(>F)`[1],
            ifelse(levene_test$`Pr(>F)`[1] > 0.05, "Homogeneous", "Heterogeneous")))
cat("==================================\n\n")

# === NON-PARAMETRIC TESTS ===
kw_test <- kruskal.test(Mobile_pct ~ Fungicide, data = mobile_data)
kw_dat <- kruskal.test(Mobile_pct ~ DAT, data = mobile_data)

# Report Kruskal-Wallis results
cat(sprintf("Kruskal-Wallis test (Fungicide): H = %.3f, df = %d, p = %.4f\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))
cat(sprintf("Kruskal-Wallis test (DAT): H = %.3f, df = %d, p = %.4f\n\n",
            kw_dat$statistic, kw_dat$parameter, kw_dat$p.value))

if (kw_test$p.value < 0.05) {
  dunn_result <- dunn.test(mobile_data$Mobile_pct, mobile_data$Fungicide, 
                           method = "bonferroni", alpha = 0.05, list = FALSE)
  
  fungicide_levels <- levels(mobile_data$Fungicide)
  n_fungi <- length(fungicide_levels)
  total_comparisons <- n_fungi * (n_fungi - 1) / 2
  sig_matrix <- matrix(FALSE, n_fungi, n_fungi,
                       dimnames = list(fungicide_levels, fungicide_levels))
  
  # Create comparison results dataframe
  comparisons <- data.frame(
    Comparison = dunn_result$comparisons,
    Z_score = dunn_result$Z,
    P_unadj = dunn_result$P,
    P_adj = dunn_result$P.adjusted,
    Significant = dunn_result$P.adjusted < 0.05
  )
  
  # Print summary information
  cat(sprintf("FUNGICIDE POST-HOC ANALYSIS:\n"))
  cat(sprintf("Number of groups: %d\n", n_fungi))
  cat(sprintf("Total pairwise comparisons: %d\n", total_comparisons))
  cat("Multiple comparison correction: Bonferroni\n\n")
  
  # Count and report significant differences
  n_significant <- sum(comparisons$Significant)
  cat(sprintf("Significant pairwise differences (p < 0.05): %d out of %d comparisons\n\n", 
              n_significant, total_comparisons))
  
  # Print detailed results for significant comparisons
  if(n_significant > 0) {
    cat("Detailed results for significant comparisons:\n")
    sig_comparisons <- comparisons[comparisons$Significant, ]
    
    for(i in 1:nrow(sig_comparisons)) {
      cat(sprintf("  %s: Z = %.3f, p_unadj = %.4f, p_adj = %.4f\n", 
                  sig_comparisons$Comparison[i],
                  sig_comparisons$Z_score[i],
                  sig_comparisons$P_unadj[i],
                  sig_comparisons$P_adj[i]))
    }
    cat("\n")
  } else {
    cat("No significant pairwise differences found\n\n")
  }
  
  # Fill significance matrix for letter assignment
  for (i in seq_along(dunn_result$comparisons)) {
    comp <- strsplit(dunn_result$comparisons[i], " - ")[[1]]
    if (dunn_result$P.adjusted[i] < 0.05) {
      sig_matrix[comp[1], comp[2]] <- TRUE
      sig_matrix[comp[2], comp[1]] <- TRUE
    }
  }
  
  median_order <- mobile_data %>%
    group_by(Fungicide) %>%
    summarise(Median = median(Mobile_pct), .groups = "drop") %>%
    arrange(desc(Median)) %>%
    pull(Fungicide) %>%
    as.character()
  
  letters_assigned <- setNames(rep("", length(fungicide_levels)), fungicide_levels)
  current_letter <- 1
  
  for (fungi in median_order) {
    if (letters_assigned[fungi] == "") {
      group <- names(which(!sig_matrix[fungi, ]))
      group <- group[letters_assigned[group] == ""]
      letters_assigned[group] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  fungicide_letters <- data.frame(
    Fungicide = names(letters_assigned),
    Letters = letters_assigned,
    stringsAsFactors = FALSE
  )
} else {
  fungicide_letters <- data.frame(Fungicide = fungicide_order, Letters = "a")
  cat("Kruskal-Wallis test (Fungicide) not significant - no post-hoc testing performed\n\n")
}

# Optional: DAT post-hoc analysis if significant
if (kw_dat$p.value < 0.05) {
  cat("DAT effect is significant - consider post-hoc analysis if needed\n\n")
}
fungicide_letters$Fungicide <- factor(fungicide_letters$Fungicide, levels = fungicide_order)

# === LONG FORMAT FOR STACKED BAR ===
bar_data <- mobile_data %>%
  select(Fungicide, DAT, Repetition, Mobile_pct, Immobile_pct) %>%
  pivot_longer(cols = c(Mobile_pct, Immobile_pct),
               names_to = "Mobility",
               values_to = "Percentage") %>%
  mutate(
    Mobility = ifelse(Mobility == "Mobile_pct", "Mobile", "Immobile"),
    Fungicide = factor(Fungicide, levels = fungicide_order)
  )

# === MEAN FOR STACKED BARS ===
bar_summary <- bar_data %>%
  group_by(Fungicide, DAT, Mobility) %>%
  summarise(Mean_Percentage = mean(Percentage), .groups = "drop") %>%
  left_join(fungicide_letters, by = "Fungicide") %>%
  mutate(Fungicide = factor(Fungicide, levels = fungicide_order))

# === CUSTOM COLORS ===
custom_colors <- c("Mobile" = "green3", "Immobile" = "darkorange")

# === FINAL PLOT ===
plot2 <- ggplot() +
  geom_bar(
    data = bar_summary,
    aes(x = as.factor(DAT), y = Mean_Percentage, fill = Mobility),
    stat = "identity"
  ) +
  geom_jitter(
    data = mobile_data,
    aes(x = as.factor(DAT), y = Mobile_pct),
    color = "black",
    width = 0.2,
    size = 1.5,
    alpha = 0.8
  ) +
  geom_text(
    data = bar_summary %>% filter(Mobility == "Mobile" & DAT == 14),
    aes(x = as.factor(DAT), y = 105, label = Letters),
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~ Fungicide, ncol = 15, scales = "fixed", drop = FALSE) +
  labs(
    title = "Mobile vs Immobile % Uptake Across DATs (Mix Treatment Only)",
    subtitle = paste0("Normality: p = ", signif(shapiro_test$p.value, 3),
                      " | Homoscedasticity: p = ", signif(levene_test$`Pr(>F)`[1], 3)),
    x = "Days After Treatment (DAT)",
    y = "Percentage of Total Uptake",
    fill = "Mobility"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

# === DISPLAY PLOT ===
print(plot2)

# === COMBINE ALL PLOTS IN 1 ROW, 3 COLUMNS ===
library(patchwork)
final_plot <- plot1 + plot2 + plot3 + plot_layout(nrow = 3, ncol=1)
final_plot


library (svglite)
ggsave("Figure 3 compare uptake mobility distribution of all fungicides across 3 7 14 DAT stacked boxplots v4.svg", plot = final_plot, width = 12, height = 10.5, units = "in", dpi = 300)

