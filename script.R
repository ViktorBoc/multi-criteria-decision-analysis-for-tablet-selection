# Multi-Criteria Decision Analysis for Tablet Selection
# Implements WSA, PROMETHEE, and TOPSIS methods with results validation

library(psych)
library(MCDA)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Data Preparation

load_data <- function() {
  # Load and prepare the performance table
  performance_table <- read_excel("tablety_data.xlsx", 
                                  sheet = "tablety_data", 
                                  col_names = FALSE) %>% 
    as.matrix()
  
  rownames(performance_table) <- paste0("tablet", 1:11)
  colnames(performance_table) <- c("Price", "Design", "RAM", "Display_res", 
                                   "Display_size", "Weight")
  
  # Transform cost criteria (lower is better)
  performance_table[, "Price"] <- max(performance_table[, "Price"]) - performance_table[, "Price"]
  performance_table[, "Weight"] <- max(performance_table[, "Weight"]) - performance_table[, "Weight"]
  
  return(performance_table)
}

# Saaty's Method for Weight Calculation

calculate_weights <- function() {
  # Pairwise comparison matrix
  saaty_matrix <- rbind(
    c(1,   2,   4,   5,   6,   9),
    c(1/2, 1,   2,   3,   3,   8),
    c(1/4, 1/2, 1,   3,   3,   5),
    c(1/5, 1/3, 1/3, 1,   2,   5),
    c(1/6, 1/3, 1/3, 1/2, 1,   3),
    c(1/9, 1/8, 1/5, 1/5, 1/3, 1)
  )
  
  # Calculate geometric means
  geom_means <- apply(saaty_matrix, 1, geometric.mean)
  
  # Normalize to get weights
  weights <- geom_means / sum(geom_means)
  names(weights) <- colnames(load_data())
  
  # Calculate consistency ratio (validation)
  lambda_max <- Re(eigen(saaty_matrix)$values[1])
  n <- nrow(saaty_matrix)
  ci <- (lambda_max - n) / (n - 1)
  ri <- 1.24  # Random index for n=6
  cr <- ci / ri
  
  if (cr > 0.1) {
    warning(paste("Consistency ratio", round(cr, 2), "> 0.10 - consider revising pairwise comparisons"))
  }
  
  return(list(weights = weights, consistency_ratio = cr))
}

# WSA (Weighted Sum Approach)


perform_wsa <- function(performance_table, weights) {
  # Normalize performance table
  normalized_table <- normalizePerformanceTable(
    performance_table,
    rep("rescaling", ncol(performance_table)))
  
  # Calculate utilities
  utilities <- weightedSum(normalized_table, weights)
  
  # Create score decomposition
  scores <- normalized_table * t(replicate(nrow(normalized_table), weights))
  
  # Prepare visualization data
  viz_data <- scores %>% 
    as.data.frame() %>% 
    rownames_to_column("Tablet") %>% 
    pivot_longer(-Tablet, names_to = "Criterion", values_to = "Score")
  
  # Plot
  p <- ggplot(viz_data, aes(x = Tablet, y = Score, fill = Criterion)) +
    geom_bar(stat = "identity") +
    labs(title = "WSA - Score Decomposition by Criterion",
         x = "Tablet", y = "Normalized Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  return(list(
    normalized_table = normalized_table,
    utilities = utilities,
    ranking = sort(utilities, decreasing = TRUE)
  ))
}

# 4. PROMETHEE Method

perform_promethee <- function(performance_table, weights) {
  # Define parameters
  preference_function <- c("U-shape", "V-shape", "V-shape-Indiff", 
                           "Level", "Usual", "Gaussian")
  
  preference_threshold <- c(50, 1, 2, 2, 1, 50)
  names(preference_threshold) <- colnames(performance_table)
  
  indifference_threshold <- c(20, 0, 1, 1, 0.5, 20)
  names(indifference_threshold) <- colnames(performance_table)
  
  gauss_parameter <- c(0, 0, 0, 0, 0, 10)
  names(gauss_parameter) <- colnames(performance_table)
  
  criteria_direction <- c("min", "max", "max", "max", "max", "min")
  names(criteria_direction) <- colnames(performance_table)
  
  # Perform analysis
  prometheeI <- PROMETHEEI(performance_table, preference_function,
                           preference_threshold, indifference_threshold,
                           gauss_parameter, weights, criteria_direction)
  
  prometheeII <- PROMETHEEII(performance_table, preference_function,
                             preference_threshold, indifference_threshold,
                             gauss_parameter, weights, criteria_direction)
  
  flows <- PROMETHEEOutrankingFlows(performance_table, preference_function,
                                    preference_threshold, indifference_threshold,
                                    gauss_parameter, weights, criteria_direction)
  
  net_flows <- flows$outrankingFlowsPos - flows$outrankingFlowsNeg
  
  return(list(
    prometheeI = prometheeI,
    prometheeII = prometheeII,
    net_flows = net_flows,
    ranking = sort(net_flows, decreasing = TRUE)
  ))
}

# TOPSIS Method


perform_topsis <- function(performance_table, weights) {
  # All criteria are now benefit criteria (after transformation)
  criteria_direction <- rep("max", length(weights))
  
  # Perform TOPSIS
  topsis_scores <- TOPSIS(performance_table, weights, criteria_direction)
  
  return(list(
    scores = topsis_scores,
    ranking = sort(topsis_scores, decreasing = TRUE)
  ))
}

# 6. Results Validation

validate_results <- function(wsa, promethee, topsis) {
  # Create comparison table
  comparison <- data.frame(
    Tablet = names(wsa$ranking),
    WSA_Rank = 1:length(wsa$ranking),
    WSA_Score = wsa$ranking,
    PROMETHEE_Rank = match(names(wsa$ranking), names(promethee$ranking)),
    PROMETHEE_Score = promethee$ranking[match(names(wsa$ranking), names(promethee$ranking))],
    TOPSIS_Rank = match(names(wsa$ranking), names(topsis$ranking)),
    TOPSIS_Score = topsis$ranking[match(names(wsa$ranking), names(topsis$ranking))]
  )
  
  # Calculate rank correlations
  cor_wsa_prom <- cor(comparison$WSA_Rank, comparison$PROMETHEE_Rank, method = "spearman")
  cor_wsa_topsis <- cor(comparison$WSA_Rank, comparison$TOPSIS_Rank, method = "spearman")
  cor_prom_topsis <- cor(comparison$PROMETHEE_Rank, comparison$TOPSIS_Rank, method = "spearman")
  
  # Visual comparison
  viz_data <- comparison %>%
    select(Tablet, WSA_Rank, PROMETHEE_Rank, TOPSIS_Rank) %>%
    pivot_longer(-Tablet, names_to = "Method", values_to = "Rank")
  
  p <- ggplot(viz_data, aes(x = Tablet, y = Rank, color = Method, group = Method)) +
    geom_point(size = 3) +
    geom_line() +
    scale_y_reverse() +  # Lower rank = better
    labs(title = "Comparison of Rankings Across Methods",
         x = "Tablet", y = "Rank (1 = Best)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  return(list(
    comparison_table = comparison,
    correlations = data.frame(
      Comparison = c("WSA-PROMETHEE", "WSA-TOPSIS", "PROMETHEE-TOPSIS"),
      Spearman = c(cor_wsa_prom, cor_wsa_topsis, cor_prom_topsis)
    )
  ))
}

# Main Analysis

# Load and prepare data
performance_data <- load_data()

# Calculate weights using Saaty's method
weight_results <- calculate_weights()
weights <- weight_results$weights
print(paste("Consistency Ratio:", round(weight_results$consistency_ratio, 3)))

# Perform WSA analysis
wsa_results <- perform_wsa(performance_data, weights)

# Perform PROMETHEE analysis
promethee_results <- perform_promethee(performance_data, weights)

# Perform TOPSIS analysis
topsis_results <- perform_topsis(performance_data, weights)

# Validate and compare results
validation <- validate_results(wsa_results, promethee_results, topsis_results)

# Print key results
print("Final Rankings:")
print(validation$comparison_table)

print("Rank Correlations:")
print(validation$correlations)
