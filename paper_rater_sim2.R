# ============================================================
# Simulation of Grading: Random vs LP-based Matched Assignment
# ============================================================

library(tidyverse)
library(dplyr)
library(fixest)
library(lpSolve)

# ---------------------------
# 1. Define Simulation Parameters
# ---------------------------

n_papers <- 1000
min_papers <- 30
max_papers <- 150

grader_profiles <- list(
  good = list(dist_type = "normal", params = list(bias = 0, sd = 0.5)),
  bad  = list(dist_type = "gamma",  params = list(shape = 2, scale = 1))
)
grader_proportions <- c(good = 0.7, bad = 0.3)
rating_probs <- c(0.3, 0.2, 0.15, 0.12, 0.1, 0.08, 0.05)

# Determine number of graders (using a load-based heuristic)
n_graders <- ceiling(2 * (1.5 * n_papers) / max_papers)
if(n_graders * min_papers > 2 * n_papers) {
  stop("Parameters violate minimum load constraints.")
}
cat("n_graders =", n_graders, "\n")

# ---------------------------
# 2. Define Helper Functions
# ---------------------------

generate_papers <- function(n_papers, rating_probs) {
  true_ratings <- sample(1:7, size = n_papers, replace = TRUE, prob = rating_probs)
  data.frame(paper_id = 1:n_papers, true_rating = true_ratings)
}

generate_graders <- function(n_graders, grader_profiles, grader_proportions) {
  types <- names(grader_profiles)
  assigned_types <- sample(types, size = n_graders, replace = TRUE, prob = grader_proportions)
  dist_types <- sapply(assigned_types, function(x) grader_profiles[[x]]$dist_type)
  dist_params <- lapply(assigned_types, function(x) grader_profiles[[x]]$params)
  graders <- data.frame(grader_id = 1:n_graders,
                        grader_type = assigned_types,
                        dist_type = dist_types,
                        stringsAsFactors = FALSE)
  graders$dist_params <- I(dist_params)
  return(graders)
}

assign_papers <- function(n_papers, n_graders, min_papers, max_papers) {
  total_assignments <- 2 * n_papers
  m <- n_graders
  if(total_assignments < m * min_papers || total_assignments > m * max_papers)
    stop("No feasible assignment with these parameters.")
  
  loads <- rep(min_papers, m)
  remaining <- total_assignments - sum(loads)
  max_extra <- max_papers - min_papers
  for(i in 1:remaining) {
    candidates <- which(loads - min_papers < max_extra)
    selected <- sample(candidates, 1)
    loads[selected] <- loads[selected] + 1
  }
  if(sum(loads) != total_assignments) stop("Load assignment error.")
  
  grader_vector <- rep(1:m, times = loads)
  grader_vector <- sample(grader_vector)
  
  assignments <- data.frame(
    paper_id = rep(1:n_papers, each = 2),
    grader_id = grader_vector
  )
  for(i in 1:n_papers) {
    rows_i <- which(assignments$paper_id == i)
    if(assignments$grader_id[rows_i[1]] == assignments$grader_id[rows_i[2]]) {
      for(j in 1:n_papers) {
        if(j == i) next
        rows_j <- which(assignments$paper_id == j)
        temp <- assignments$grader_id[rows_i[1]]
        assignments$grader_id[rows_i[1]] <- assignments$grader_id[rows_j[1]]
        assignments$grader_id[rows_j[1]] <- temp
        if(length(unique(assignments$grader_id[rows_i])) == 2 &&
           length(unique(assignments$grader_id[rows_j])) == 2) break
        else {
          temp <- assignments$grader_id[rows_i[1]]
          assignments$grader_id[rows_i[1]] <- assignments$grader_id[rows_j[1]]
          assignments$grader_id[rows_j[1]] <- temp
        }
      }
    }
  }
  if(any(sapply(split(assignments$grader_id, assignments$paper_id), function(x) length(unique(x)) != 2)))
    stop("Some papers have duplicate graders.")
  return(assignments)
}

simulate_score <- function(true_rating, dist_type, params) {
  if(dist_type == "normal") {
    error <- rnorm(1, mean = params$bias, sd = params$sd)
    score <- true_rating + error
  } else if(dist_type == "gamma") {
    error <- rgamma(1, shape = params$shape, scale = params$scale) - (params$shape * params$scale)
    score <- true_rating + error
  } else stop("Unknown distribution type.")
  score <- round(score)
  score <- min(max(score, 1), 7)
  return(score)
}

simulate_grading <- function(papers, assignments, graders) {
  sim_data <- merge(assignments, papers, by = "paper_id")
  sim_data <- merge(sim_data, graders, by = "grader_id")
  sim_data <- sim_data %>%
    rowwise() %>%
    mutate(simulated_score = simulate_score(true_rating, dist_type, dist_params)) %>%
    ungroup()
  return(sim_data)
}

compute_benchmarks <- function(sim_data) {
  sim_data <- sim_data %>%
    mutate(error = simulated_score - true_rating,
           abs_error = abs(error),
           sq_error = error^2)
  
  mae <- mean(sim_data$abs_error)
  mse <- mean(sim_data$sq_error)
  classification_rate <- mean(sim_data$simulated_score == sim_data$true_rating)
  
  benchmarks <- list(mae = mae,
                     mse = mse,
                     classification_rate = classification_rate)
  return(benchmarks)
}


# ----- LP-based Matching Assignment Function -----

compute_loads <- function(n_papers, n_graders, min_papers, max_papers) {
  total_assignments <- 2 * n_papers
  loads <- rep(min_papers, n_graders)
  remaining <- total_assignments - sum(loads)
  max_extra <- max_papers - min_papers
  for (i in 1:remaining) {
    candidates <- which(loads - min_papers < max_extra)
    selected <- sample(candidates, 1)  # Sample once and store it.
    loads[selected] <- loads[selected] + 1
  }
  return(loads)
}


assign_papers_matched_lp <- function(n_papers, graders, min_papers, max_papers) {
  m <- nrow(graders)
  loads <- compute_loads(n_papers, m, min_papers, max_papers)
  if(sum(loads) != 2 * n_papers) stop("Load mismatch.")
  
  # Build list of all unique pairs (a, b) with a < b.
  pair_list <- list()
  pair_index <- 1
  
  # Create a progress bar for the outer loop (over 'a').
  pb <- txtProgressBar(min = 1, max = m - 1, style = 3)
  for(a in 1:(m - 1)) {
    for(b in (a + 1):m) {
      pair_list[[pair_index]] <- c(a, b)
      pair_index <- pair_index + 1
    }
    setTxtProgressBar(pb, a)
  }
  close(pb)
  
  N <- length(pair_list)  # Number of pair variables.
  
  # Objective: cost for pair (a, b) is |loads[a] - loads[b]|
  cost <- sapply(pair_list, function(p) { abs(loads[p[1]] - loads[p[2]]) })
  
  # Build constraint matrix.
  A <- matrix(0, nrow = m, ncol = N)
  for (j in 1:N) {
    pair <- pair_list[[j]]
    A[pair[1], j] <- 1
    A[pair[2], j] <- 1
  }
  b_vec <- loads
  
  # Total papers constraint: sum(y) = n_papers.
  A_total <- matrix(1, nrow = 1, ncol = N)
  b_total <- n_papers
  
  A_full <- rbind(A, A_total)
  b_full <- c(b_vec, b_total)
  
  lp_result <- lp("min", cost, A_full, rep("=", length(b_full)), b_full, all.int = TRUE)
  if(lp_result$status != 0) stop("LP did not find an optimal solution.")
  
  y <- lp_result$solution  # LP solution.
  
  # Build assignment pairs: for each pair (a, b), replicate it y times.
  assignment_pairs <- do.call(rbind, lapply(1:N, function(j) {
    count <- y[j]
    if(count > 0) {
      matrix(pair_list[[j]], ncol = 2, nrow = count, byrow = TRUE)
    }
  }))
  
  assignment_pairs <- as.data.frame(assignment_pairs)
  names(assignment_pairs) <- c("grader_id1", "grader_id2")
  
  if(nrow(assignment_pairs) != n_papers) {
    warning("LP matching produced an unexpected number of papers; filling remaining assignments randomly.")
    remaining <- n_papers - nrow(assignment_pairs)
    extra_assignments <- assign_papers(remaining, m, 0, max(loads))
    extra_pairs <- extra_assignments %>%
      group_by(paper_id) %>%
      summarise(grader_id1 = first(grader_id),
                grader_id2 = last(grader_id)) %>%
      ungroup()
    assignment_pairs <- rbind(assignment_pairs, extra_pairs)
  }
  
  assignment_pairs <- assignment_pairs[sample(nrow(assignment_pairs)), ]
  assignment_pairs$paper_id <- 1:nrow(assignment_pairs)
  
  assignments_long <- data.frame(
    paper_id = rep(assignment_pairs$paper_id, each = 2),
    grader_id = c(assignment_pairs$grader_id1, assignment_pairs$grader_id2)
  )
  
  return(assignments_long)
}


# ---------------------------
# 3. Main Simulation Functions for Two Scenarios
# ---------------------------

run_simulation_random <- function(n_papers, n_graders, min_papers, max_papers, rating_probs,
                                  grader_profiles, grader_proportions) {
  papers <- generate_papers(n_papers, rating_probs)
  graders <- generate_graders(n_graders, grader_profiles, grader_proportions)
  assignments <- assign_papers(n_papers, n_graders, min_papers, max_papers)
  sim_data <- simulate_grading(papers, assignments, graders)
  benchmarks <- compute_benchmarks(sim_data)
  list(papers = papers,
       graders = graders,
       assignments = assignments,
       sim_data = sim_data,
       benchmarks = benchmarks)
}

run_simulation_matched <- function(n_papers, n_graders, min_papers, max_papers, rating_probs,
                                   grader_profiles, grader_proportions) {
  papers <- generate_papers(n_papers, rating_probs)
  graders <- generate_graders(n_graders, grader_profiles, grader_proportions)
  assignments <- assign_papers_matched_lp(n_papers, graders, min_papers, max_papers)
  sim_data <- simulate_grading(papers, assignments, graders)
  benchmarks <- compute_benchmarks(sim_data)
  list(papers = papers,
       graders = graders,
       assignments = assignments,
       sim_data = sim_data,
       benchmarks = benchmarks)
}

# ---------------------------
# 4. Run Both Simulations
# ---------------------------
sim_random <- run_simulation_random(n_papers, n_graders, min_papers, max_papers,
                                    rating_probs, grader_profiles, grader_proportions)
sim_matched  <- run_simulation_matched(n_papers, n_graders, min_papers, max_papers,
                                       rating_probs, grader_profiles, grader_proportions)

cat("Random Assignment Benchmarks:\n")
print(sim_random$benchmarks)
cat("\nMatched Assignment Benchmarks:\n")
print(sim_matched$benchmarks)

# ---------------------------
# 5. Bias Estimation & Correction (Same for Both Scenarios)
# ---------------------------
estimate_biases <- function(sim_data) {
  est <- feols(simulated_score ~ 1 | paper_id + grader_id, data = sim_data)
  grader_biases <- fixef(est)$grader_id
  grader_biases <- grader_biases - mean(grader_biases)
  sim_data$estimated_grader_bias <- grader_biases[as.character(sim_data$grader_id)]
  sim_data$corrected_score <- sim_data$simulated_score - sim_data$estimated_grader_bias
  sim_data$corrected_score <- round(pmin(pmax(sim_data$corrected_score, 1), 7))
  return(sim_data)
}

sim_random$sim_data <- estimate_biases(sim_random$sim_data)
sim_matched$sim_data  <- estimate_biases(sim_matched$sim_data)

bias_threshold <- 0.075
apply_selective_correction <- function(sim_data, bias_threshold) {
  sim_data <- sim_data %>%
    mutate(correction_term = ifelse(abs(estimated_grader_bias) >= bias_threshold,
                                    - estimated_grader_bias,
                                    0))
  return(sim_data)
}

sim_random$sim_data <- apply_selective_correction(sim_random$sim_data, bias_threshold)
sim_matched$sim_data  <- apply_selective_correction(sim_matched$sim_data, bias_threshold)

aggregate_paper_scores <- function(sim_data) {
  sim_data %>%
    group_by(paper_id, true_rating) %>%
    summarise(
      avg_raw_score = round(mean(simulated_score)),
      avg_correction = mean(correction_term),
      avg_corrected_score = round(mean(simulated_score) + mean(correction_term)),
      .groups = "drop"
    )
}

paper_results_random  <- aggregate_paper_scores(sim_random$sim_data)
paper_results_matched <- aggregate_paper_scores(sim_matched$sim_data)

overall_rate_random_raw <- mean(paper_results_random$avg_raw_score == paper_results_random$true_rating)
overall_rate_random_corr <- mean(paper_results_random$avg_corrected_score == paper_results_random$true_rating)
overall_rate_matched_raw <- mean(paper_results_matched$avg_raw_score == paper_results_matched$true_rating)
overall_rate_matched_corr <- mean(paper_results_matched$avg_corrected_score == paper_results_matched$true_rating)

cat("Random Assignment - Paper-Averaged Classification Rate (Raw):", overall_rate_random_raw, "\n")
cat("Random Assignment - Paper-Averaged Classification Rate (Corrected):", overall_rate_random_corr, "\n")
cat("Matched Assignment - Paper-Averaged Classification Rate (Raw):", overall_rate_matched_raw, "\n")
cat("Matched Assignment - Paper-Averaged Classification Rate (Corrected):", overall_rate_matched_corr, "\n")

# ---------------------------
# 6. Compute Macro F1 Score (Multinomial)
# ---------------------------
compute_multinomial_F1 <- function(true, pred) {
  classes <- sort(unique(true))
  f1s <- c()
  for (c in classes) {
    tp <- sum(true == c & pred == c)
    fp <- sum(true != c & pred == c)
    fn <- sum(true == c & pred != c)
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    f1s <- c(f1s, f1)
  }
  macro_f1 <- mean(f1s)
  list(macro_f1 = macro_f1, per_class_f1 = f1s)
}

F1_random_raw  <- compute_multinomial_F1(paper_results_random$true_rating, paper_results_random$avg_raw_score)
F1_random_corr <- compute_multinomial_F1(paper_results_random$true_rating, paper_results_random$avg_corrected_score)
F1_matched_raw  <- compute_multinomial_F1(paper_results_matched$true_rating, paper_results_matched$avg_raw_score)
F1_matched_corr <- compute_multinomial_F1(paper_results_matched$true_rating, paper_results_matched$avg_corrected_score)

cat("Random Assignment Macro F1 (Raw):", F1_random_raw$macro_f1, "\n")
cat("Random Assignment Macro F1 (Corrected):", F1_random_corr$macro_f1, "\n")
cat("Matched Assignment Macro F1 (Raw):", F1_matched_raw$macro_f1, "\n")
cat("Matched Assignment Macro F1 (Corrected):", F1_matched_corr$macro_f1, "\n")

# ---------------------------
# 7. (Optional) Visualize with Heatmaps at Paper Level
# ---------------------------
plot_heatmap <- function(paper_results, score_var, title_label) {
  score_counts <- paper_results %>%
    group_by(true_rating, !!sym(score_var)) %>%
    summarise(count = n(), .groups = "drop")
  all_combinations <- expand.grid(true_rating = 1:7, score = 1:7)
  names(all_combinations)[2] <- score_var
  score_counts <- left_join(all_combinations, score_counts, by = c("true_rating", score_var))
  score_counts$count[is.na(score_counts$count)] <- 0
  
  ggplot(score_counts, aes(x = factor(true_rating), y = factor(!!sym(score_var)), fill = count)) +
    geom_tile(color = "grey") +
    geom_text(aes(label = count), size = 4) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(x = "True Rating", y = paste("Average", score_var),
         title = title_label) +
    theme_minimal()
}

heatmap_random_raw  <- plot_heatmap(paper_results_random, "avg_raw_score", 
                                    "Heatmap: True Ratings vs Average Raw Scores (Random)")
heatmap_random_corr <- plot_heatmap(paper_results_random, "avg_corrected_score", 
                                    "Heatmap: True Ratings vs Average Corrected Scores (Random)")
heatmap_matched_raw  <- plot_heatmap(paper_results_matched, "avg_raw_score", 
                                     "Heatmap: True Ratings vs Average Raw Scores (Matched)")
heatmap_matched_corr <- plot_heatmap(paper_results_matched, "avg_corrected_score", 
                                     "Heatmap: True Ratings vs Average Corrected Scores (Matched)")

print(heatmap_random_raw)
print(heatmap_random_corr)
print(heatmap_matched_raw)
print(heatmap_matched_corr)
