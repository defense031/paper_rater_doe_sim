# ============================================================
# Simulation of Grading: Random vs. Approximate LP-based Matched Assignment 
# with Sensitivity Testing for Grader Quality Composition
# ============================================================

# Load required packages
library(tidyverse)
library(dplyr)
library(fixest)
library(scales)
library(reticulate)

# Import Python's PuLP package
pulp <- import("pulp")

# ---------------------------
# 1. Define Simulation Parameters
# ---------------------------
n_papers <- 5000
min_papers <- 30
max_papers <- 150

# Matching threshold for grader loads
diff_threshold <- 15

# Define six grader profiles.
grader_profiles <- list(
  good = list(dist_type = "normal", params = list(bias = 0, sd = 0.5)),
  overly_positive = list(dist_type = "normal", params = list(bias = 0.5, sd = 0.5)),
  overly_negative = list(dist_type = "normal", params = list(bias = -0.5, sd = 0.5)),
  apathy = list(dist_type = "uniform", params = list(min = -0.5, max = 0.5)),
  erratic = list(dist_type = "normal", params = list(bias = 0, sd = 2)),
  bad = list(dist_type = "gamma", params = list(shape = 2, scale = 1))
)

# Default non-good proportions when good = 0.40.
default_non_good <- c(overly_positive = 0.15,
                      overly_negative = 0.15,
                      apathy = 0.10,
                      erratic = 0.10,
                      bad = 0.10)
default_non_good_sum <- sum(default_non_good)  # 0.60

# We'll vary the proportion of good graders: 80%, 70%, 60%, 50%, 40%.
good_prop_values <- c(0.80, 0.70, 0.60, 0.50, 0.40)

# True rating distribution probabilities.
rating_probs <- c(0.3, 0.2, 0.15, 0.12, 0.1, 0.08, 0.05)

# Determine number of graders (using a load-based heuristic).
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
  } else if(dist_type == "uniform") {
    error <- runif(1, min = params$min, max = params$max)
    score <- true_rating + error
  } else {
    stop("Unknown distribution type.")
  }
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
  list(mae = mae, mse = mse, classification_rate = classification_rate)
}

# ----- Define Python LP Solver using PuLP via reticulate -----
solve_integer_lp_with_pulp <- function(cost, A, b, time_limit = 600) {
  # Convert R objects to Python objects
  py_cost <- r_to_py(as.numeric(cost))
  py_A <- r_to_py(as.matrix(A))
  py_b <- r_to_py(as.numeric(b))
  
  # Define and run the Python function using reticulate.
  # We use pulp.PULP_CBC_CMD with a time limit, a relative gap tolerance (7.5%% here), and a maximum node limit.
  py_run_string(sprintf("
import pulp

def solve_integer_lp(cost, A, b, time_limit):
    num_vars = len(cost)
    num_constraints = len(A)
    prob = pulp.LpProblem('LP', pulp.LpMinimize)
    # Create integer variables: y_0, y_1, ... y_(num_vars-1)
    y_vars = [pulp.LpVariable('y_{}'.format(j), lowBound=0, cat='Integer') for j in range(num_vars)]
    # Objective: minimize sum(cost[j] * y_vars[j])
    prob += pulp.lpSum([cost[j] * y_vars[j] for j in range(num_vars)])
    # Add equality constraints: for each row i, sum(A[i][j] * y_vars[j]) == b[i]
    for i in range(num_constraints):
        prob += pulp.lpSum([A[i][j] * y_vars[j] for j in range(num_vars)]) == b[i]
    # Set up the CBC solver with a time limit, a relative gap tolerance (7.5%% here), and a maximum node limit.
    solver = pulp.PULP_CBC_CMD(timeLimit=%d, gapRel=0.075, maxNodes=1000000)
    prob.solve(solver)
    status = pulp.LpStatus[prob.status]
    if status not in ['Optimal', 'Not Solved']:
        return None
    solution = [pulp.value(var) for var in y_vars]
    return solution
", time_limit))
  
  # Call the Python function with the provided time limit.
  solution <- py$solve_integer_lp(py_cost, py_A, py_b, time_limit)
  if (is.null(solution))
    stop("Python LP solver did not return a feasible solution within the time limit.")
  return(as.numeric(solution))
}


# ----- LP-based Approximate Matching Assignment Function with Python LP -----
compute_loads <- function(n_papers, n_graders, min_papers, max_papers) {
  total_assignments <- 2 * n_papers
  loads <- rep(min_papers, n_graders)
  remaining <- total_assignments - sum(loads)
  max_extra <- max_papers - min_papers
  for(i in 1:remaining) {
    candidates <- which(loads - min_papers < max_extra)
    selected <- sample(candidates, 1)
    loads[selected] <- loads[selected] + 1
  }
  return(loads)
}


assign_papers_matched_approx <- function(n_papers, graders, min_papers, max_papers, diff_threshold = diff_threshold) {
  m <- nrow(graders)
  loads <- compute_loads(n_papers, m, min_papers, max_papers)
  if(sum(loads) != 2 * n_papers) stop("Load mismatch.")
  
  # Build candidate pair list: include pairs with load difference <= diff_threshold.
  pair_list <- list()
  pair_index <- 1
  for(a in 1:(m-1)) {
    for(b in (a+1):m) {
      if(abs(loads[a] - loads[b]) <= diff_threshold) {
        pair_list[[pair_index]] <- c(a, b)
        pair_index <- pair_index + 1
      }
    }
  }
  # If no candidate pairs meet the threshold, use all pairs.
  if(length(pair_list) == 0) {
    for(a in 1:(m-1)) {
      for(b in (a+1):m) {
        pair_list[[pair_index]] <- c(a, b)
        pair_index <- pair_index + 1
      }
    }
  }
  
  N <- length(pair_list)
  cost <- sapply(pair_list, function(p) { abs(loads[p[1]] - loads[p[2]]) })
  
  A <- matrix(0, nrow = m, ncol = N)
  for(j in 1:N) {
    pair <- pair_list[[j]]
    A[pair[1], j] <- 1
    A[pair[2], j] <- 1
  }
  b_vec <- loads
  
  A_total <- matrix(1, nrow = 1, ncol = N)
  b_total <- n_papers
  
  A_full <- rbind(A, A_total)
  b_full <- c(b_vec, b_total)
  
  # Solve the LP using the Python solver
  y <- tryCatch({
    solve_integer_lp_with_pulp(cost, A_full, b_full)
  }, error = function(e) {
    warning("Python LP solver failed in primary attempt. Trying fallback approach.")
    return(NULL)
  })
  
  # Fallback: use all candidate pairs if the first attempt fails.
  if(is.null(y)) {
    pair_list <- list()
    pair_index <- 1
    for(a in 1:(m-1)) {
      for(b in (a+1):m) {
        pair_list[[pair_index]] <- c(a, b)
        pair_index <- pair_index + 1
      }
    }
    N <- length(pair_list)
    cost <- sapply(pair_list, function(p) { abs(loads[p[1]] - loads[p[2]]) })
    A <- matrix(0, nrow = m, ncol = N)
    for(j in 1:N) {
      pair <- pair_list[[j]]
      A[pair[1], j] <- 1
      A[pair[2], j] <- 1
    }
    b_vec <- loads
    A_total <- matrix(1, nrow = 1, ncol = N)
    b_total <- n_papers
    A_full <- rbind(A, A_total)
    b_full <- c(b_vec, b_total)
    y <- solve_integer_lp_with_pulp(cost, A_full, b_full)
    if(is.null(y))
      stop("Fallback LP did not find an optimal solution.")
  }
  
  # Process solution vector y to construct assignment pairs.
  assignment_pairs <- do.call(rbind, lapply(1:N, function(j) {
    count <- y[j]
    if(count > 0) {
      matrix(pair_list[[j]], ncol = 2, nrow = count, byrow = TRUE)
    }
  }))
  
  assignment_pairs <- as.data.frame(assignment_pairs)
  names(assignment_pairs) <- c("grader_id1", "grader_id2")
  
  if(nrow(assignment_pairs) != n_papers) {
    warning("Approximate matching produced an unexpected number of papers; filling remaining assignments randomly.")
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
  assignments <- assign_papers_matched_approx(n_papers, graders, min_papers, max_papers, diff_threshold = diff_threshold)
  sim_data <- simulate_grading(papers, assignments, graders)
  benchmarks <- compute_benchmarks(sim_data)
  list(papers = papers,
       graders = graders,
       assignments = assignments,
       sim_data = sim_data,
       benchmarks = benchmarks)
}

# Bias estimation & selective correction functions.
estimate_biases <- function(sim_data) {
  est <- feols(simulated_score ~ 1 | paper_id + grader_id, data = sim_data)
  grader_biases <- fixef(est)$grader_id
  grader_biases <- grader_biases - mean(grader_biases)
  sim_data$estimated_grader_bias <- grader_biases[as.character(sim_data$grader_id)]
  sim_data$corrected_score <- sim_data$simulated_score - sim_data$estimated_grader_bias
  sim_data$corrected_score <- round(pmin(pmax(sim_data$corrected_score, 1), 7))
  return(sim_data)
}

apply_selective_correction <- function(sim_data, bias_threshold) {
  sim_data$correction_term <- ifelse(abs(sim_data$estimated_grader_bias) >= bias_threshold,
                                     - sim_data$estimated_grader_bias,
                                     0)
  return(sim_data)
}

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

# ---------------------------
# Sensitivity Testing: Vary Good Grader Proportion (80% to 40%)
# ---------------------------
sensitivity_results <- data.frame()
bias_threshold <- 0.075

for(good_prop in good_prop_values) {
  # Scale non-good types proportionally.
  non_good_total <- 1 - good_prop
  new_non_good <- default_non_good / default_non_good_sum * non_good_total
  new_grader_proportions <- c(good = good_prop, new_non_good)
  names(new_grader_proportions)[-1] <- names(default_non_good)
  
  # Run both simulations.
  sim_rand <- run_simulation_random(n_papers, n_graders, min_papers, max_papers,
                                    rating_probs, grader_profiles, new_grader_proportions)
  sim_match <- run_simulation_matched(n_papers, n_graders, min_papers, max_papers,
                                      rating_probs, grader_profiles, new_grader_proportions)
  
  # Estimate biases and apply selective correction.
  sim_rand$sim_data <- estimate_biases(sim_rand$sim_data)
  sim_match$sim_data <- estimate_biases(sim_match$sim_data)
  sim_rand$sim_data <- apply_selective_correction(sim_rand$sim_data, bias_threshold)
  sim_match$sim_data <- apply_selective_correction(sim_match$sim_data, bias_threshold)
  
  # Aggregate paper-level scores.
  paper_rand <- aggregate_paper_scores(sim_rand$sim_data)
  paper_match <- aggregate_paper_scores(sim_match$sim_data)
  
  # Compute performance metrics.
  cr_rand <- mean(paper_rand$avg_corrected_score == paper_rand$true_rating)
  cr_match <- mean(paper_match$avg_corrected_score == paper_match$true_rating)
  
  # Compute macro F1 for raw scores.
  macro_f1_rand_raw <- compute_multinomial_F1(paper_rand$true_rating, paper_rand$avg_raw_score)$macro_f1
  if(length(macro_f1_rand_raw) == 0) macro_f1_rand_raw <- NA
  
  macro_f1_match_raw <- compute_multinomial_F1(paper_match$true_rating, paper_match$avg_raw_score)$macro_f1
  if(length(macro_f1_match_raw) == 0) macro_f1_match_raw <- NA
  
  # Compute macro F1 for corrected scores.
  F1_rand <- compute_multinomial_F1(paper_rand$true_rating, paper_rand$avg_corrected_score)$macro_f1
  F1_match <- compute_multinomial_F1(paper_match$true_rating, paper_match$avg_corrected_score)$macro_f1
  
  sensitivity_results <- rbind(sensitivity_results,
                               data.frame(good_prop = good_prop, method = "Random", version = "Raw",
                                          classification_rate = mean(paper_rand$avg_raw_score == paper_rand$true_rating),
                                          macro_f1 = macro_f1_rand_raw),
                               data.frame(good_prop = good_prop, method = "Random", version = "Corrected",
                                          classification_rate = cr_rand, macro_f1 = F1_rand),
                               data.frame(good_prop = good_prop, method = "Matched", version = "Raw",
                                          classification_rate = mean(paper_match$avg_raw_score == paper_match$true_rating),
                                          macro_f1 = macro_f1_match_raw),
                               data.frame(good_prop = good_prop, method = "Matched", version = "Corrected",
                                          classification_rate = cr_match, macro_f1 = F1_match)
  )
}


print(sensitivity_results)

# ---------------------------
# Visualization: 4 Lines per Graph
# ---------------------------
p_cr <- ggplot(sensitivity_results, aes(x = good_prop, y = classification_rate, color = method, linetype = version)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(title = "Paper-Averaged Classification Rate vs. Good Grader Proportion",
       x = "Proportion of Good Graders", y = "Classification Rate",
       color = "Method", linetype = "Score Version") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_linetype_manual(values = c("Raw" = "solid", "Corrected" = "dashed")) +
  theme_minimal()

p_F1 <- ggplot(sensitivity_results, aes(x = good_prop, y = macro_f1, color = method, linetype = version)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(title = "Macro F1 Score vs. Good Grader Proportion",
       x = "Proportion of Good Graders", y = "Macro F1 Score",
       color = "Method", linetype = "Score Version") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_linetype_manual(values = c("Raw" = "solid", "Corrected" = "dashed")) +
  theme_minimal()


print(p_cr)
print(p_F1)

# ---------------------------
# (Optional) Heatmaps for one selected scenario (e.g., when good_prop = 0.60, Random Corrected)
# ---------------------------
non_good_total <- 1 - 0.60
new_non_good <- default_non_good / default_non_good_sum * non_good_total
new_grader_proportions <- c(good = 0.60, new_non_good)
names(new_grader_proportions)[-1] <- names(default_non_good)

sim_demo <- run_simulation_random(n_papers, n_graders, min_papers, max_papers,
                                  rating_probs, grader_profiles, new_grader_proportions)
sim_demo$sim_data <- estimate_biases(sim_demo$sim_data)
sim_demo$sim_data <- apply_selective_correction(sim_demo$sim_data, bias_threshold)
paper_demo <- aggregate_paper_scores(sim_demo$sim_data)
# Pre-aggregate counts for the heatmap.
paper_counts <- paper_demo %>%
  group_by(true_rating, avg_corrected_score) %>%
  summarise(count = n(), .groups = "drop")

# Ensure that every combination of true_rating and avg_corrected_score (1 through 7) is present.
all_combs <- expand.grid(true_rating = 1:7, avg_corrected_score = 1:7)
paper_counts <- left_join(all_combs, paper_counts, 
                          by = c("true_rating", "avg_corrected_score"))
paper_counts$count[is.na(paper_counts$count)] <- 0

# Plot the heatmap using geom_tile.
heatmap_demo <- ggplot(paper_counts, aes(x = factor(true_rating), y = factor(avg_corrected_score), fill = count)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = count), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Heatmap: True Ratings vs Average Corrected Scores (Random, 60% Good)",
       x = "True Rating", y = "Average Corrected Score") +
  theme_minimal()

print(heatmap_demo)
