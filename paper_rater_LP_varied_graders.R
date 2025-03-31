
#########################
# LOAD REQUIRED PACKAGES
#########################
library(tidyverse)
library(dplyr)
library(fixest)
library(scales)
library(reticulate)
library(glmnet)

# Import Python's PuLP package
pulp <- import("pulp")

#########################
# 1. SIMULATION PARAMETERS
#########################
n_papers <- 3500
min_papers <- 30
max_papers <- 150

num_readers <- 3   # set to desired number of graders per paper


# Matching threshold for grader loads
diff_threshold <- 25

#Bias correction threshold
bias_threshold<- 0.075

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
# (This remains as before; you may wish to adjust it.)
n_graders <- ceiling(2 * (1.5 * n_papers) / max_papers)
if(n_graders * min_papers > 2 * n_papers) {
  stop("Parameters violate minimum load constraints.")
}
cat("n_graders =", n_graders, "\n")

#########################
# 2. HELPER FUNCTIONS
#########################

# (a) Generate paper true ratings.
generate_papers <- function(n_papers, rating_probs) {
  true_ratings <- sample(1:7, size = n_papers, replace = TRUE, prob = rating_probs)
  data.frame(paper_id = 1:n_papers, true_rating = true_ratings)
}

# (b) Generate grader profiles.
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

# (c) Random assignment of papers to graders.
assign_papers <- function(n_papers, n_graders, min_papers, max_papers, num_readers = 2) {
  total_assignments <- num_readers * n_papers
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
    paper_id = rep(1:n_papers, each = num_readers),
    grader_id = grader_vector
  )
  
  # Ensure each paper has num_readers unique graders.
  max_attempts <- 100
  attempt <- 1
  while(any(sapply(split(assignments$grader_id, assignments$paper_id), 
                   function(x) length(unique(x)) != num_readers)) && attempt <= max_attempts) {
    dup_papers <- which(sapply(split(assignments$grader_id, assignments$paper_id), 
                               function(x) length(unique(x)) != num_readers))
    for(i in dup_papers) {
      assignments$grader_id[assignments$paper_id == i] <- sample(1:m, num_readers)
    }
    attempt <- attempt + 1
  }
  if(any(sapply(split(assignments$grader_id, assignments$paper_id), 
                function(x) length(unique(x)) != num_readers)))
    stop("Some papers have duplicate graders after multiple attempts.")
  return(assignments)
}

# (d) Simulate paper grading.
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

# (e) LP solver for the pair-based (num_readers == 2) case.
solve_integer_lp_with_pulp <- function(cost, A, b, time_limit = 600) {
  py_cost <- r_to_py(as.numeric(cost))
  py_A <- r_to_py(as.matrix(A))
  py_b <- r_to_py(as.numeric(b))
  
  py_run_string(sprintf("
import pulp

def solve_integer_lp(cost, A, b, time_limit):
    num_vars = len(cost)
    num_constraints = len(A)
    prob = pulp.LpProblem('LP', pulp.LpMinimize)
    y_vars = [pulp.LpVariable('y_{}'.format(j), lowBound=0, cat='Integer') for j in range(num_vars)]
    prob += pulp.lpSum([cost[j] * y_vars[j] for j in range(num_vars)])
    for i in range(num_constraints):
        prob += pulp.lpSum([A[i][j] * y_vars[j] for j in range(num_vars)]) == b[i]
    solver = pulp.PULP_CBC_CMD(timeLimit=%d, gapRel=0.075, maxNodes=1000000)
    prob.solve(solver)
    status = pulp.LpStatus[prob.status]
    if status not in ['Optimal', 'Not Solved']:
        return None
    solution = [pulp.value(var) for var in y_vars]
    return solution
", time_limit))
  
  solution <- py$solve_integer_lp(py_cost, py_A, py_b, time_limit)
  if(is.null(solution))
    stop("Python LP solver did not return a feasible solution within the time limit.")
  return(as.numeric(solution))
}

# (f) LP solver for the general (num_readers > 2) case.
solve_integer_lp_with_pulp_general <- function(candidate_groups, loads, n_papers, time_limit = 600) {
  m <- length(loads)
  N <- length(candidate_groups)
  cost <- sapply(candidate_groups, function(g) { max(loads[g]) - min(loads[g]) })
  
  A <- matrix(0, nrow = m, ncol = N)
  for(j in 1:N) {
    group <- candidate_groups[[j]]
    A[group, j] <- 1
  }
  A_total <- rbind(A, rep(1, N))
  b_total <- c(loads, n_papers)
  
  py_cost <- r_to_py(as.numeric(cost))
  py_A <- r_to_py(as.matrix(A_total))
  py_b <- r_to_py(as.numeric(b_total))
  
  py_run_string(sprintf("
import pulp

def solve_integer_lp_general(cost, A, b, time_limit):
    num_vars = len(cost)
    num_constraints = len(b)
    prob = pulp.LpProblem('LP', pulp.LpMinimize)
    y_vars = [pulp.LpVariable('y_{}'.format(j), lowBound=0, cat='Integer') for j in range(num_vars)]
    prob += pulp.lpSum([cost[j] * y_vars[j] for j in range(num_vars)])
    for i in range(num_constraints):
        prob += pulp.lpSum([A[i][j] * y_vars[j] for j in range(num_vars)]) == b[i]
    solver = pulp.PULP_CBC_CMD(timeLimit=%d, gapRel=0.075, maxNodes=1000000)
    prob.solve(solver)
    status = pulp.LpStatus[prob.status]
    if status not in ['Optimal', 'Not Solved']:
        return None
    solution = [pulp.value(var) for var in y_vars]
    return solution
", time_limit))
  
  solution <- py$solve_integer_lp_general(py_cost, py_A, py_b, time_limit)
  if(is.null(solution))
    stop("General Python LP solver did not return a feasible solution within the time limit.")
  return(as.numeric(solution))
}


# --- Updated compute_loads: note no self-referential default for num_readers ---
compute_loads <- function(n_papers, n_graders, min_papers, max_papers, num_readers) {
  total_assignments <- num_readers * n_papers
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

# --- Updated LP-based matched assignment function ---
assign_papers_matched_approx <- function(n_papers, graders, min_papers, max_papers, 
                                         diff_threshold = diff_threshold, num_readers, time_limit = 600) {
  m <- nrow(graders)
  # Pass num_readers into compute_loads:
  loads <- compute_loads(n_papers, m, min_papers, max_papers, num_readers = num_readers)
  if(sum(loads) != n_papers * num_readers) stop("Load mismatch.")
  
  if(num_readers == 2) {
    # Use the original pair-based formulation.
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
    
    y <- tryCatch({
      solve_integer_lp_with_pulp(cost, A_full, b_full, time_limit)
    }, error = function(e) {
      warning("Python LP solver failed in primary attempt. Trying fallback approach.")
      return(NULL)
    })
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
      y <- solve_integer_lp_with_pulp(cost, A_full, b_full, time_limit)
      if(is.null(y))
        stop("Fallback LP did not find an optimal solution.")
    }
    
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
      extra_assignments <- assign_papers(remaining, m, 0, max(loads), num_readers = 2)
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
    
  } else {
    # For num_readers > 2, generate candidate groups.
    candidate_groups <- combn(1:m, num_readers, simplify = FALSE)
    # Filter candidate groups based on load closeness.
    candidate_groups <- Filter(function(g) { max(loads[g]) - min(loads[g]) <= diff_threshold }, candidate_groups)
    if(length(candidate_groups) == 0)
      candidate_groups <- combn(1:m, num_readers, simplify = FALSE)
    
    y <- tryCatch({
      solve_integer_lp_with_pulp_general(candidate_groups, loads, n_papers, time_limit)
    }, error = function(e) {
      warning("General Python LP solver failed. Fallback not implemented for num_readers > 2.")
      return(NULL)
    })
    if(is.null(y))
      stop("No feasible solution found for the general assignment LP.")
    
    # Construct the assignment from the solution.
    assignment_list <- list()
    for(j in 1:length(candidate_groups)) {
      count <- y[j]
      if(count > 0) {
        for(k in 1:count) {
          assignment_list[[length(assignment_list) + 1]] <- candidate_groups[[j]]
        }
      }
    }
    if(length(assignment_list) != n_papers) {
      warning("General LP did not produce exactly n_papers groups; using random assignment for remaining papers.")
      remaining <- n_papers - length(assignment_list)
      for(i in 1:remaining) {
        assignment_list[[length(assignment_list) + 1]] <- sample(1:m, num_readers)
      }
    }
    assignments_long <- do.call(rbind, lapply(1:length(assignment_list), function(i) {
      data.frame(paper_id = i, grader_id = assignment_list[[i]])
    }))
    assignments_long <- tidyr::unnest(assignments_long, cols = c(grader_id))
    return(assignments_long)
  }
}



# Get biases

estimate_biases_ridge <- function(sim_data, lambda = 0.1) {
  # Create a design matrix for grader dummies (without intercept)
  X <- model.matrix(~ 0 + factor(grader_id), data = sim_data)
  y <- sim_data$simulated_score
  
  # Fit ridge regression (alpha = 0 is ridge)
  fit <- glmnet(X, y, alpha = 0, lambda = lambda)
  
  # Extract coefficients (they are in a sparse matrix format)
  coef_ridge <- as.vector(coef(fit))[-1]  # Remove intercept if present
  names(coef_ridge) <- levels(factor(sim_data$grader_id))
  
  # Center the bias estimates (so that they have zero mean)
  grader_bias <- coef_ridge - mean(coef_ridge)
  
  # Map each observation its estimated bias.
  sim_data$estimated_grader_bias <- grader_bias[as.character(sim_data$grader_id)]
  
  # Compute the corrected score.
  sim_data$corrected_score <- sim_data$simulated_score - sim_data$estimated_grader_bias
  sim_data$corrected_score <- round(pmin(pmax(sim_data$corrected_score, 1), 7))
  
  return(sim_data)
}

apply_selective_correction_ridge <- function(sim_data, bias_threshold) {
  # Only adjust when the absolute estimated bias exceeds the threshold.
  sim_data$correction_term <- ifelse(abs(sim_data$estimated_grader_bias) >= bias_threshold,
                                     -sim_data$estimated_grader_bias, 0)
  sim_data$corrected_score <- sim_data$simulated_score + sim_data$correction_term
  sim_data$corrected_score <- round(pmin(pmax(sim_data$corrected_score, 1), 7))
  return(sim_data)
}



















#########################
# Simulation Epochs and Aggregation
#########################
#########################

run_simulation_random <- function(n_papers, n_graders, min_papers, max_papers, rating_probs,
                                  grader_profiles, grader_proportions, num_readers) {
  papers <- generate_papers(n_papers, rating_probs)
  graders <- generate_graders(n_graders, grader_profiles, grader_proportions)
  assignments <- assign_papers(n_papers, n_graders, min_papers, max_papers, num_readers)
  sim_data <- simulate_grading(papers, assignments, graders)
  benchmarks <- compute_benchmarks(sim_data)
  list(papers = papers,
       graders = graders,
       assignments = assignments,
       sim_data = sim_data,
       benchmarks = benchmarks)
}

run_simulation_matched <- function(n_papers, n_graders, min_papers, max_papers, rating_probs,
                                   grader_profiles, grader_proportions, num_readers, time_limit = 600) {
  papers <- generate_papers(n_papers, rating_probs)
  graders <- generate_graders(n_graders, grader_profiles, grader_proportions)
  assignments <- assign_papers_matched_approx(n_papers, graders, min_papers, max_papers,
                                              diff_threshold = diff_threshold, num_readers, time_limit = time_limit)
  sim_data <- simulate_grading(papers, assignments, graders)
  benchmarks <- compute_benchmarks(sim_data)
  list(papers = papers,
       graders = graders,
       assignments = assignments,
       sim_data = sim_data,
       benchmarks = benchmarks)
}

#########################
# Epoch-based Simulation and Aggregation (using global num_readers)
#########################
#########################
# Epoch-based Simulation and Aggregation (using global num_readers)
#########################
num_epochs <- 10

# Initialize a list to store epoch-level results.
epoch_results_list <- vector("list", num_epochs)

for (epoch in 1:num_epochs) {
  
  epoch_results <- data.frame()
  
  for (good_prop in good_prop_values) {
    
    # Scale non-good proportions.
    non_good_total <- 1 - good_prop
    new_non_good <- default_non_good / default_non_good_sum * non_good_total
    new_grader_proportions <- c(good = good_prop, new_non_good)
    names(new_grader_proportions)[-1] <- names(default_non_good)
    
    # Run both simulation types using the global num_readers.
    sim_rand <- run_simulation_random(n_papers, n_graders, min_papers, max_papers,
                                      rating_probs, grader_profiles, new_grader_proportions, 
                                      num_readers)
    sim_match <- run_simulation_matched(n_papers, n_graders, min_papers, max_papers,
                                        rating_probs, grader_profiles, new_grader_proportions, 
                                        num_readers, time_limit = 600)
    
    # Apply bias estimation:
    # For the random simulation, use the original fixed-effects correction.
    sim_rand$sim_data <- estimate_biases(sim_rand$sim_data)
    sim_rand$sim_data <- apply_selective_correction(sim_rand$sim_data, bias_threshold)
    # For the LP-based (matched) simulation:
    sim_match$sim_data <- estimate_biases_ridge(sim_match$sim_data, lambda = 0.1)
    sim_match$sim_data <- apply_selective_correction_ridge(sim_match$sim_data, bias_threshold)
    
    
    
    
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
    
    # Combine metrics for this good_prop.
    epoch_results <- rbind(epoch_results,
                           data.frame(epoch = epoch, good_prop = good_prop, method = "Random", version = "Raw",
                                      classification_rate = mean(paper_rand$avg_raw_score == paper_rand$true_rating),
                                      macro_f1 = macro_f1_rand_raw),
                           data.frame(epoch = epoch, good_prop = good_prop, method = "Random", version = "Corrected",
                                      classification_rate = cr_rand,
                                      macro_f1 = F1_rand),
                           data.frame(epoch = epoch, good_prop = good_prop, method = "Matched", version = "Raw",
                                      classification_rate = mean(paper_match$avg_raw_score == paper_match$true_rating),
                                      macro_f1 = macro_f1_match_raw),
                           data.frame(epoch = epoch, good_prop = good_prop, method = "Matched", version = "Corrected",
                                      classification_rate = cr_match,
                                      macro_f1 = F1_match)
    )
  }
  
  epoch_results_list[[epoch]] <- epoch_results
  cat("Epoch", epoch, "completed\n")
}

# Combine results from all epochs.
all_results_df <- do.call(rbind, epoch_results_list)

# Aggregate over epochs to compute means and 95% confidence intervals.
summary_results <- all_results_df %>%
  group_by(good_prop, method, version) %>%
  summarise(
    classification_rate_mean = mean(classification_rate, na.rm = TRUE),
    classification_rate_se   = sd(classification_rate, na.rm = TRUE) / sqrt(n()),
    macro_f1_mean            = mean(macro_f1, na.rm = TRUE),
    macro_f1_se              = sd(macro_f1, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

print(summary_results)

#########################
# Visualization
#########################
library(ggplot2)

# (a) Plot for classification rate with 95% CI.
p_cr <- ggplot(summary_results, aes(x = good_prop, y = classification_rate_mean, color = method, linetype = version)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = classification_rate_mean - 1.96 * classification_rate_se,
                  ymax = classification_rate_mean + 1.96 * classification_rate_se,
                  fill = method), alpha = 0.2, color = NA) +
  labs(title = "Paper-Averaged Classification Rate vs. Good Grader Proportion",
       x = "Proportion of Good Graders", y = "Classification Rate") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal()

# (b) Plot for Macro F1 score with 95% CI.
p_F1 <- ggplot(summary_results, aes(x = good_prop, y = macro_f1_mean, color = method, linetype = version)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = macro_f1_mean - 1.96 * macro_f1_se,
                  ymax = macro_f1_mean + 1.96 * macro_f1_se,
                  fill = method), alpha = 0.2, color = NA) +
  labs(title = "Macro F1 Score vs. Good Grader Proportion",
       x = "Proportion of Good Graders", y = "Macro F1 Score") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal()

print(p_cr)
print(p_F1)

# (c) Bar charts to visualize the improvement due to bias correction.
improvement_df <- summary_results %>%
  select(good_prop, method, version, classification_rate_mean, macro_f1_mean) %>%
  pivot_wider(names_from = version, values_from = c(classification_rate_mean, macro_f1_mean)) %>%
  mutate(
    classification_improvement = classification_rate_mean_Corrected - classification_rate_mean_Raw,
    macro_f1_improvement = macro_f1_mean_Corrected - macro_f1_mean_Raw
  )

# Convert good_prop to a factor for discrete x-axis.
improvement_df$good_prop <- factor(improvement_df$good_prop, levels = sort(unique(improvement_df$good_prop)))

p_improve_class_bar <- ggplot(improvement_df, aes(x = good_prop, y = classification_improvement, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Improvement in Classification Rate with Bias Correction",
       x = "Proportion of Good Graders", y = "Improvement (Corrected - Raw)") +
  scale_x_discrete(labels = function(x) scales::percent(as.numeric(x))) +
  theme_minimal()

p_improve_f1_bar <- ggplot(improvement_df, aes(x = good_prop, y = macro_f1_improvement, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Improvement in Macro F1 Score with Bias Correction",
       x = "Proportion of Good Graders", y = "Improvement (Corrected - Raw)") +
  scale_x_discrete(labels = function(x) scales::percent(as.numeric(x))) +
  theme_minimal()

print(p_improve_class_bar)
print(p_improve_f1_bar)



