library(tidyverse)
library(purrr)
library(tidyr)
library(fixest)
library(scales)
library(reticulate)
library(glmnet)
library(data.table)


# Import Python's PuLP package
pulp <- import("pulp")

### ---------------------------
### Simulation Parameters
### ---------------------------
n_papers <- 3500
min_papers <- 30
max_papers <- 150
num_readers <- 3   # Number of graders per paper

diff_threshold <- 25
bias_threshold <- 0.075

grader_profiles <- list(
  good = list(dist_type = "normal", params = list(bias = 0, sd = 0.5)),
  overly_positive = list(dist_type = "normal", params = list(bias = 1, sd = 0.5)),
  overly_negative = list(dist_type = "normal", params = list(bias = -1, sd = 0.5)),
  apathy = list(dist_type = "uniform", params = list(min = -3, max = 3)),
  erratic = list(dist_type = "normal", params = list(bias = 0, sd = 2)),
  bad = list(dist_type = "gamma", params = list(shape = 2, scale = 1))
)

default_non_good <- c(overly_positive = 0.15,
                      overly_negative = 0.15,
                      apathy = 0.10,
                      erratic = 0.10,
                      bad = 0.10)
default_non_good_sum <- sum(default_non_good)

good_prop_values <- c(0.80, 0.70, 0.60, 0.50, 0.40)
rating_probs <- c(0.3, 0.2, 0.15, 0.12, 0.1, 0.08, 0.05)

n_graders <- ceiling(2 * (1.5 * n_papers) / max_papers)
if(n_graders * min_papers > 2 * n_papers) {
  stop("Parameters violate minimum load constraints.")
}
cat("n_graders =", n_graders, "\n")

### ---------------------------
### Helper Functions
### ---------------------------

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
assign_papers<- function(n_papers, n_graders, min_papers, max_papers, num_readers) {
  total_assignments <- num_readers * n_papers
  m <- n_graders
  
  # Check feasibility.
  if(total_assignments < m * min_papers || total_assignments > m * max_papers)
    stop("No feasible assignment with these parameters.")
  
  # Compute loads in a vectorized manner.
  base_load <- rep(min_papers, m)
  extra_needed <- total_assignments - sum(base_load)
  max_extra <- max_papers - min_papers
  
  # Each grader can get up to 'max_extra' extra assignments.
  extra_pool <- rep(1:m, each = max_extra)
  extra_assignments <- sample(extra_pool, extra_needed, replace = FALSE)
  extra_counts <- tabulate(extra_assignments, nbins = m)
  
  # Final loads: each grader gets the base plus their extra assignments.
  loads <- base_load + extra_counts
  
  # Prepare an empty assignment matrix (rows: papers, columns: reviewer slots).
  assignment_mat <- matrix(NA_integer_, nrow = n_papers, ncol = num_readers)
  
  # Vector of remaining capacities for each grader.
  remaining_cap <- loads
  
  # Process papers in a random order.
  paper_order <- sample(1:n_papers)
  for (paper in paper_order) {
    # Process reviewer slots in a random order.
    slot_order <- sample(1:num_readers)
    for (slot in slot_order) {
      # Determine candidate graders:
      # They must have remaining capacity and not already be assigned to this paper.
      already_assigned <- assignment_mat[paper, ]
      candidates <- which(remaining_cap > 0 & !(1:m %in% already_assigned))
      
      if (length(candidates) == 0) {
        stop("No candidate graders available for paper ", paper)
      }
      
      # Randomly choose one candidate.
      chosen <- sample(candidates, 1)
      assignment_mat[paper, slot] <- chosen
      remaining_cap[chosen] <- remaining_cap[chosen] - 1
    }
  }
  
  # Convert the assignment matrix to a long-format data frame.
  assignments <- data.frame(
    paper_id = rep(1:n_papers, each = num_readers),
    grader_id = as.vector(t(assignment_mat))
  )
  
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
  # Convert inputs to data.table (if not already)
  DT_papers <- as.data.table(papers)
  DT_assignments <- as.data.table(assignments)
  DT_graders <- as.data.table(graders)
  
  # Merge assignments with paper ratings and grader characteristics.
  sim_data <- merge(DT_assignments, DT_papers, by = "paper_id")
  sim_data <- merge(sim_data, DT_graders, by = "grader_id")
  
  # Process simulation for "normal" distribution
  sim_data[dist_type == "normal", `:=`(
    bias = sapply(dist_params, function(x) x$bias),
    sd_val = sapply(dist_params, function(x) x$sd)
  )]
  sim_data[dist_type == "normal", simulated_score := {
    score <- true_rating + rnorm(.N, mean = bias, sd = sd_val)
    score <- round(score)
    pmin(pmax(score, 1), 7)
  }]
  
  # Process simulation for "gamma" distribution
  sim_data[dist_type == "gamma", `:=`(
    shape_val = sapply(dist_params, function(x) x$shape),
    scale_val = sapply(dist_params, function(x) x$scale)
  )]
  sim_data[dist_type == "gamma", simulated_score := {
    score <- true_rating + (rgamma(.N, shape = shape_val, scale = scale_val) - (shape_val * scale_val))
    score <- round(score)
    pmin(pmax(score, 1), 7)
  }]
  
  # Process simulation for "uniform" distribution
  sim_data[dist_type == "uniform", `:=`(
    min_val = sapply(dist_params, function(x) x$min),
    max_val = sapply(dist_params, function(x) x$max)
  )]
  sim_data[dist_type == "uniform", simulated_score := {
    score <- true_rating + runif(.N, min = min_val, max = max_val)
    score <- round(score)
    pmin(pmax(score, 1), 7)
  }]
  
  # Optionally, remove temporary columns created for simulation.
  sim_data[, c("bias", "sd_val", "shape_val", "scale_val", "min_val", "max_val") := NULL]
  
  # Return a data.frame if needed
  return(as.data.frame(sim_data))
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

# (e) LP solvers (kept as in your original code)
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

# (f) Compute loads.
compute_loads <- function(n_papers, n_graders, min_papers, max_papers, num_readers) {
  total_assignments <- num_readers * n_papers
  base_load <- rep(min_papers, n_graders)
  extra_needed <- total_assignments - sum(base_load)
  max_extra <- max_papers - min_papers
  
  # Create a pool where each grader can receive up to max_extra extra assignments
  extra_pool <- rep(1:n_graders, each = max_extra)
  
  # Randomly select extra_needed assignments from the pool without replacement
  extra_assignment <- sample(extra_pool, extra_needed, replace = FALSE)
  
  # Count how many extra assignments each grader gets
  extra_counts <- tabulate(extra_assignment, nbins = n_graders)
  
  # Final loads for each grader
  loads <- base_load + extra_counts
  return(loads)
}

# (g) LP-based matched assignment function.
assign_papers_matched_approx <- function(n_papers, graders, min_papers, max_papers, 
                                         diff_threshold = diff_threshold, num_readers, time_limit = 600) {
  m <- nrow(graders)
  loads <- compute_loads(n_papers, m, min_papers, max_papers, num_readers)
  if(sum(loads) != n_papers * num_readers) stop("Load mismatch.")
  
  if(num_readers == 2) {
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
      warning("Primary LP failed. Trying fallback.")
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
      if(is.null(y)) stop("Fallback LP did not produce a solution.")
    }
    
    assignment_pairs <- do.call(rbind, lapply(1:N, function(j) {
      count <- y[j]
      if(count > 0) matrix(pair_list[[j]], ncol = 2, nrow = count, byrow = TRUE)
    }))
    assignment_pairs <- as.data.frame(assignment_pairs)
    names(assignment_pairs) <- c("grader_id1", "grader_id2")
    if(nrow(assignment_pairs) != n_papers) {
      warning("Unexpected number of assignments; filling randomly.")
      remaining <- n_papers - nrow(assignment_pairs)
      extra_assignments <- assign_papers(remaining, m, 0, max(loads), num_readers = 2)
      extra_pairs <- extra_assignments %>% 
        group_by(paper_id) %>% 
        summarise(grader_id1 = first(grader_id), grader_id2 = last(grader_id)) %>% 
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
    candidate_groups <- combn(1:m, num_readers, simplify = FALSE)
    candidate_groups <- Filter(function(g) { max(loads[g]) - min(loads[g]) <= diff_threshold }, candidate_groups)
    if(length(candidate_groups) == 0)
      candidate_groups <- combn(1:m, num_readers, simplify = FALSE)
    y <- tryCatch({
      solve_integer_lp_with_pulp_general(candidate_groups, loads, n_papers, time_limit)
    }, error = function(e) {
      warning("General LP failed.")
      return(NULL)
    })
    if(is.null(y)) stop("No feasible solution for general LP.")
    assignment_list <- list()
    for(j in 1:length(candidate_groups)) {
      count <- y[j]
      if(count > 0) {
        for(k in 1:count) assignment_list[[length(assignment_list) + 1]] <- candidate_groups[[j]]
      }
    }
    if(length(assignment_list) != n_papers) {
      warning("General LP did not produce exactly n_papers groups; using random assignment for extras.")
      remaining <- n_papers - length(assignment_list)
      for(i in 1:remaining) assignment_list[[length(assignment_list) + 1]] <- sample(1:m, num_readers)
    }
    assignments_long <- do.call(rbind, lapply(1:length(assignment_list), function(i) {
      data.frame(paper_id = i, grader_id = assignment_list[[i]])
    }))
    assignments_long <- tidyr::unnest(assignments_long, cols = c(grader_id))
    return(assignments_long)
  }
}

# (h) Overlap improvement heuristic with incremental update.
improve_overlap <- function(assignments, n_graders, num_readers, load_targets, 
                            min_overlap_fraction = 0.5, min_papers, max_iter = 1000, tol = 0.01) {
  compute_overlap <- function(assignments, n_graders) {
    overlap <- matrix(0, nrow = n_graders, ncol = n_graders)
    papers <- unique(assignments$paper_id)
    for (p in papers) {
      graders_p <- sort(unique(assignments$grader_id[assignments$paper_id == p]))
      if (length(graders_p) >= 2) {
        for (i in 1:(length(graders_p) - 1)) {
          for (j in (i + 1):length(graders_p)) {
            overlap[graders_p[i], graders_p[j]] <- overlap[graders_p[i], graders_p[j]] + 1
            overlap[graders_p[j], graders_p[i]] <- overlap[graders_p[j], graders_p[i]] + 1
          }
        }
      }
    }
    return(overlap)
  }
  
  overlap <- compute_overlap(assignments, n_graders)
  iter <- 0
  improved <- TRUE
  while (improved && iter < max_iter) {
    iter <- iter + 1
    improved <- FALSE
    deficits <- matrix(0, nrow = n_graders, ncol = n_graders)
    for (j in 1:(n_graders - 1)) {
      for (k in (j + 1):n_graders) {
        req_overlap <- min_overlap_fraction * min(load_targets[j], load_targets[k])
        deficits[j, k] <- req_overlap - overlap[j, k]
      }
    }
    max_deficit <- max(deficits)
    if (max_deficit <= tol) break
    idx <- which(deficits == max_deficit, arr.ind = TRUE)[1, ]
    j <- idx[1]
    k <- idx[2]
    papers_j <- assignments$paper_id[assignments$grader_id == j]
    papers_k <- assignments$paper_id[assignments$grader_id == k]
    candidate_papers <- setdiff(papers_j, papers_k)
    swap_done <- FALSE
    if (length(candidate_papers) > 0) {
      for (p in candidate_papers) {
        graders_p <- sort(unique(assignments$grader_id[assignments$paper_id == p]))
        candidate_i <- setdiff(graders_p, j)
        if (length(candidate_i) == 0) next
        current_loads <- as.numeric(table(factor(assignments$grader_id, levels = 1:n_graders)))
        for (i in candidate_i) {
          if ((current_loads[i] - 1) < min_papers) next
          if (!(k %in% graders_p)) {
            old_graders <- graders_p
            new_graders <- union(setdiff(old_graders, i), k)
            assignments <- assignments[!(assignments$paper_id == p & assignments$grader_id == i), ]
            assignments <- rbind(assignments, data.frame(paper_id = p, grader_id = k))
            for (g in setdiff(old_graders, i)) {
              overlap[i, g] <- overlap[i, g] - 1
              overlap[g, i] <- overlap[g, i] - 1
            }
            for (g in setdiff(new_graders, k)) {
              overlap[k, g] <- overlap[k, g] + 1
              overlap[g, k] <- overlap[g, k] + 1
            }
            swap_done <- TRUE
            improved <- TRUE
            break
          }
        }
        if (swap_done) break
      }
    }
  }
  cat("Overlap improvement heuristic finished in", iter, "iterations.\n")
  return(assignments)
}

# (i) Bias estimation and correction functions.
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

estimate_biases_ridge <- function(sim_data, lambda = 0.1) {
  X <- model.matrix(~ 0 + factor(grader_id), data = sim_data)
  y <- sim_data$simulated_score
  fit <- glmnet(X, y, alpha = 0, lambda = lambda)
  coef_ridge <- as.vector(coef(fit))[-1]
  names(coef_ridge) <- levels(factor(sim_data$grader_id))
  grader_bias <- coef_ridge - mean(coef_ridge)
  sim_data$estimated_grader_bias <- grader_bias[as.character(sim_data$grader_id)]
  sim_data$corrected_score <- sim_data$simulated_score - sim_data$estimated_grader_bias
  sim_data$corrected_score <- round(pmin(pmax(sim_data$corrected_score, 1), 7))
  return(sim_data)
}

apply_selective_correction_ridge <- function(sim_data, bias_threshold) {
  sim_data$correction_term <- ifelse(abs(sim_data$estimated_grader_bias) >= bias_threshold,
                                     -sim_data$estimated_grader_bias, 0)
  sim_data$corrected_score <- sim_data$simulated_score + sim_data$correction_term
  sim_data$corrected_score <- round(pmin(pmax(sim_data$corrected_score, 1), 7))
  return(sim_data)
}

# (j) Aggregation functions.
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

### ---------------------------
### Simulation Functions
### ---------------------------

# run_simulation_dual: returns a list with paper-level outputs.
run_simulation_dual <- function(n_papers, n_graders, min_papers, max_papers, rating_probs,
                                grader_profiles, grader_proportions, num_readers, time_limit = 600) {
  papers <- generate_papers(n_papers, rating_probs)
  graders <- generate_graders(n_graders, grader_profiles, grader_proportions)
  assignments_random <- assign_papers(n_papers, n_graders, min_papers, max_papers, num_readers)
  assignments_lp_initial <- assign_papers_matched_approx(n_papers, graders, min_papers, max_papers,
                                                         diff_threshold = diff_threshold, num_readers, time_limit = time_limit)
  load_targets <- compute_loads(n_papers, n_graders, min_papers, max_papers, num_readers)
  assignments_lp <- improve_overlap(assignments_lp_initial, n_graders, num_readers, load_targets, 
                                    min_overlap_fraction = 0.5, min_papers = min_papers, max_iter = 1000)
  sim_data_rand <- simulate_grading(papers, assignments_random, graders)
  sim_data_lp <- simulate_grading(papers, assignments_lp, graders)
  
  list(
    papers = papers,
    graders = graders,
    assignments_random = assignments_random,
    assignments_lp = assignments_lp,
    sim_data_rand = sim_data_rand,
    sim_data_lp = sim_data_lp
  )
}

# run_simulation_for_good_prop: as defined earlier.
run_simulation_for_good_prop <- function(good_prop, 
                                         n_papers, n_graders, min_papers, max_papers,
                                         rating_probs, grader_profiles, default_non_good, 
                                         default_non_good_sum, num_readers, bias_threshold, time_limit) {
  non_good_total <- 1 - good_prop
  new_non_good <- default_non_good / default_non_good_sum * non_good_total
  new_grader_proportions <- c(good = good_prop, new_non_good)
  names(new_grader_proportions)[-1] <- names(default_non_good)
  
  sim_dual <- run_simulation_dual(n_papers, n_graders, min_papers, max_papers,
                                  rating_probs, grader_profiles, new_grader_proportions, 
                                  num_readers, time_limit)
  
  sim_dual$sim_data_rand <- apply_selective_correction(estimate_biases(sim_dual$sim_data_rand), bias_threshold)
  sim_dual$sim_data_lp   <- apply_selective_correction_ridge(estimate_biases_ridge(sim_dual$sim_data_lp, lambda = 0.1), bias_threshold)
  
  paper_rand <- aggregate_paper_scores(sim_dual$sim_data_rand)
  paper_lp   <- aggregate_paper_scores(sim_dual$sim_data_lp)
  
  cr_rand <- mean(paper_rand$avg_corrected_score == paper_rand$true_rating)
  cr_lp   <- mean(paper_lp$avg_corrected_score == paper_lp$true_rating)
  
  macro_f1_rand_raw <- compute_multinomial_F1(paper_rand$true_rating, paper_rand$avg_raw_score)$macro_f1
  macro_f1_lp_raw    <- compute_multinomial_F1(paper_lp$true_rating, paper_lp$avg_raw_score)$macro_f1
  F1_rand <- compute_multinomial_F1(paper_rand$true_rating, paper_rand$avg_corrected_score)$macro_f1
  F1_lp <- compute_multinomial_F1(paper_lp$true_rating, paper_lp$avg_corrected_score)$macro_f1
  
  list(
    metrics = tibble(
      good_prop = good_prop,
      match_type = c("Random", "Random", "LP", "LP"),
      version = c("Raw", "Corrected", "Raw", "Corrected"),
      classification_rate = c(mean(paper_rand$avg_raw_score == paper_rand$true_rating), cr_rand,
                              mean(paper_lp$avg_raw_score == paper_lp$true_rating), cr_lp),
      macro_f1 = c(macro_f1_rand_raw, F1_rand, macro_f1_lp_raw, F1_lp)
    ),
    paper_rand = paper_rand,
    paper_lp = paper_lp
  )
}

# Wrap with possibly().
run_simulation_for_good_prop_possibly <- possibly(
  run_simulation_for_good_prop,
  otherwise = list(
    metrics = tibble(good_prop = NA_real_, match_type = NA_character_, version = NA_character_,
                     classification_rate = NA_real_, macro_f1 = NA_real_),
    paper_rand = tibble(), 
    paper_lp = tibble()
  )
)

# run_one_epoch: runs over all good_prop values in one epoch.
run_one_epoch <- function(epoch_num, good_prop_values, params) {
  results_list <- map(good_prop_values, function(gp) {
    run_simulation_for_good_prop_possibly(
      gp,
      params$n_papers,
      params$n_graders,
      params$min_papers,
      params$max_papers,
      params$rating_probs,
      params$grader_profiles,
      params$default_non_good,
      params$default_non_good_sum,
      params$num_readers,
      params$bias_threshold,
      params$time_limit
    )
  })
  
  metrics_df <- map_dfr(results_list, "metrics") %>% mutate(epoch = epoch_num)
  paper_rand_df <- bind_rows(map(results_list, "paper_rand")) %>% mutate(epoch = epoch_num)
  paper_lp_df   <- bind_rows(map(results_list, "paper_lp")) %>% mutate(epoch = epoch_num)
  
  list(
    metrics = metrics_df,
    paper_rand = paper_rand_df,
    paper_lp = paper_lp_df
  )
}

# ---------------------------
# Define parameter list.
# ---------------------------
params <- list(
  n_papers = n_papers,
  n_graders = n_graders,
  min_papers = min_papers,
  max_papers = max_papers,
  rating_probs = rating_probs,
  grader_profiles = grader_profiles,
  default_non_good = default_non_good,
  default_non_good_sum = default_non_good_sum,
  num_readers = num_readers,
  bias_threshold = bias_threshold,
  time_limit = 600
)














# ---------------------------
# Run simulation over epochs.
# ---------------------------
num_epochs <- 1

epoch_results <- map(1:num_epochs, function(epoch) {
  run_one_epoch(epoch, good_prop_values, params)
})

# Combine aggregated metrics across epochs.
all_metrics_df <- map_dfr(epoch_results, "metrics")
# Combine paper-level outputs (Random matching) for confusion matrix.
all_paper_rand <- bind_rows(map(epoch_results, "paper_rand"))
all_metrics_df<-all_metrics_df[complete.cases(all_metrics_df),]

print(all_metrics_df)

# ---------------------------
# Compute paper_counts for the confusion matrix.
# ---------------------------

# Diagnostic: check how many paper records we have.
total_observed <- nrow(all_paper_rand)
cat("Total observed papers:", total_observed, "\n")

# Compute confusion matrix counts.
paper_counts <- all_paper_rand %>%
  group_by(true_rating, avg_corrected_score) %>%
  summarise(count = n(), .groups = "drop")

# Compute percentage relative to the actual total.
paper_counts <- paper_counts %>%
  mutate(percentage = count / total_observed * 100)

# Ensure every combination (1:7 x 1:7) is present.
all_combs <- expand.grid(true_rating = 1:7, avg_corrected_score = 1:7)
paper_counts <- left_join(all_combs, paper_counts, by = c("true_rating", "avg_corrected_score"))
paper_counts$percentage[is.na(paper_counts$percentage)] <- 0

# Plot the heatmap.
heatmap_conf <- ggplot(paper_counts, aes(x = factor(avg_corrected_score), y = factor(true_rating), fill = percentage)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = round(percentage, 1)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Macro Confusion Matrix (Observed %):\nTrue Ratings vs. Predicted Scores (Random)",
       x = "Predicted (Avg Corrected) Score", y = "True Rating") +
  theme_minimal()

print(heatmap_conf)










#########################
# Visualizations
#########################
# Aggregate over epochs to compute means and 95% confidence intervals.
summary_results <- all_metrics_df %>%
  group_by(good_prop, match_type, version) %>%
  summarise(
    classification_rate_mean = mean(classification_rate, na.rm = TRUE),
    classification_rate_se   = sd(classification_rate, na.rm = TRUE) / sqrt(n()),
    macro_f1_mean            = mean(macro_f1, na.rm = TRUE),
    macro_f1_se              = sd(macro_f1, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# (a) Plot for Classification Rate with 95% CI.
p_cr <- ggplot(summary_results, aes(x = good_prop, y = classification_rate_mean, 
                                    color = match_type, linetype = version)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = classification_rate_mean - 1.96 * classification_rate_se,
                  ymax = classification_rate_mean + 1.96 * classification_rate_se,
                  fill = match_type), alpha = 0.2, color = NA) +
  labs(title = "Paper-Averaged Classification Rate vs. Good Grader Proportion",
       x = "Proportion of Good Graders", y = "Classification Rate") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal()

# (b) Plot for Macro F1 Score with 95% CI.
p_F1 <- ggplot(summary_results, aes(x = good_prop, y = macro_f1_mean, 
                                    color = match_type, linetype = version)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = macro_f1_mean - 1.96 * macro_f1_se,
                  ymax = macro_f1_mean + 1.96 * macro_f1_se,
                  fill = match_type), alpha = 0.2, color = NA) +
  labs(title = "Macro F1 Score vs. Good Grader Proportion",
       x = "Proportion of Good Graders", y = "Macro F1 Score") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal()

print(p_cr)
print(p_F1)

# (c) Bar Charts for Improvement (Corrected - Raw)
# Pivot the summary_results so that "Raw" and "Corrected" are separate columns.
improvement_df <- summary_results %>%
  select(good_prop, match_type, version, classification_rate_mean, macro_f1_mean) %>%
  pivot_wider(names_from = version, values_from = c(classification_rate_mean, macro_f1_mean)) %>%
  mutate(
    classification_improvement = classification_rate_mean_Corrected - classification_rate_mean_Raw,
    macro_f1_improvement = macro_f1_mean_Corrected - macro_f1_mean_Raw
  )

# Convert good_prop to a factor for discrete x-axis.
improvement_df$good_prop <- factor(improvement_df$good_prop, levels = sort(unique(improvement_df$good_prop)))

p_improve_class_bar <- ggplot(improvement_df, aes(x = good_prop, y = classification_improvement, fill = match_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Improvement in Classification Rate with Bias Correction",
       x = "Proportion of Good Graders", y = "Improvement (Corrected - Raw)") +
  scale_x_discrete(labels = function(x) scales::percent(as.numeric(x))) +
  theme_minimal()

p_improve_f1_bar <- ggplot(improvement_df, aes(x = good_prop, y = macro_f1_improvement, fill = match_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Improvement in Macro F1 Score with Bias Correction",
       x = "Proportion of Good Graders", y = "Improvement (Corrected - Raw)") +
  scale_x_discrete(labels = function(x) scales::percent(as.numeric(x))) +
  theme_minimal()

print(p_improve_class_bar)
print(p_improve_f1_bar)


























# Helper: run simulation for a given number of readers.
run_simulation_for_readers <- function(num_readers_val, num_epochs, good_prop_values, params_base) {
  # Update the parameter list with the new num_readers value.
  params <- params_base
  params$num_readers <- num_readers_val
  
  # Run the simulation for each epoch.
  epoch_results <- map(1:num_epochs, function(epoch) {
    run_one_epoch(epoch, good_prop_values, params)
  })
  
  # Combine aggregated metrics from all epochs and add a column for num_readers.
  metrics_df <- map_dfr(epoch_results, "metrics") %>%
    mutate(num_readers = num_readers_val)
  
  return(metrics_df)
}

# Define your base parameters (as in your current script).
params <- list(
  n_papers = n_papers,
  n_graders = n_graders,
  min_papers = min_papers,
  max_papers = max_papers,
  rating_probs = rating_probs,
  grader_profiles = grader_profiles,
  default_non_good = default_non_good,
  default_non_good_sum = default_non_good_sum,
  num_readers = num_readers,  # This will be overridden.
  bias_threshold = bias_threshold,
  time_limit = 600
)

# Set the number of epochs for each simulation.
num_epochs <- 5

# Specify the different numbers of readers you want to test.
readers_to_test <- c(2, 3, 4)

# Run the simulation for each value of num_readers and combine the results.
all_metrics_readers <- map_dfr(readers_to_test, function(nr) {
  run_simulation_for_readers(nr, num_epochs, good_prop_values, params)
})

print(all_metrics_readers)
