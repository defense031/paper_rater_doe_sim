# This script will create a simulation that applies a DOE concept to test the final screening under multiple 
## paper to rater matching methods

library(tidyverse)
library(dplyr)

# Simulation Shell for Grading 

# Simulation Shell for Grading Assignments in R with Modular Grader Profiles

# Load necessary libraries
library(dplyr)

# ---------------------------
# 1. Define Simulation Parameters
# ---------------------------

# Number of papers to simulate
n_papers <- 1000

# Grader load constraints
min_papers <- 30
max_papers <- 150

# Define grader profiles.
# Each element is a list with two elements:
#  - dist_type: the distribution used to simulate grading error.
#  - params: a list of parameters for that distribution.
#
# In this proof-of-concept we have two types:
# For example, "good" graders use a normal distribution with bias 0 and sd 0.5;
# "bad" graders use a gamma distribution with shape 2 and scale 0.5.
grader_profiles <- list(
  good = list(dist_type = "normal", params = list(bias = 0, sd = 0.5)),
  bad = list(dist_type = "gamma", params = list(shape = 2, scale = 0.5))
)

# Proportions for each grader type (order should match names(grader_profiles))
grader_proportions <- c(good = 0.9, bad = 0.1)

# True rating distribution probabilities for ratings 1 (most likely) to 7 (least likely)
rating_probs <- c(0.3, 0.2, 0.15, 0.12, 0.1, 0.08, 0.05)

# Determine number of graders.
# For feasibility, we require: m * min_papers <= 2*n_papers <= m * max_papers.
n_graders <- ceiling(2 * n_papers / max_papers)
if(n_graders * min_papers > 2 * n_papers) {
  stop("The chosen parameters do not allow a feasible assignment (minimum load constraint violated).")
}

# ---------------------------
# 2. Define Helper Functions
# ---------------------------

# Function to generate papers with a "true" rating
generate_papers <- function(n_papers, rating_probs) {
  true_ratings <- sample(1:7, size = n_papers, replace = TRUE, prob = rating_probs)
  papers <- data.frame(paper_id = 1:n_papers, true_rating = true_ratings)
  return(papers)
}

# Function to generate grader profiles based on modular grader definitions.
# Each grader is assigned a type according to the specified proportions.
# The function stores the distribution type and a list of parameters.
generate_graders <- function(n_graders, grader_profiles, grader_proportions) {
  types <- names(grader_profiles)
  assigned_types <- sample(types, size = n_graders, replace = TRUE, prob = grader_proportions)
  # Retrieve the distribution type for each assigned grader type.
  dist_types <- sapply(assigned_types, function(x) grader_profiles[[x]]$dist_type)
  # Retrieve the parameter list for each assigned grader type.
  dist_params <- lapply(assigned_types, function(x) grader_profiles[[x]]$params)
  
  graders <- data.frame(grader_id = 1:n_graders,
                        grader_type = assigned_types,
                        dist_type = dist_types,
                        stringsAsFactors = FALSE)
  # Create a list column for distribution parameters.
  graders$dist_params <- I(dist_params)
  return(graders)
}

# Function to assign papers to graders.
# Each paper is assigned to two distinct graders.
# A rejection sampling approach ensures that each grader’s load is between min_papers and max_papers.
assign_papers <- function(n_papers, n_graders, min_papers, max_papers) {
  valid_assignment <- FALSE
  attempt <- 1
  max_attempts <- 10000
  
  while(!valid_assignment && attempt <= max_attempts) {
    # Create an assignment data frame: 2 rows per paper (one for each grader)
    assignments <- data.frame(paper_id = rep(1:n_papers, each = 2),
                              grader_id = NA)
    
    # For each paper, randomly select two distinct graders
    for (i in 1:n_papers) {
      assignments$grader_id[(2*i - 1):(2*i)] <- sample(1:n_graders, size = 2, replace = FALSE)
    }
    
    # Calculate the number of assignments per grader
    grader_loads <- table(assignments$grader_id)
    
    if(all(grader_loads >= min_papers & grader_loads <= max_papers)) {
      valid_assignment <- TRUE
    }
    attempt <- attempt + 1
  }
  
  if(!valid_assignment) {
    stop("Failed to generate a valid assignment within the maximum number of attempts.")
  }
  
  return(assignments)
}

# Function to simulate the grader's score given the paper's true rating.
# This function uses the grader's distribution type and parameters to generate an error,
# which is then added to the true rating.
simulate_score <- function(true_rating, dist_type, params) {
  if(dist_type == "normal") {
    # For a normal distribution, sample an error with given bias and sd.
    error <- rnorm(1, mean = params$bias, sd = params$sd)
    score <- true_rating + error
  } else if(dist_type == "gamma") {
    # For a gamma distribution, sample an error and center it by subtracting its mean.
    error <- rgamma(1, shape = params$shape, scale = params$scale) - (params$shape * params$scale)
    score <- true_rating + error
  } else {
    stop("Distribution type not recognized.")
  }
  # Ensure the score is an integer between 1 and 7.
  score <- round(score)
  score <- min(max(score, 1), 7)
  return(score)
}

# Function to simulate grading.
# It merges paper, assignment, and grader information and uses simulate_score to generate the score.
simulate_grading <- function(papers, assignments, graders) {
  # Merge assignments with paper true ratings and grader profiles
  sim_data <- merge(assignments, papers, by = "paper_id")
  sim_data <- merge(sim_data, graders, by = "grader_id")
  
  # Apply the simulate_score function row-wise.
  sim_data <- sim_data %>%
    rowwise() %>%
    mutate(simulated_score = simulate_score(true_rating, dist_type, dist_params)) %>%
    ungroup()
  
  return(sim_data)
}

# Function to compute benchmarking metrics.
# Calculates mean absolute error (MAE), mean squared error (MSE), and classification accuracy.
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

# ---------------------------
# 3. Main Simulation Function
# ---------------------------

run_simulation <- function(n_papers, n_graders, min_papers, max_papers, rating_probs,
                           grader_profiles, grader_proportions) {
  papers <- generate_papers(n_papers, rating_probs)
  graders <- generate_graders(n_graders, grader_profiles, grader_proportions)
  assignments <- assign_papers(n_papers, n_graders, min_papers, max_papers)
  sim_data <- simulate_grading(papers, assignments, graders)
  benchmarks <- compute_benchmarks(sim_data)
  
  return(list(papers = papers,
              graders = graders,
              assignments = assignments,
              sim_data = sim_data,
              benchmarks = benchmarks))
}

# ---------------------------
# 4. Run the Simulation
# ---------------------------

simulation_results <- run_simulation(n_papers, n_graders, min_papers, max_papers, rating_probs,
                                     grader_profiles, grader_proportions)

# Display benchmark results
print(simulation_results$benchmarks)






# Assume simulation_results$sim_data is available from the simulation.
# Create a count table for each (true_rating, simulated_score) pair.
score_counts <- simulation_results$sim_data %>%
  group_by(true_rating, simulated_score) %>%
  summarise(count = n(), .groups = "drop")

# Ensure all 7 x 7 combinations are represented (fill missing combinations with 0)
all_combinations <- expand.grid(true_rating = 1:7, simulated_score = 1:7)
score_counts <- left_join(all_combinations, score_counts, by = c("true_rating", "simulated_score"))
score_counts$count[is.na(score_counts$count)] <- 0

# Create the heatmap with ggplot2
heatmap_plot <- ggplot(score_counts, aes(x = factor(true_rating), y = factor(simulated_score), fill = count)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = count), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "True Rating", y = "Simulated Score",
       title = "Heatmap of True Ratings vs Simulated Scores") +
  theme_minimal()

# Display the plot
print(heatmap_plot)







