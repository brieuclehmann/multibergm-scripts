
# Load multibergm (install from github)
# devtools::install_github("brieuclehmann/multibergm")
library(multibergm)
source("scripts/generate_networks.R")
dir.create("output", showWarnings=FALSE)
set.seed(1)

# Set parameters
decay <- 0.9
n_nets  <- 10 # Vary between 10, 20, 50
n_batches <- 2 # number of available CPU cores (for parallelisation)
burn_in <- 2000
main_iters <- 10000
aux_iters <- 1000
tot_iters <- burn_in + main_iters
n_nodes <- 30
ergm_formula <- ~ edges + nodematch("hemisphere") + gwesp(decay, fixed = TRUE)
n_terms <- length(attr(terms(ergm_formula), "term.labels"))

# Generate ERGM coefficients
coef_mean <- c(-3, 0.5, 0.5)
coef_cov <- matrix(0, n_terms, n_terms)
diag(coef_cov)  <- c(1, 0.5, 0.5)
coef_cov[1,2] <- coef_cov[2,1] <- -0.5
coef_cov_scaled <- coef_cov / 50

coef_multi <- mvtnorm::rmvnorm(n_nets, coef_mean, coef_cov_scaled)

# Generate networks
nets_multi <- generate_networks(ergm_formula, n_nodes, coef_multi)

# Fit multibergm
multi_formula <- update(ergm_formula, nets_multi ~ .)
control <- control_multibergm(multi_formula, 
                              proposal_update_max = burn_in, 
                              aux_iters = aux_iters,
                              n_batches = n_batches)
fit_multi <- multibergm(multi_formula, main_iters = tot_iters, control = control)

# Save output
out_file <- file.path("output", paste0("single_n", n_nets, ".RDS"))
saveRDS(fit_multi, out_file)
