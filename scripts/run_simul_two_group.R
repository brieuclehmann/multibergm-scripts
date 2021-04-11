# Load multibergm (install from github)
# devtools::install_github("brieuclehmann/multibergm")
library(multibergm)
source("scripts/generate_networks.R")
dir.create("output", showWarnings=FALSE)
set.seed(1)

# Set parameters
decay <- 0.9
n_nets  <- 10 # number of nets per group, vary between 10, 20, 50
n_batches <- 2 # number of available CPU cores (for parallelisation)
burn_in <- 2000
main_iters <- 20000
aux_iters <- 1000
tot_iters <- burn_in + main_iters
n_nodes <- 30
ergm_formula <- ~ edges + nodematch("hemisphere") + gwesp(decay, fixed = TRUE)
n_terms <- length(attr(terms(ergm_formula), "term.labels"))

# Generate ERGM coefficients
coef_mean_a <- c(-3,   0.5, 0.5)
coef_mean_b <- c(-2.6, 0.5, 0.2)
coef_cov <- matrix(0, n_terms, n_terms)
diag(coef_cov)  <- c(1, 0.5, 0.5)
coef_cov[1,2] <- coef_cov[2,1] <- -0.5
coef_cov_scaled <- coef_cov / 50

coef_multi  <- rbind(rmvnorm(n_nets, coef_mean_a, coef_cov_scaled),
                     rmvnorm(n_nets, coef_mean_b, coef_cov_scaled))

# Generate networks
nets_multi <- generate_networks(ergm_formula, n_nodes, coef_multi)

# Fit multibergm
group_ind <- rep(c(1,2), each = n_nets)
multi_formula <- update(ergm_formula, nets_multi ~ .)
control <- control_multibergm(multi_formula,
                              groups = group_ind,
                              proposal_update_max = burn_in, 
                              aux_iters = aux_iters,
                              n_batches = n_batches)

fit_multi <- multibergm(multi_formula, 
                        groups = group_ind,
                        main_iters = tot_iters, 
                        control = control)

# Save output
out_file <- file.path("output", paste0("twogrp_n", n_nets, ".RDS"))
saveRDS(fit_multi, out_file)
