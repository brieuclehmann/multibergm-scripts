
library(multibergm)
source("scripts/generate_networks.R")
set.seed(1)

# Set parameters
decay <- as.double(Sys.getenv("decay"))
n_nets  <- as.integer(Sys.getenv("n_nets"))
n_batches <- as.integer(Sys.getenv("n_batches"))
burn_in <- 2000
main_iters <- 20000
aux_iters <- 1000
tot_iters <- burn_in + main_iters
n_nodes <- 30
ergm_formula <- ~ edges + nodematch("hemisphere") + gwesp(decay, fixed = TRUE)
n_terms <- length(attr(terms(ergm_formula), "term.labels"))

n_nodes  <- as.integer(Sys.getenv("n_nodes"))
if (is.na(n_nodes)) n_nodes <- 30

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
#nets_common <- generate_networks(ergm_formula, n_nodes, coef_common)
#get_net_stats.list(nets_common, ergm_formula, "model")

nets_multi <- generate_networks(ergm_formula, n_nodes, coef_multi)
group_ind <- rep(c(0, 1), each = n_nets)
multi_formula <- update(ergm_formula, nets_multi ~ .)
nets_multi <- lapply(seq(length(nets_multi)), 
                     function(x) set.network.attribute(nets_multi[[x]], 
                                                       "group", 
                                                       group_ind[x]))

get_net_stats.list(nets_multi[1:n_nets], ergm_formula, "model") %>% apply(2, quantile, c(0, 0.5, 1))
get_net_stats.list(nets_multi[-(1:n_nets)], ergm_formula, "model") %>% apply(2, quantile, c(0, 0.5, 1))

# Fit multibergm
model_formula <- ~ 1 + group
model_matrix <- get_model_matrix(multi_formula, model_formula)
model_matrix[ ,1] <- model_matrix[ ,1] - model_matrix[ ,2]

control <- control_multibergm(multi_formula,
                              mod_mat = model_matrix,
                              proposal_update_max = burn_in, 
                              aux_iters = aux_iters,
                              n_batches = n_batches)

fit_multi <- multibergm(multi_formula,
                        model_matrix = model_matrix,
                        groups = group_ind,
                        main_iters = tot_iters, 
                        control = control)

# Save output
out_file <- file.path("output", paste0("twogrp_n", n_nets, ".RDS"))
if (n_nodes != 30) {
    out_file <- file.path("output", paste0("two_grp_n", n_nets, "_k", n_nodes, ".RDS"))
}

saveRDS(fit_multi, out_file)
