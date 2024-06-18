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

# Generate ERGM coefficients
coef_mean_a <- c(-3,   0.5, 0.5)
coef_mean_b <- c(-2.6, 0.5, 0.2)
coef_cov <- matrix(0, n_terms, n_terms)
diag(coef_cov)  <- c(1, 0.5, 0.5)
coef_cov[1,2] <- coef_cov[2,1] <- -0.5
coef_cov_scaled <- coef_cov / 50


coef_multi  <- - sweep(matrix(coef_mean_b - coef_mean_a) %*% t(matrix((n_nets - seq(n_nets)) / (n_nets - 1))), 
                       1, coef_mean_b)

# Add some noise
coef_multi <- t(coef_multi) + rmvnorm(n_nets, rep(0, n_terms), coef_cov_scaled)

# Generate networks
#nets_common <- generate_networks(ergm_formula, n_nodes, coef_common)
#get_net_stats.list(nets_common, ergm_formula, "model")

nets_multi <- generate_networks(ergm_formula, n_nodes, coef_multi)
multi_formula <- update(ergm_formula, nets_multi ~ .)
x <- (seq(n_nets) - 1) / (n_nets - 1)
nets_multi <- lapply(seq(length(nets_multi)), 
                     function(i) set.network.attribute(nets_multi[[i]], 
                                                       "x", 
                                                       x[i]))

#get_net_stats.list(nets_multi[1:n_nets], ergm_formula, "model") %>% apply(2, quantile, c(0, 0.5, 1))
#get_net_stats.list(nets_multi[-(1:n_nets)], ergm_formula, "model") %>% apply(2, quantile, c(0, 0.5, 1))

# Fit multibergm
model_formula <- ~ 1 + x
model_matrix <- get_model_matrix(multi_formula, model_formula)

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
out_file <- file.path("output", paste0("cts_n", n_nets, ".RDS"))
saveRDS(fit_multi, out_file)
