
library(multibergm)
source("scripts/generate_networks.R")
set.seed(1)

# Set parameters
params <- commandArgs(trailingOnly = TRUE)
decay <- 0.9
n_nets  <- as.integer(Sys.getenv("n_nets"))
n_batches <- as.integer(Sys.getenv("n_batches"))
burn_in <- 2000
main_iters <- 10000

aux_iters <- as.integer(Sys.getenv("aux"))
if (is.na(aux_iters)) aux_iters <- 1000
tot_iters <- burn_in + main_iters

n_nodes <- as.integer(Sys.getenv("n_nodes"))
if (is.na(n_nodes)) n_nodes <- 30


ergm_formula <- nets_multi ~ edges + nodematch("hemisphere") + gwesp(decay, fixed = TRUE)
n_terms <- length(attr(terms(ergm_formula), "term.labels"))

# Generate ERGM coefficients
coef_mean <- c(-3, 0.5, 0.5)
coef_cov <- matrix(0, n_terms, n_terms)
diag(coef_cov)  <- c(1, 0.5, 0.5)
coef_cov[1,2] <- coef_cov[2,1] <- -0.5
coef_cov_scaled <- coef_cov / 50

coef_multi <- rmvnorm(n_nets, coef_mean, coef_cov_scaled)
# coef_common <- matrix(rep(coef_mean, each = n_nets), ncol = 3)


# Generate networks
nets_multi <- generate_networks(ergm_formula, n_nodes, coef_multi)
get_net_stats.list(nets_multi, ergm_formula, "model") %>% apply(2, range)

# nets_common <- generate_networks(ergm_formula, n_nodes, coef_common)
# get_net_stats.list(nets_common, ergm_formula, "model") %>% apply(2, range)

# Fit multibergm
multi_formula <- update(ergm_formula, nets_multi ~ .)
model_formula <- ~ 1
model_matrix <- get_model_matrix(ergm_formula, model_formula)

control <- control_multibergm(multi_formula,
                              mod_mat = model_matrix,
                              proposal_update_max = burn_in, 
                              aux_iters = aux_iters,
                              n_batches = n_batches)
fit_multi <- multibergm(multi_formula, main_iters = tot_iters, control = control)

# Save output
out_file <- file.path("output", paste0("single_n", n_nets, ".RDS"))
if (Sys.getenv("n_nodes") != "") {
    out_file <- file.path("output", paste0("single_n", n_nets, "_k", n_nodes, ".RDS"))
    if (aux_iters != 1000) {
        out_file <- file.path("output", paste0("single_n", n_nets, "_k", n_nodes, "_aux5k.RDS"))
    }
}

saveRDS(fit_multi, out_file)
