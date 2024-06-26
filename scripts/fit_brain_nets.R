
library(multibergm)
source("scripts/generate_networks.R")
set.seed(1)

# Set parameters
n_nets  <- as.integer(Sys.getenv("n_nets"))
n_batches <- as.integer(Sys.getenv("n_batches"))
decay <- as.double(Sys.getenv("decay"))

burn_in <- 2000
main_iters <- 20000
aux_iters <- 10000
tot_iters <- burn_in + main_iters
ergm_formula <- ~ edges + 
  nodematch("hemisphere") + 
  nodematch("homotopy") + 
  gwesp(decay, fixed = TRUE)

# Load networks
nets <- readRDS("data/constK3_n587.RDS")
net_ind <- c(1:n_nets, rev(length(nets):(length(nets) - n_nets + 1)))
nets <- nets[net_ind]
n_nodes <- network.size(nets[[1]])
node_info <- function(x) {
  n_nodes <- network.size(x)
  x %v% "hemisphere" <- rep(c("left", "right"), n_nodes / 2)
  x %v% "homotopy" <- rep(1:(n_nodes / 2), each = 2)
  x
}
nets <- lapply(nets, node_info)
# nets_common <- generate_networks(ergm_formula, n_nodes, coef_common)
# get_net_stats.list(nets_common, ergm_formula, "model") %>% apply(2, range)

# Fit multibergm
group_ind <- rep(c(0, 1), each = n_nets)
multi_formula <- update(ergm_formula, nets ~ .)
nets <- lapply(seq(length(nets)), 
               function(x) set.network.attribute(nets[[x]], 
                                                 "group", 
                                                 group_ind[x]))
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
out_file <- file.path("output", paste0("brain_n", n_nets, ".RDS"))
saveRDS(fit_multi, out_file)
