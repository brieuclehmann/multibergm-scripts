# Load multibergm (install from github)
# devtools::install_github("brieuclehmann/multibergm")
library(multibergm)
set.seed(1)
dir.create("output", showWarnings=FALSE)

# Set parameters
n_nets <- 100
n_batches <- 2 # number of available CPU cores (for parallelisation)
decay <- 0.9
burn_in <- 2000
main_iters <- 20000
aux_iters <- 10000
tot_iters <- burn_in + main_iters
ergm_formula <- ~ edges + 
  nodematch("hemisphere") + 
  nodematch("homotopy") + 
  gwesp(decay, fixed = TRUE)

# Load networks
nets <- readRDS("data/brain_nets.RDS")
net_ind <- c(1:n_nets, rev(length(nets):(length(nets) - n_nets + 1)))
nets <- nets[net_ind]
n_nodes <- network.size(nets[[1]])
node_info <- function(x) {
  n_nodes <- network.size(x)
  x %v% "hemisphere" <- rep(c("left", "right"), n_nodes / 2)
  x %v% "homotopy" <- rep(1:(n_nodes / 2), each = 2)
  x
}
nets <- lapply(nets, node_info) # Add covariate information

# Fit multibergm
group_ind <- rep(c(1, 2), each = n_nets)
multi_formula <- update(ergm_formula, nets ~ .)
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
out_file <- file.path("output", "brain.RDS")
saveRDS(fit_multi, out_file)
