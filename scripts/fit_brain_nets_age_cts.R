
library(multibergm)
source("scripts/generate_networks.R")
set.seed(1)

# Set parameters
n_nets  <- as.integer(Sys.getenv("n_nets"))
n_batches <- as.integer(Sys.getenv("n_batches"))
decay <- as.double(Sys.getenv("decay"))
covar <- Sys.getenv("covar")
model <- Sys.getenv("model")
print(model)

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
ages <- scan("data/ages.txt")
cattell <- scan("data/cattell.txt")
ages <- scale(ages)
#ages <- ages - min(ages, na.rm = TRUE)
cattell[is.nan(cattell)] <- NA
cattell <- scale(cattell)
#cattell <- cattell - min(cattell, na.rm = TRUE)

nets <- lapply(seq(length(nets)), 
               function(i) {
                 set.network.attribute(nets[[i]], "age", ages[i])
                 set.network.attribute(nets[[i]], "cattell", cattell[i])
                 }
               )

if (covar == "age") {
  net_ind <- round(seq(1, 587, length.out = n_nets))
} else {
  na_cattell <- which(is.na(cattell))
  net_ind <- order(cattell)[round(seq(1, 587 - length(na_cattell), length.out = n_nets))]
}

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

multi_formula <- update(ergm_formula, nets ~ .)

# if (covar == "age") {
#   model_formula <- ~ 1 + age
# } else {
#   model_formula <- ~ 1 + cattell
# }
model_formula <- ~ 1 + age*cattell
if (model == "quadratic") {
  model_formula <- ~ 1 + poly(age, degree = 2, raw = TRUE) + poly(cattell, degree = 2, raw = TRUE) + age:cattell
}
print(model_formula)

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
out_file <- file.path("output", paste0("brain_cts_n", n_nets, ".RDS"))
if (model == "quadratic") out_file <- file.path("output", paste0("brain_cts_n", n_nets, "_quadratic.RDS"))

saveRDS(fit_multi, out_file)
