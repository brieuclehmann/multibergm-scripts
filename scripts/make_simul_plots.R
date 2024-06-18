library(multibergm)
source("scripts/plot_utils.R")
library(patchwork)
theme_set(theme_minimal(base_size = 12))
set.seed(1)

### SIMULATED DATA
# One group
n_nets_range <- c(10, 20, 50)
burn_in <- 2000
decay <- 0.9

out_df <- data.frame(iteration = integer(), var= integer(), stat = character(), 
                     estimate = double(), n_nets = integer())
for (n_nets in n_nets_range) {
  
  fit_file <- paste0("output/single_n", n_nets, ".RDS")
  fit <- readRDS(fit_file)
  
  post_iters      <- seq(burn_in + 1L, fit$main_iters)
  samples <- subset(fit$params, iters = post_iters)
  output <- get("mu", samples)
  output <- abind::adrop(unclass(output), 1)
  #if (dim(output)[2] == 1)
  #  output <- abind::adrop(output, 2)
  
  ergm_terms <- param_names(fit$control$model)
  ergm_terms <- c("edges", "gwesp", "hemisphere")
  
  model_terms <- c("Intercept")
  
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- ergm_terms
  dim_names[[2]] <- model_terms
  dimnames(output) <- dim_names
  truth_df <- data.frame(stat = ergm_terms, var = model_terms, 
                         mean = c(-3, 0.5, 0.5))
  truth_df$facet <- paste(truth_df$stat, truth_df$var, sep = ": ")
  
  var_names <- c("iteration", "var", "stat")
  this_df <- melt(output, varnames = var_names, value.name = "estimate")
  this_df$n_nets <- n_nets
  out_df <- rbind(out_df, this_df)
  
  p1 <- densityplot2(output) +
    geom_vline(data = truth_df, aes(xintercept = mean), color = "red")
  p2 <- traceplot2(output)
  p3 <- autocorrplot2(output)
  
  p_out <- cowplot::plot_grid(p1, p2, p3, nrow = 1)
  #p_out <- p1 + p2 + p3
  plot_file <- paste0("plots/single_n", n_nets, ".png")
  ggsave(plot_file, p_out, width = 7, height = 5)
  
  
  ### GOODNESS-OF-FIT ###
  sample_size <- 100
  n_iters <- dim(output)[1]
  coefs   <- output[sample(n_iters, sample_size),1,]
  p_gof <- gof(fit, coefs = coefs)
  gof_file <- paste0("plots/single_n", n_nets, "gof.png")
  ggsave(gof_file, p_gof, width = 7, height = 5)
  
}
out_df$n_nets <- as.factor(out_df$n_nets)
out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")

p_combi <- ggplot(out_df, aes(x = estimate, fill = n_nets, group = n_nets)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = truth_df, aes(xintercept = mean), color = "red") +
  facet_wrap("facet", ncol = 1, scales = "free") +
  ylab("density") +
  scale_fill_viridis_d()
combi_file <- paste0("plots/single_combi.png")
ggsave(combi_file, p_combi, width = 6, height = 4)


### Continuous plot ###
out_df <- data.frame(iteration = integer(), var= integer(), stat = character(), 
                     estimate = double(), n_nets = integer())
for (n_nets in n_nets_range) {
  
  fit_file <- paste0("output/cts_n", n_nets, "_0.9.RDS")
  fit <- readRDS(fit_file)
  
  post_iters      <- seq(burn_in + 1L, fit$main_iters)
  samples <- subset(fit$params, iters = post_iters)
  output <- get("mu", samples)
  output <- abind::adrop(unclass(output), 1)

  
  ergm_terms <- param_names(fit$control$model)
  ergm_terms <- c("edges", "gwesp", "hemisphere")
  
  model_terms <- c("Intercept", "x")
  
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- ergm_terms
  dim_names[[2]] <- model_terms
  dimnames(output) <- dim_names
  truth_df <- data.frame(stat = rep(ergm_terms, length(model_terms)),
                         var = rep(model_terms, each = length(ergm_terms)), 
                         mean = c(-3, 0.5, 0.5,
                                  0.4, 0, -0.3))
  truth_df$facet <- paste(truth_df$stat, truth_df$var, sep = ": ")
  
  var_names <- c("iteration", "var", "stat")
  this_df <- melt(output, varnames = var_names, value.name = "estimate")
  this_df$n_nets <- n_nets
  out_df <- rbind(out_df, this_df)
  
  p1 <- densityplot2(output) +
    geom_vline(data = truth_df, aes(xintercept = mean), color = "red")
  p2 <- traceplot2(output)
  p3 <- autocorrplot2(output)
  
}
out_df$n_nets <- as.factor(out_df$n_nets)
out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")

p_combi <- ggplot(out_df, aes(x = estimate, fill = n_nets, group = n_nets)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = truth_df, aes(xintercept = mean), color = "red") +
  facet_wrap("facet", ncol = 2, scales = "free") +
  ylab("density") +
  scale_fill_viridis_d()
combi_file <- paste0("plots/cts_combi.png")
ggsave(combi_file, p_combi, width = 6, height = 4)


# TWO GROUP

burn_in = 2000
out_df <- data.frame(iteration = integer(), group = integer(), 
                     stat = character(), estimate = double(), 
                     n_nets = integer())
decay <- 0.9
for (n_nets in n_nets_range) {
  
  fit_file <- paste0("output/twogrp_n", n_nets, "_", decay, ".RDS")
  fit <- readRDS(fit_file)
  
  post_iters      <- seq(burn_in + 1L, fit$main_iters)
  samples <- subset(fit$params, iters = post_iters)
  output <- get("mu", samples)
  output <- abind::adrop(unclass(output), 1)
  
  ergm_terms <- param_names(fit$control$model)
  ergm_terms <- c("edges", "gwesp", "hemisphere")
  
  model_terms <- c("Group 1", "Group 2")
  
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- ergm_terms
  dim_names[[2]] <- model_terms
  dimnames(output) <- dim_names
  truth_df <- data.frame(stat = rep(ergm_terms, 2), 
                         mean = c(-3, 0.5, 0.5, -2.6, 0.5, 0.2),
                         var = as.factor(rep(c("Group 1", "Group 2"), each = 3)))
  truth_df$facet <- paste(truth_df$stat, truth_df$var, sep = ": ")
  
  var_names <- c("iteration", "var", "stat")
  this_df <- melt(output, varnames = var_names, value.name = "estimate")
  this_df$n_nets <- n_nets
  out_df <- rbind(out_df, this_df)
  
}

out_df$n_nets <- as.factor(out_df$n_nets)
out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")

p_combi <- ggplot(out_df, aes(x = estimate, fill = n_nets, group = n_nets)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = truth_df, aes(xintercept = mean), color = "red") +
  facet_wrap("facet", ncol = 2, scales = "free") +
  ylab("density") +
  scale_fill_viridis_d()
combi_file <- paste0("plots/twogrp_combi.png")
ggsave(combi_file, p_combi, width = 6, height = 4)


