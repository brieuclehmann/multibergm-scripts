# Load multibergm (install from github)
# devtools::install_github("brieuclehmann/multibergm")
library(multibergm)
library(patchwork)
library(ggplot2)
library(dplyr)
theme_set(theme_minimal(base_size = 12))
set.seed(1)
dir.create("plots", showWarnings = FALSE)

### SIMULATED DATA
# One group
n_nets_range <- c(10, 20, 50)
burn_in <- 2000

out_df <- data.frame(iteration = integer(), stat = character(), 
                     estimate = double(), n_nets = integer())
for (n_nets in n_nets_range) {
  
  fit_file <- paste0("output/single_n", n_nets, ".RDS")
  fit <- readRDS(fit_file)
  
  post_iters      <- seq(burn_in + 1L, fit$main_iters)
  samples <- subset(fit$params, iters = post_iters)
  output <- get("mu_pop", samples)
  output <- abind::adrop(unclass(output), 1)
  
  model_terms <- param_names(fit$control$model)
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- model_terms
  dimnames(output) <- dim_names
  truth_df <- data.frame(stat = model_terms, mean = c(-3, 0.5, 0.5))
  
  var_names <- c("iteration", "stat")
  this_df <- reshape2::melt(output, varnames = var_names, value.name = "estimate")
  this_df$n_nets <- n_nets
  out_df <- rbind(out_df, this_df)
  
  if (n_nets == 10) {
    p1 <- multibergm:::densityplot(output) +
      geom_vline(data = truth_df, aes(xintercept = mean), color = "red")
    p2 <- multibergm:::traceplot(output)
    p3 <- multibergm:::autocorrplot(output)
    
    p_out <- p1 + p2 + p3
    plot_file <- paste0("plots/single_n", n_nets, ".png")
    ggsave(plot_file, p_out, width = 7, height = 5)
    
    p_gof <- gof(fit, burn_in = burn_in, aux_iters = 5000)
    gof_file <- paste0("plots/single_n", n_nets, "gof.png")
    ggsave(gof_file, p_gof, width = 7, height = 5)
  }
  
}
out_df$n_nets <- as.factor(out_df$n_nets)
p_combi <- ggplot(out_df, aes(x = estimate, fill = n_nets, group = n_nets)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = truth_df, aes(xintercept = mean), color = "red") +
  facet_wrap("stat", ncol = 1, scales = "free") +
  ylab("density") +
  scale_fill_viridis_d()
combi_file <- paste0("plots/single_combi.png")
ggsave(combi_file, p_combi, width = 6, height = 4)


# TWO GROUP

burn_in = 2000
out_df <- data.frame(iteration = integer(), group = integer(), 
                     stat = character(), estimate = double(), 
                     n_nets = integer())
for (n_nets in n_nets_range) {
  
  fit_file <- paste0("output/twogrp_n", n_nets, ".RDS")
  fit <- readRDS(fit_file)
  
  post_iters      <- seq(burn_in + 1L, fit$main_iters)
  samples <- subset(fit$params, iters = post_iters)
  output <- get("mu_group", samples)
  output <- abind::adrop(unclass(output), 1)
  
  model_terms <- param_names(fit$control$model)
  n_dim <- length(dim(output))
  dim_names <- vector("list", n_dim)
  dim_names[[n_dim]] <- model_terms
  dimnames(output) <- dim_names
  truth_df <- data.frame(stat = rep(model_terms, 2), 
                         mean = c(-3, 0.5, 0.5, -2.6, 0.5, 0.2),
                         group = as.factor(rep(c(1,2), each = 3)))
  
  var_names <- c("iteration", "group", "stat")
  this_df <- reshape2::melt(output, varnames = var_names, value.name = "estimate")
  this_df$n_nets <- n_nets
  out_df <- rbind(out_df, this_df)
}
out_df$n_nets <- as.factor(out_df$n_nets)
out_df$group  <- as.factor(out_df$group)
edge_out_df <- filter(out_df, stat == "edges") %>% mutate(stat = droplevels(stat))
edge_truth_df <- filter(truth_df, stat == "edges") %>% mutate(stat = droplevels(stat))
p_edge <- ggplot(edge_out_df, aes(x = estimate, fill = group, group = group)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = edge_truth_df, aes(xintercept = mean, colour = group)) +
  facet_wrap(c("stat", "n_nets"), 
             labeller = labeller(stat = label_value, n_nets = label_both)) +
  ylab("") + xlab("")

hem_out_df <- filter(out_df, stat == "nodematch.hemisphere") %>% mutate(stat = droplevels(stat))
hem_truth_df <- filter(truth_df, stat == "nodematch.hemisphere") %>% mutate(stat = droplevels(stat))
p_hem <- ggplot(hem_out_df, aes(x = estimate, fill = group, group = group)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = hem_truth_df, aes(xintercept = mean, colour = group)) +
  facet_wrap(c("stat", "n_nets"), 
             labeller = labeller(stat = label_value, n_nets = label_both)) +
  ylab("density") + xlab("")

gwesp_out_df <- filter(out_df, stat == "gwesp.fixed.0.9") %>% mutate(stat = droplevels(stat))
gwesp_truth_df <- filter(truth_df, stat == "gwesp.fixed.0.9") %>% mutate(stat = droplevels(stat))
p_gwesp <- ggplot(gwesp_out_df, aes(x = estimate, fill = group, group = group)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = gwesp_truth_df, aes(xintercept = mean, colour = group)) +
  facet_wrap(c("stat", "n_nets"), 
             labeller = labeller(stat = label_value, n_nets = label_both)) +
  ylab("")

p_combi <- p_edge / p_hem / p_gwesp +
  plot_layout(guides = 'collect')
combi_file <- paste0("plots/twogrp_combi.png")
ggsave(combi_file, p_combi, width = 7, height = 5)

### BRAIN NETWORKS

fit_file <- "output/brain.RDS"
fit <- readRDS(fit_file)

burn_in <- 2000
post_iters      <- seq(burn_in + 1L, fit$main_iters)
samples <- subset(fit$params, iters = post_iters)
output <- get("mu_group", samples)
output <- abind::adrop(unclass(output), 1)

model_terms <- param_names(fit$control$model)
n_dim <- length(dim(output))
dim_names <- vector("list", n_dim)
dim_names[[n_dim]] <- model_terms
dimnames(output) <- dim_names

p1 <- multibergm:::densityplot(output) + 
  scale_fill_discrete(labels = c("young", "old"))
p2 <- multibergm:::traceplot(output) + 
  theme(legend.position = 'none')
p3 <- multibergm:::autocorrplot(output) + 
  theme(legend.position = 'none')

p_out <- p1 + p2 + p3 + plot_layout(guides = "collect")
plot_file <- "plots/brain.png"
ggsave(plot_file, p_out, width = 8, height = 5)


### Goodness of fit
sample_size <- 200
aux_iters <- 20000

post_iters     <- seq(burn_in + 1L, fit$main_iters, 1)
fit$params  <- subset(fit$params, iters = post_iters)

grp_label <- c("Young", "Old")
fit$groups <- factor(grp_label[fit$groups], levels = grp_label)
fit$networks <- mapply(function(x,y) set.network.attribute(x, "group", y),
                       fit$networks, fit$groups, SIMPLIFY = FALSE)
obs_df <- get_net_stats(fit$networks, fit$formula, "gof")

# Get posterior samples
param <- "mu_group"
output <- get(param, fit$params)

n_iters <- mcmcr::niters(fit$params$theta)
coefs   <- subset(output, iters = sample(n_iters, sample_size))
coefs   <- abind::adrop(unclass(coefs), 1)

sim_df <- obs_df[0, ]
for (j in seq_along(unique(fit$groups))) {
  g <- unique(fit$groups)[j]
  group_nets <- fit$networks[fit$groups == g]
  
  for (i in seq_len(sample_size)) {
    y         <- group_nets[[sample(length(group_nets), 1)]]
    myformula <- statnet.common:::nonsimp_update.formula(fit$formula, y ~.,
                                                         from.new = "y")
    
    if (length(dim(coefs)) == 2) {
      this_coef <- coefs[i,]
    } else {
      this_coef <- coefs[i,j,]
    }
    
    net_sim <- ergm::simulate_formula(myformula,
                                      coef = this_coef,
                                      nsim = 1,
                                      constraints = fit$constraints,
                                      control = control.simulate.formula(
                                        MCMC.burnin = aux_iters
                                      )
    )
    
    this_df <- get_net_stats(net_sim, myformula, "gof")
    this_df$group <- g
    sim_df  <- rbind(sim_df, this_df)
  }
}

sim_df <- sim_df %>%
  group_by(Stat, n, group) %>%
  summarise(Group = quantile(Value, 0.5),
            Lower = quantile(Value, 0.05),
            Upper = quantile(Value, 0.95))

deg_obs <- filter(obs_df, Stat == "degree")
deg_sim <- filter(sim_df, Stat == "degree")
p1 <- multibergm:::plot_dist_gof(deg_obs, deg_sim, "Degree") +
  guides(fill = "none", color = "none") + ylab("")

geodist_obs <- filter(obs_df, Stat == "distance" & !is.na(n))
geodist_sim <- filter(sim_df, Stat == "distance" & !is.na(n))
n_max_sim <- max(1L, geodist_sim$n[geodist_sim$Upper > 0], na.rm = TRUE)
n_max_obs <- max(1L, geodist_obs$n[geodist_obs$Value > 0], na.rm = TRUE)
n_max <- max(n_max_sim, n_max_obs)

geodist_sim <- geodist_sim %>%
  filter(n <= n_max)
geodist_obs <- geodist_obs %>%
  filter(n <= n_max)

inf_obs_df <- obs_df %>%
  filter(Stat == "distance" & is.na(n)) %>%
  mutate(n = "Inf")
inf_sim_df <- sim_df %>%
  filter(Stat == "distance" & is.na(n)) %>%
  mutate(n = "Inf")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
blue <- gg_color_hue(2)[2]

p2a <- ggplot(geodist_sim) +
  geom_blank() +
  geom_boxplot(data = filter(geodist_obs, group == "Young"), aes(x = n, y = Value, group = n)) +
  geom_line(data = filter(geodist_sim, group == "Young"), aes(y = Group, x = n, color = group),
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_ribbon(data = filter(geodist_sim, group == "Young"),
              aes(ymin = Lower, ymax = Upper, x = n,
                  fill = group),
              alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
  scale_linetype_discrete(name = NULL, labels = "Observed data") +
  ylab("Probability") + xlab("") + ylim(c(0, 0.6)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "null"))

p2b <- ggplot(filter(inf_obs_df, group == "Young"), aes(x = n, y = Value, group = n)) +
  geom_boxplot() +
  geom_point(data = filter(inf_sim_df, group == "Young"), aes(y = Group, x = n, color = group),
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_linerange(data = filter(inf_sim_df, group == "Young"),
              aes(ymin = Lower, ymax = Upper, x = n,
                  colour = group),
              alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE, size = 10) +
  scale_linetype_discrete(name = NULL, labels = "Observed data") +
  ylab("") + xlab("     Geodesic distance") + ylim(c(0,1)) + 
  theme(plot.margin = unit(c(0, 10, 0, 0), "null"))
  

p2c <- ggplot(filter(geodist_obs, group == "Old"), aes(x = n, y = Value, group = n)) +
  geom_boxplot() +
  geom_line(data = filter(geodist_sim, group == "Old"), aes(y = Group, x = n, color = group),
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_ribbon(data = filter(geodist_sim, group == "Old"),
              aes(ymin = Lower, ymax = Upper, x = n,
                  fill = group),
              alpha = 0.4, inherit.aes = FALSE, show.legend = FALSE) +
  scale_color_discrete(h = c(195, 15), name = NULL) +
  scale_fill_discrete(h = c(195, 15), name = NULL) +
  scale_linetype_discrete(name = NULL, labels = "Observed data") +
  ylab("") + xlab("") + ylim(c(0, 0.6)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "null"))

p2d <- ggplot(filter(inf_obs_df, group == "Old"), aes(x = n, y = Value, group = n)) +
  geom_boxplot() +
  geom_point(data = filter(inf_sim_df, group == "Old"), aes(y = Group, x = n, color = group),
             inherit.aes = FALSE, show.legend = FALSE) +
  geom_linerange(data = filter(inf_sim_df, group == "Old"),
                 aes(ymin = Lower, ymax = Upper, x = n,
                     colour = group),
                 alpha = 0.4, inherit.aes = FALSE, size = 10, show.legend = FALSE) +
  scale_linetype_discrete(name = NULL, labels = "Observed data") +
  scale_color_discrete(h = c(195, 15), name = NULL) +
  ylab("") + xlab("") + ylim(c(0,1)) + 
  theme(plot.margin = unit(c(0, 10, 0, 0), "null"))

p2 <- p2a + p2b + p2c + p2d + 
  plot_layout(widths = c(20, 1, 20, 1), guides = "collect")

esp_obs <- filter(obs_df, Stat == "esp")
esp_sim <- filter(sim_df, Stat == "esp")
p3 <- multibergm:::plot_dist_gof(deg_obs, deg_sim, "Edgewise shared partners") +
  ylab("") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

p_gof <- (p1 / p2 / p3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
gof_file <- paste0("plots/brain_gof.png")
ggsave(gof_file, p_gof, width = 6, height = 7)
