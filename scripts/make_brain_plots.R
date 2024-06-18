library(multibergm)
library(patchwork)
theme_set(theme_minimal(base_size = 9))

source("scripts/plot_utils.R")
set.seed(1)

### BRAIN NETWORKS

n_nets <- 100
fit_file <- paste0("output/brain_n", n_nets, ".RDS")
fit <- readRDS(fit_file)

burn_in <- 2000
main_iters <- 20000
aux_iters <- 10000
decay <- 0.9

post_iters      <- seq(burn_in + 1L, fit$main_iters)
samples <- subset(fit$params, iters = post_iters)
output <- get("mu", samples)
output <- abind::adrop(unclass(output), 1)

ergm_terms <- param_names(fit$control$model)
ergm_terms <- c("edges","hemisphere", "homotopy", "gwesp")

model_terms <- c("Young", "Old")

n_dim <- length(dim(output))
dim_names <- vector("list", n_dim)
dim_names[[n_dim]] <- ergm_terms
dim_names[[2]] <- model_terms
dimnames(output) <- dim_names

p1 <- densityplot_group(output) + theme(legend.position = 'none')
p2 <- traceplot_group(output) + theme(legend.position = 'none')
p3 <- autocorrplot_group(output)

p_out <- cowplot::plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, 1, 1.5))
#p_out <- p1 + p2 + p3 + plot_layout(guides = "collect")
plot_file <- "plots/brain_n100.png"
ggsave(plot_file, p_out, width = 8, height = 5)


### Goodness of fit
sample_size <- 200
aux_iters <- 20000

post_iters     <- seq(burn_in + 1L, fit$main_iters, 1)
samples <- subset(fit$params, iters = post_iters)
output <- get("mu", samples)

# Get posterior samples
coefs   <- subset(output, iters = sample(main_iters, sample_size))
coefs   <- abind::adrop(unclass(coefs), 1)

grp_label <- c("Young", "Old")
obs_young_df <- get_net_stats(fit$networks[1:100], fit$formula, "gof")
obs_young_df$group <- "Young"
obs_old_df <- get_net_stats(fit$networks[101:200], fit$formula, "gof")
obs_old_df$group <- "Old"

obs_df <- bind_rows(obs_young_df, obs_old_df)

sim_df <- obs_df[0, ]
for (j in c(1, 2)) {
  group_nets <- fit$networks[1:100]
  if (j == 2) group_nets <- fit$networks[101:200]
  
  for (i in seq_len(sample_size)) {
    y         <- group_nets[[sample(length(group_nets), 1)]]
    myformula <- nonsimp_update.formula(fit$formula, y ~.,
                                        from.new = "y")
    
    this_coef <- coefs[i,j,]
    net_sim <- ergm::simulate_formula(myformula,
                                      coef = this_coef,
                                      nsim = 1,
                                      constraints = fit$constraints,
                                      control = control.simulate.formula(
                                        MCMC.burnin = aux_iters
                                      )
    )
    
    this_df <- get_net_stats(net_sim, myformula, "gof")
    this_df$group <- grp_label[j]
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
p1 <- plot_dist_gof(deg_obs, deg_sim, "Degree") +
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
p3 <- plot_dist_gof(deg_obs, deg_sim, "Edgewise shared partners") +
  ylab("") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

p_gof <- (p1 / p2 / p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
gof_file <- paste0("plots/brain_n100gof.png")
ggsave(gof_file, p_gof, width = 6, height = 7)


### AGE ###

n_nets <- 100
fit_file <- paste0("output/brain_age_n", n_nets, ".RDS")
fit <- readRDS(fit_file)

samples <- subset(fit$params, iters = post_iters)
output <- get("mu", samples)
output <- abind::adrop(unclass(output), 1)

model_terms <- c("Intercept", "age")

n_dim <- length(dim(output))
dim_names <- vector("list", n_dim)
dim_names[[n_dim]] <- ergm_terms
dim_names[[2]] <- model_terms
dimnames(output) <- dim_names

zero_df <- data.frame(stat = rep(ergm_terms, length(model_terms)),
                      var = rep(model_terms, each = length(ergm_terms)), 
                      mean = c(NA, NA, NA, NA,
                               0, 0, 0, 0))

zero_df$facet <- paste(zero_df$stat, zero_df$var, sep = ": ")
zero_df$facet <- factor(zero_df$facet, unique(zero_df$facet), ordered = TRUE)

p1 <- densityplot2(output) +
  geom_vline(data = zero_df, aes(xintercept = mean), color = "red")

out_file <- paste0("plots/age_density.png")
ggsave(out_file, p1, width = 3, height = 4)


### AGE ###

n_nets <- 100
fit_file <- paste0("output/brain_cattell_n", n_nets, ".RDS")
fit <- readRDS(fit_file)

samples <- subset(fit$params, iters = post_iters)
output <- get("mu", samples)
output <- abind::adrop(unclass(output), 1)

model_terms <- c("Intercept", "Cattell")

n_dim <- length(dim(output))
dim_names <- vector("list", n_dim)
dim_names[[n_dim]] <- ergm_terms
dim_names[[2]] <- model_terms
dimnames(output) <- dim_names

zero_df <- data.frame(stat = rep(ergm_terms, length(model_terms)),
                       var = rep(model_terms, each = length(ergm_terms)), 
                       mean = c(NA, NA, NA, NA,
                                0, 0, 0, 0))

zero_df$facet <- paste(zero_df$stat, zero_df$var, sep = ": ")
zero_df$facet <- factor(zero_df$facet, unique(zero_df$facet), ordered = TRUE)

p1 <- densityplot2(output) +
  geom_vline(data = zero_df, aes(xintercept = mean), color = "red")

out_file <- paste0("plots/cattell_density.png")
ggsave(out_file, p1, width = 3, height = 4)
