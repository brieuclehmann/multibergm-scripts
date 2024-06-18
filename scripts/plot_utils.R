
traceplot2 <- function(output) {
  
  var_names <- c("iteration", "var", "stat")
  out_df <- melt(output, varnames = var_names, value.name = "estimate")
  out_df$iteration <- as.integer(out_df$iteration)
  out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")
  ncol = length(unique(out_df$var))
  p <- ggplot(out_df, aes(x = .data$iteration, y = .data$estimate)) +
    geom_line(colour = "black") +
    facet_wrap("facet", ncol = ncol, scales = "free_y")
  
  p
}

traceplot_group <- function(output) {
  
  var_names <- c("iteration", "var", "stat")
  out_df <- melt(output, varnames = var_names, value.name = "estimate")
  out_df$iteration <- as.integer(out_df$iteration)
  ncol = length(unique(out_df$var))
  p <- ggplot(out_df, aes(x = .data$iteration, y = .data$estimate)) +
    geom_line(aes(colour = var)) +
    facet_wrap("stat", ncol = 1, scales = "free_y")
  
  p
}



densityplot2 <- function(output) {
  
  var_names <- c("iteration", "var", "stat")
  out_df <- melt(output, varnames = var_names, value.name = "estimate")
  out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")

  out_df$facet <- factor(out_df$facet, unique(out_df$facet), ordered = TRUE)
  
  ncol = length(unique(out_df$var))
  p <- ggplot(out_df, aes(x = .data$estimate)) +
    geom_density(fill = "black", alpha = 0.5, colour = NA) +
    facet_wrap(c("facet"), ncol = ncol, scales = "free") +
    ylab("density")
  
  p 
}

densityplot_group <- function(output, param) {
  
  var_names <- c("iteration", "var", "stat")
  out_df <- melt(output, varnames = var_names, value.name = "estimate")
  out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")
  p <- ggplot(out_df, aes(x = .data$estimate)) +
    geom_density(aes(fill = .data$var), alpha = 0.5, colour = NA) +
    facet_wrap(c("stat"), ncol = 1, scales = "free") +
    ylab("density")
  
  p 
}


autocorrplot2 <- function(output, lagmax = 40) {
  
  # Get autocorrelation for each model term
  autocorr <- apply(output, c(2,3),
                    function(x) acf(x, lag.max = lagmax, plot = F)$acf)
  
  var_names <- c("lag", "var", "stat")
  out_df <- melt(autocorr, varnames = var_names, value.name = "autocorrelation")
  out_df$lag <- out_df$lag - 1
  
  out_df$facet <- paste(out_df$stat, out_df$var, sep = ": ")
  ncol = length(unique(out_df$var))
  
  ggplot(out_df, aes(x = .data$lag, y = .data$autocorrelation)) +
    geom_bar(stat = "identity", width = 0.3) +
    facet_wrap("facet", ncol = ncol, scales = "free_y")
}

autocorrplot_group <- function(output, lagmax = 40) {
  
  # Get autocorrelation for each model term
  autocorr <- apply(output, c(2,3),
                    function(x) acf(x, lag.max = lagmax, plot = F)$acf)
  
  var_names <- c("lag", "var", "stat")
  out_df <- melt(autocorr, varnames = var_names, value.name = "autocorrelation")
  out_df$lag <- out_df$lag - 1
  
  ggplot(out_df, aes(x = .data$lag, y = .data$autocorrelation)) +
    geom_bar(aes(fill = .data$var), 
             position = "dodge", stat = "identity", width = 0.4) +
    facet_wrap("stat", ncol = 1, scales = "free_y")
}
