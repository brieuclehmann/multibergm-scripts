generate_networks <- function(formula, n_nodes, coefs, seed = NULL) {
  set.seed(seed)
  
  basis_net <- network(n_nodes, 
                       directed=FALSE, 
                       vertex.attr = list(hemisphere = rep(c("left", "right"), 
                                                           each = n_nodes / 2),
                                          homotopy = rep(1:(n_nodes / 2), 2)))
  
  apply(coefs, 1, function(x) simulate_formula(formula,
                                               coef = x,
                                               basis = basis_net))
}