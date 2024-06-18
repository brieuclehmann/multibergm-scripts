## R scripts to generate results in 'A Bayesian multilevel model for populations of
networks using exponential-family random graphs'

First, install `multibergm` package:

```
devtools::install_github("brieuclehmann/multibergm")
```

Then, run the following scripts.

`scripts/run_simul.R` simulates a single group of networks and fits the model to this group (c.f. Figures 2-4).  

`scripts/run_simul_cts.R` simulates a single group of networks with a single continuous covariate and fits the model to this group (c.f. Figure 5).

`scripts/run_simul_two_group.R` simulates two groups of networks and fits the model to both groups (c.f. Figure 6).

`scripts/fit_brain_nets.R` fits a two-group model to the networks stored in `data/brain_nets.RDS` (c.f. Figures 7-8).  

`scripts/fit_brain_nets_age_cts.R` fits a model with age and Cattell score to the networks stored in `data/brain_nets.RDS` (c.f. Figures 9).

`scripts/make_simul_plots.R` creates the figures related to the simulated data (Figures 2-6).  

`scripts/make_brain_plots.R` creates the figures related to the brain networks (Figures 7-9).
