## R scripts to generate results in 'Bayesian exponential random graph models for populations of networks'

First, install `multibergm` package:

```
devtools::install_github("brieuclehmann/github")
```

Then, run the following scripts.

`scripts/run_simul.R` simulates a single group of networks and fits the model to this group.

`scripts/run_simul_two_group.R` simulates two groups of networks and fits the model to both groups.

`scripts/fit_brain_nets.R` fits the model to the networks stored in `data/brain_nets.RDS`.

`scripts/make_plots.R` creates the figures in the manuscript.