
DECAY = 0.9

output_single: output/single_n10.RDS output/single_n20.RDS output/single_n50.RDS output/single_n100.RDS

output/single_n10.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=10" -v "n_batches=4" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul.sh

output/single_n20.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=20" -v "n_batches=10" -pe smp 10 -l h_rt=2:0:0 scripts/run_simul.sh

output/single_n50.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_batches=10" -pe smp 10 -l h_rt=4:0:0 scripts/run_simul.sh

output/single_n100.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=100" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul.sh

# Continuous covariate

output_cts: output/cts_n10.RDS output/cts_n20.RDS output/cts_n50.RDS

output/cts_n10.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=10" -v "n_batches=4" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_cts.sh

output/cts_n20.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=20" -v "n_batches=10" -pe smp 10 -l h_rt=2:0:0 scripts/run_simul_cts.sh

output/cts_n50.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_batches=10" -pe smp 10 -l h_rt=4:0:0 scripts/run_simul_cts.sh

# Increasing network size

output/single_n50_k60.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_nodes=60" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul.sh

output/single_n50_k90.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_nodes=90" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul.sh

output/single_n50_k120.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_nodes=120" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul.sh

output/single_n50_k150.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_nodes=150" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul.sh

# Two groups
	
output_two: output/twogrp_n10.RDS output/twogrp_n20.RDS output/twogrp_n50.RDS

output/twogrp_n10.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=10" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul_two_group.sh

output/twogrp_n20.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=20" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul_two_group.sh

output/twogrp_n50.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=50" -v "n_batches=10" -pe smp 10 -l h_rt=6:0:0 scripts/run_simul_two_group.sh

# Brain networks

output/brain_n100.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=100" -v "n_batches=20" -pe smp 16 -l h_rt=12:0:0 scripts/fit_brain_nets.sh

output/brain_cts_n100.RDS:
	qsub -v "decay=${DECAY}" -v "n_nets=100" -v "n_batches=20" -pe smp 16 -l h_rt=12:0:0 scripts/fit_brain_nets_age_cts.sh

output/brain_cts_n100_quadratic.RDS:
	qsub -v "decay=${DECAY}" -v "model=quadratic" -v "n_nets=100" -v "n_batches=20" -pe smp 16 -l h_rt=12:0:0 scripts/fit_brain_nets_age_cts.sh

# Prior sensitivity analysis
output/single_n10_prior_mucov10.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "cov_scale=10" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_mucov1.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "cov_scale=1"  -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_mu10_muinf.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "cov_scale=10" -v "mu0=inf" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_mu100_muinf.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "cov_scale=100" -v "mu0=inf" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_mu1_muinf.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "cov_scale=1" -v "mu0=inf" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_nu5.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "nu0=5" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_nu10.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "nu0=10" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh

output/single_n10_prior_nu50.RDS:
	qsub -v "n_nets=10" -v "n_batches=4" -v "nu0=50" -pe smp 4 -l h_rt=2:0:0 scripts/run_simul_prior.sh
