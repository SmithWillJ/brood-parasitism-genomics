#!/bin/bash -l
#SBATCH --job-name=Test_BayeScEnv
#SBATCH --account=project_xxxx
#SBATCH --output=Test_bayescenv_output.txt
#SBATCH --error=Test_bayescenv_errors.txt

./bayescenv-1.1/bin/linux64/bayescenv test.bsc -env env.txt -o TRIAL_BAYESCENV -nbp 20 -pilot 5000 -thin 10 -n 5000 burn 50000 -pr_jump 0.5 -pr_pref 0.5
