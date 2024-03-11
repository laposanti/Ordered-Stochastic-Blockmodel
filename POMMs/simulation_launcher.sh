#! /bin/bash -l
#SBATCH --job-name="DED"
# speficity number of nodes 
#SBATCH -N 10

# specify number of tasks/cores per node required
#SBATCH --ntasks-per-node 21

# specify the walltime e.g 20 mins
#SBATCH -t 330:00:00

# set to email at start,end and failed jobs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lapo.santi@ucdconnect.ie

# run from current directory
cd $SLURM_SUBMIT_DIR
module load R/4.3.1

# command to use
Rscript Simulation_launcher.R
 