#! /bin/bash -l
#SBATCH --job-name="DED"
# speficity number of nodes 
#SBATCH -N 1

# specify number of tasks/cores per node required
#SBATCH --ntasks-per-node 24

# specify the walltime e.g 20 mins
#SBATCH -t 330:00:00

# set to email at start,end and failed jobs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lapo.santi@ucdconnect.ie

# run from current directory
cd $SLURM_SUBMIT_DIR
module load R/4.4.0 
module load python/3.3.6

# command to use
git pull
Rscript Model_selectionK3_Simple.R &
Rscript Model_selectionK4_Simple.R &
Rscript Model_selectionK5_Simple.R &
Rscript Model_selectionK6_Simple.R &

wait;





