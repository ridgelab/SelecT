#! /bin/bash

#SBATCH --time=5:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=3000M  # memory per CPU core
#SBATCH -J calcStats
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
#SBATCH --mail-user=example@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
cd "$SLURM_SUBMIT_DIR"

workspace_dir="${SLURM_SUBMIT_DIR}/example/SelecT_workspace"
simulation_dir="${SLURM_SUBMIT_DIR}/example/sim"
chromosome=$1
window=$2

time java -Xmx3000m -jar ${SLURM_SUBMIT_DIR}/StatsCalc.jar \
	${workspace_dir} \
	${simulation_dir} \
	${chromosome} \
	${window}

exit 0
