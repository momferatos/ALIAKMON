#!/bin/bash -l

#-----------------------------------------------------------------
# GPU job , using 80 procs on 4 nodes ,
# with 2 gpus per node, 1 procs per node  and 20 threads per MPI task.
#-----------------------------------------------------------------

#SBATCH --job-name=gpujob # Job name
#SBATCH --partition=gpu # ARIS partition 
#SBATCH --account=pa201102 # Change to your account number
#SBATCH --output=gpujob.out # Stdout (%j expands to jobId)
#SBATCH --error=gpujob.err # Stderr (%j expands to jobId)
#SBATCH --ntasks=2 # Total number of tasks
#SBATCH --gres=gpu:2 # GPUs per node
#SBATCH --nodes=1 # Total number of nodes requested
#SBATCH --ntasks-per-node=2 # Tasks per node
#SBATCH --cpus-per-task=1  # Threads per task
#SBATCH --mem=56000 # Memory per job in MB
#SBATCH -t 00:10:00 

# Load any necessary modules

module purge
module load gnu/8
module load intel/18
module load intelmpi/2018
module load cuda/10.1.168

export LD_LIBRARY_PATH=/users/pr008/gmomfer/Packages/heffte/2.0.0/lib/:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Launch the executable
srun /users/pr008/gmomfer/codes/aliakmon/fft_cuda.exe


