#!/bin/bash
#SBATCH --job-name=filter_and_subgraph
#SBATCH --mail-type=ALL
#SBATCH --array=0
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=100G
#SBATCH --ntasks-per-node=2
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1

module load gcc/12.2.0
module load cuda/11.8
module load cudnn/11.3.1-8.6.0.162
module load anaconda3
source /apps/anaconda3/anaconda3-2019.03/etc/profile.d/conda.sh
conda activate graphsage

# Define array of samples
samples=(A-01_post-treatment A-03_post-treatment A-03_pre-treatment A-05_post-treatment A-05_pre-treatment O-01_post-treatment O-01_pre-treatment1 O-01_pre-treatment2 O-01_pre-treatment3)

# Extract sample based on SLURM_ARRAY_TASK_ID
sample_name=${samples[$SLURM_ARRAY_TASK_ID]}

# Define other parameters and make directory to store outputs
dmax="3.0"
scratch="/graphsage/tmp"
full_graph="${scratch}/${sample_name}_dmax${dmax}.pickle"
comp_min_n=10  # minimum component size for filtering
nroots=5000  # number of root nodes to sample

mkdir -p ${scratch}/subgraph3NB${nroots}

gpickle="${scratch}/subgraph3NB${nroots}/${sample_name}_dmax${dmax}.subgraph3NB.pickle"
sampledRootsID_csv="${scratch}/subgraph3NB${nroots}/${sample_name}_dmax${dmax}.subgraph3NB.rootNodesID.csv"
log_file="/analysis/output/graphsage/${sample_name}_filter_subgraph.log"

# Run the Python script
python -u /analysis/scripts/GraphSAGE/filter_and_subgraph.py $sample_name $full_graph $comp_min_n $nroots $gpickle $sampledRootsID_csv $log_file