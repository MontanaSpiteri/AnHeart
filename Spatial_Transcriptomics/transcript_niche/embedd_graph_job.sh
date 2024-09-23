#!/bin/bash
#SBATCH --job-name=embedd_graph
#SBATCH --mail-type=ALL
#SBATCH --array=0-8
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=250G
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

# Define parameters
dmax="3.0"
nroots=5000
scratch="/graphsage/tmp/subgraph3NB${nroots}"
output_npy="${scratch}/${sample_name}_50_embeddings.npy"
node_meta_min10="${scratch}/${sample_name}_node_meta.csv"
trained_emModel="${scratch}/embmodel/"
full_g="/graphsage/tmp/${sample_name}_dmax${dmax}.pickle"
full_g_csv="/graphsage/tmp/${sample_name}_dmax${dmax}.node_meta.csv"
number_of_samples1=20
number_of_samples2=10
min_comp=10
log_file="/analysis/output/graphsage/${sample_name}_embedd_graph.log"

# Run the Python script
python -u /analysis/scripts/GraphSAGE/embedd_graph.py $sample_name $output_npy $node_meta_min10 $trained_emModel $full_g $full_g_csv $number_of_samples1 $number_of_samples2 $min_comp $log_file
