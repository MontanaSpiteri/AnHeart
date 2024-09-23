#!/bin/bash
#SBATCH --job-name=build_graph
#SBATCH --mail-type=ALL
#SBATCH --array=0-8
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

# Define other parameters
dmax="3.0"
dir="/data/ST"
scratch="/graphsage/tmp"
input_csv="/data_processed/ST/${sample_name}/transcripts.csv.gz"
gene_panel="${dir}/data/meta/xenium_gene_panel.csv"
output_gpickle="${scratch}/${sample_name}_dmax${dmax}.pickle"
output_meta="${scratch}/${sample_name}_dmax${dmax}.node_meta.csv"
output_comps="${scratch}/${sample_name}_dmax${dmax}.components.csv"
log_file="${dir}/analysis/output/graphsage/${sample_name}_build_graph.log"

# Run the Python script
python -u /analysis/scripts/GraphSAGE/build_graph.py $sample_name $input_csv $gene_panel $output_meta $output_comps $dmax $output_gpickle $log_file
