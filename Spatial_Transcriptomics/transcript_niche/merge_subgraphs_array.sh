#!/bin/bash
#SBATCH --job-name=merge_subgraphs
#SBATCH --mail-type=ALL
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

# Define other parameters
dmax="3.0"
nroots=5000
scratch="/graphsage/tmp/subgraph3NB${nroots}"
gpickles=""
rootID_csvs=""
for sample_name in "${samples[@]}"; do
gpickles="${gpickles}${scratch}/${sample_name}_dmax${dmax}.subgraph3NB.pickle,"
rootID_csvs="${rootID_csvs}${scratch}/${sample_name}_dmax${dmax}.subgraph3NB.rootNodesID.csv,"
done
gpickles=${gpickles%,}
rootID_csvs=${rootID_csvs%,}
sg_merge="${scratch}/merged_int_subgraph_sg_dmax${dmax}.pickle"
rootID_reindex="${scratch}/merged_int_subgraph_sg_dmax${dmax}_rootIDs.csv"
log_file="/analysis/output/graphsage/merge_subgraphs.log"

# Run the Python script
python -u /analysis/scripts/GraphSAGE/merge_subgraphs.py $sg_merge $gpickles $rootID_csvs $rootID_reindex $log_file