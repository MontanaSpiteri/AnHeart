#!/bin/bash
#SBATCH --job-name=train_subgraphs
#SBATCH --mail-type=ALL
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=100G
#SBATCH --ntasks-per-node=2
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A100:1 

module load gcc/12.2.0
module load cuda/11.8
module load cudnn/11.3.1-8.6.0.162
module load anaconda3
source /apps/anaconda3/anaconda3-2019.03/etc/profile.d/conda.sh
conda activate graphsage

# Define parameters
dmax="3.0"
nroots=5000
scratch="/graphsage/tmp/subgraph3NB${nroots}"
number_of_walks=5
num_length=2
# neighbours 20 and 10
number_of_samples1=20
number_of_samples2=10
gpickle="${scratch}/merged_int_subgraph_sg_dmax${dmax}.pickle"
number_epoch=10
trained_model="${scratch}/model"
trained_emModel="${scratch}/embmodel"

mkdir -p ${trained_model}
mkdir -p ${trained_emModel}

rootids="${scratch}/merged_int_subgraph_sg_dmax${dmax}_rootIDs.csv"
log_file="/analysis/output/graphsage/train_subgraphs.log"

# Run the Python script
python -u /analysis/scripts/GraphSAGE/train_subgraphs.py $number_of_walks $num_length $number_of_samples1 $number_of_samples2 $gpickle $number_epoch $trained_model $trained_emModel $rootids $log_file
