#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --job-name=cytospace
#SBATCH --time=12:00:00
#SBATCH -n 10
#SBATCH --mem=150G
#SBATCH --array=1-6

module load anaconda3
cd /cytospace
conda activate cytospace_v1.1.0
module load gcc

# declare sample names
declare -a PAIRS=(
  "O-01_pre-treatment O-01_pre-treatment"
  "O-01_post-treatment O-01_post-treatment"
  "A-03_pre-treatment A-03_pre-treatment"
  "A-03_post-treatment A-03_post-treatment"
  "A-05_pre-treatment A-05_pre-treatment"
  "A-05_post-treatment A-05_post-treatment"
)

# Get the index of the sample pair from the SLURM array task ID
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))
# Extract the sample names for this task
SAMPLE_NAMES=(${PAIRS[$INDEX]})
SAMPLENAME=${SAMPLE_NAMES[0]}
SAMPLENAME2=${SAMPLE_NAMES[1]}
INPUT_DIR="/data/process/cytospace"

cytospace --single-cell \
--scRNA-path ${INPUT_DIR}/${SAMPLENAME}_scrna_counts.txt \
--cell-type-path ${INPUT_DIR}/${SAMPLENAME}_scrna_annos.txt \
--st-path ${INPUT_DIR}/${SAMPLENAME2}_spatial_counts.txt \
--coordinates-path ${INPUT_DIR}/${SAMPLENAME2}_spatial_coords.txt \
--st-cell-type-path ${INPUT_DIR}/${SAMPLENAME2}_spatial_annos.txt \
--output-prefix ${SAMPLENAME2}_ \
--output-folder /analysis/output/cytospace \
--number-of-processors 6 \
--number-of-selected-spots 10000 \
--solver-method lap_CSPR