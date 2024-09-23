library(SingleCellExperiment)
library(scran)
library(scater)
library(scuttle)
library(tidyverse)
library(numbat)
library(dplyr)
library(HDF5Array)
library(scistreer)

# Numbat allele preparation (pileup & phasing) from BAM files performed for each participant individually (Pre-treatment & Post-treatment samples combined). Code example below:
#Rscript pileup_and_phase.R \
#--label "A-01" \
#--samples "A-01_pre-treatment","A-01_post-treatment" \
#--bams /SNRNA/A-01_pre-treatment/outs/possorted_genome_bam.bam,/SNRNA/A-01_post-treatment/outs/possorted_genome_bam.bam \
#--barcodes /SNRNA/A-01_pre-treatment/outs/filtered_feature_bc_matrix/barcodes.tsv,/data_processed/SNRNA/A-01_post-treatment/outs/filtered_feature_bc_matrix/barcodes.tsv \
#--outdir /data/process/numbat_prep \
#--gmap /software/external/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
#--eagle /software/external/Eagle_v2.4.1/eagle \
#--snpvcf /ref_db/human/hg38/numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
#--paneldir /ref_db/human/hg38/numbat/1000G_hg38 \
#--ncores 6


# Process Numbat outputs
# Define directories
input_dir <- "/data/process/numbat_prep/"
out_dir <- "/analysis/output/numbat_output/"

# Load the integrated SummarizedExperiment object
sce_all <- loadHDF5SummarizedExperiment("/data/process/", prefix="sce_integrated_")
counts(sce_all) <- as(counts(sce_all), 'dgCMatrix')  

# Define the samples list
samples_list <- list(
  list(patient="A-01", pre="A-01_pre-treatment", post="A-01_post-treatment"),
  list(patient="A-02", pre="A-02_pre-treatment", post="A-02_post-treatment"),
  list(patient="A-03", pre="A-03_pre-treatment", post="A-03_post-treatment"),
  list(patient="A-04", pre="A-04_pre-treatment", post="A-04_post-treatment"),
  list(patient="A-05", pre="A-05_pre-treatment", post="A-05_post-treatment"),
  list(patient="A-06", pre="A-06_pre-treatment", post="A-06_post-treatment"),
  list(patient="A-07", pre="A-07_pre-treatment", post="A-07_post-treatment"),
  list(patient="O-01", pre="O-01_pre-treatment", post="O-01_post-treatment")
  list(patient="O-02", pre="O-02_pre-treatment", post="O-02_post-treatment"),
  list(patient="O-03", pre="O-03_pre-treatment", post="O-03_post-treatment")
)

# Determine the task index from SLURM array
index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Extract the relevant sample pair based on the job array index
sample_pair <- samples_list[[index]]

# Specific operations for the sample pair
patient <- sample_pair$patient
pre <- sample_pair$pre
post <- sample_pair$post

# Subset the SingleCellExperiment object for selected samples
sce <- sce_all[, sce_all$sample %in% c(pre, post)]

# Prepare data for analysis
count_mat <- counts(sce)
count_mat <- count_mat[!duplicated(rownames(count_mat)), ]
count_mat <- count_mat[, !duplicated(colnames(count_mat))]

# Load allele count data specific to each sample
df_allele_pre <- read_tsv(paste0(input_dir, pre, "_allele_counts.tsv"))
df_allele_post <- read_tsv(paste0(input_dir, post, "_allele_counts.tsv"))
df_allele <- rbind(df_allele_pre, df_allele_post)
df_allele <- df_allele[df_allele$cell %in% colnames(count_mat), ]

rm(df_allele_pre, df_allele_post)

# Annotations and cell grouping for analysis
sce$Barcode <- colnames(sce)
cell_annot <- data.frame(cell = sce$Barcode, group = sce$annotation)

# Aggregate counts for reference (Immune cells as reference)
ref_internal <- aggregate_counts(count_mat, cell_annot)
ref_internal <- ref_internal[!duplicated(rownames(ref_internal)), ]
ref_internal <- as.matrix(ref_internal[, "Immune", drop = FALSE])

# Run the analysis using numbat
results <- run_numbat(
  count_mat, 
  ref_internal, 
  df_allele, 
  genome = "hg38",
  t = 1e-5,
  ncores = 8,
  ncores_nni = 8,
  max_entropy = 0.8,
  plot = TRUE,
  out_dir = paste0(out_dir, patient)
)