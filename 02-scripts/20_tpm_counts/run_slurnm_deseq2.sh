#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --output=deseq2.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=990G
#SBATCH --time=30:00:00
#SBTACH --gres=shard:12
sbatch 03_DESeq2_DGEA_neutro.R --count_mat=/data/projects/2022/CRCA/results/v0.1/crc-atlas-dataset/latest/ds_analyses/liana_cell2cell/20_05_tmp/neutrophil_subclusters/tumor_blood/02_pseudobulk/Neutrophil_count_mat.csv --colData=/data/projects/2022/CRCA/results/v0.1/crc-atlas-dataset/latest/ds_analyses/liana_cell2cell/20_05_tmp/neutrophil_subclusters/tumor_blood/02_pseudobulk/Neutrophil_colData.csv --rowData=/data/projects/2022/CRCA/results/v0.1/crc-atlas-dataset/latest/ds_analyses/liana_cell2cell/20_05_tmp/neutrophil_subclusters/tumor_blood/02_pseudobulk/Neutrophil_rowData.csv --prefix=tmp_neutro  --sample_col=sample_col --cond_col=sample_type --sum2zero=FALSE c1="tumor" c2="blood" --cpus=8 covariate_formula = "patient_id +"
