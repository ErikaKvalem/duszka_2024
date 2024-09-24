#!/bin/bash
#SBATCH --job-name=nb_liver_deseq2
#SBATCH --output=nb_liver_deseq2.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=06:00:00
sbatch 30_runDESeq2_ICBI.R /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/coldata_liver.csv /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/counts_liver.tsv --result_dir=/data/projects/2024/duszka/NB_project/out/40_tpm_counts_240924/liver_NB_WB/CR_Ctrl --c1="CR.Ctrl" --c2="ad.lib.Ctrl" --condition_col="sample_type" --sample_col="sequencingID" --prefix="liver_CR_Ctrl" --organism=mouse --skip_gsea

sbatch 30_runDESeq2_ICBI.R /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/coldata_liver.csv /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/counts_liver.tsv --result_dir=/data/projects/2024/duszka/NB_project/out/40_tpm_counts_240924/liver_NB_WB/adlibNB --c1="ad.lib.NB" --c2="ad.lib.Ctrl" --condition_col="sample_type" --sample_col="sequencingID" --prefix="liver_adlibNB" --organism=mouse --skip_gsea

sbatch 30_runDESeq2_ICBI.R /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/coldata_liver.csv /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/counts_liver.tsv --result_dir=/data/projects/2024/duszka/NB_project/out/40_tpm_counts_240924/liver_NB_WB/CR_NB --c1="CR.NB" --c2="CR.Ctrl" --condition_col="sample_type" --sample_col="sequencingID" --prefix="liver_CR_NB" --organism=mouse --skip_gsea

sbatch 30_runDESeq2_ICBI.R /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/coldata_liver.csv /data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/tables/counts_liver.tsv --result_dir=/data/projects/2024/duszka/NB_project/out/40_tpm_counts_240924/liver_NB_WB/NB --c1="CR.NB" --c2="ad.lib.NB" --condition_col="sample_type" --sample_col="sequencingID" --prefix="liver_NB" --organism=mouse --skip_gsea


