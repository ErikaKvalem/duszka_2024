library(bigutilsr)
library(ggplot2)
#theme_set(bigstatsr::theme_bigstatsr(0.8))

path ="/data/projects/2024/duszka/NB_project/data/results_counts_raw_tpm_ann/"

file ="Galaxy3603-Column_Join_count_data_tpm.tsv"

counts_tpm <- read.csv(paste0(path,file), sep = "\t")

pca <- prcomp(counts_tpm, scale. = TRUE, rank. = 10)
U <- pca$x

qplot(U[, 1], U[, 2]) + coord_equal()

df_out <- apply(U, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) ))
