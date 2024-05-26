#!/usr/bin/env Rscript
'
Usage:
03_DESeq2_DGEA.R --count_mat=<count_mat> --colData=<colData> --rowData=<rowData> --prefix=<prefix>  --sample_col=<sample_col> --cond_col=<cond_col> --sum2zero=<sum2zero> --cpus=<cpus> [options]

Mandatory arguments:
  --count_mat=<count_mat>
  --colData=<colData>
  --rowData=<rowData>
  --prefix=<prefix>
  --sample_col=<sample_col>
  --cond_col=<cond_col>


  --sum2zero=<sum2zero>
  --cpus=<cpus>


Optional arguments:
  --covariate_formula=<covariate_formula>
  --c1=<c1>
  --c2=<c2>
  --resDir=<resDir>

' -> doc

# load required packages
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

suppressPackageStartupMessages({
library(BiocParallel)
library(conflicted)
library(readr)
library(tibble)
library(dplyr)
library(stringr)
library(forcats)
library(DESeq2)
library(IHW)
library(limma)
})

# Load parameters
prefix <- arguments$prefix
count_mat <- read_csv(arguments$count_mat)
colData <- read_csv(arguments$colData)
rowData <- read_csv(arguments$rowData)
sample_col <- arguments$sample_col
cond_col <- arguments$cond_col
covariate_formula <- arguments$covariate_formula

sum2zero <-  as.logical(arguments$sum2zero)
if(!sum2zero) {
  c1 <- arguments$c1
  c2 <- arguments$c2
} else {
  c1 = NULL
  c2 = NULL
}

n_cpus <- arguments$cpus

resDir = "/data/projects/2022/CRCA/results/v0.1/crc-atlas-dataset/latest/ds_analyses/liana_cell2cell/core_atlas/tumor_normal/epithelial_cancer/03_deseq2/"

#covariate_formula = ""
design_formula <- as.formula(paste0("~", covariate_formula, " ", cond_col))

register(MulticoreParam(workers = n_cpus))

countData = count_mat |> column_to_rownames(var = "gene_id") |> ceiling()
colData = colData |> column_to_rownames(var = sample_col)



if (length(unique(colData[[cond_col]])) < 2 ) {
  print(paste0("Categories in cond col",length(unique(colData[[cond_col]])) ))
  quit()
}

if(sum2zero) {
  design_mat = model.matrix(design_formula, contrasts.arg = structure(as.list("contr.sum"), names=cond_col), data=colData)  
} else {
  design_mat = design_formula
}

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  rowData = rowData,
  design = design_formula
)
# define reference level (not really necessary when uisng contrasts)
if(!sum2zero) {
dds[[cond_col]] <- relevel(dds[[cond_col]], ref = c2)
} 

## keep only genes where we have >= 10 reads per samplecondition in at least 2 samples
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 20
dds <- dds[keep, ]
#rownames(dds) <- rowData(dds)$GeneSymbol

# save normalized filtered count file
#norm_mat <- counts(dds, normalized = TRUE) |> as_tibble(rownames = "gene_id")
#write_tsv(norm_mat,  paste0(prefix, "_NormalizedCounts.tsv"))

# save normalized batch corrected filtered count file
#vst <- vst(dds, blind = FALSE)
#vst <- vst(dds, nsub=sum( rowMeans( counts(dds, normalized=TRUE)) > 5 ), blind = FALSE) # to avoid this error " less than 'nsub' rows with mean normalized count > 5, " 

#if(!sum2zero) {
#batch <- gsub("\\+", "", covariate_formula) |> str_squish()
#assay(vst) <- limma::removeBatchEffect(x = assay(vst), batch = vst[[batch]])
#write_tsv(assay(vst) |> as_tibble(rownames = "gene_id"), paste0(prefix, "_vst_batch_corrected_NormalizedCounts.tsv"))
#}
# run DESeq
dds <- estimateDispersions(dds)
ddsClean <- dds[which(mcols(dds)$betaConv),]
#dds <- nbinomWaldTest(dds, maxit=500)

dds <- DESeq(dds, parallel = (n_cpus > 1))

if(sum2zero) {
  # order needs to be the one of the levels of the factor (same as for contrast matrix)
  unique_conditions = levels(as.factor(colData[[cond_col]]))
  n_unique = length(unique_conditions)
  # with sum2zero we test that a single coefficient != 0
  # a coefficient corresponds to the difference from the overall mean
  # the intercept correponds to the overall mean
  contr_mat = diag(n_unique - 1) 
  # make list with one contrast per item
  contrasts = lapply(seq_len(n_unique - 1), function(i) { contr_mat[, i] }) 
  # the above added n-1 comparisons, we need to construct the last (redundant) one manually
  contrasts = append(contrasts, list(-apply(contr_mat, MARGIN = 1, sum) / (n_unique - 1)))
  # pad end of vector with zeros (required if there are covariates in the design).
  # we can assume that the "condition columns" always come at the front since
  # it is the first argument of the formula
  contrasts = lapply(contrasts, function(x) {
    c(0, x, rep.int(0, length(resultsNames(dds)) - n_unique))
  })
  # set names of contrasts
  names(contrasts) = unique_conditions
} else {
  contrasts = list(c(cond_col, c1, c2))
  names(contrasts) = sprintf("%s_vs_%s", c1, c2)
}
#
# set names of contrasts
#contrasts <- list(c(cond_col, c1, c2))
#names(contrasts) <- sprintf("%s_vs_%s", c1, c2)

## IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW <- lapply(names(contrasts), function(name) {
  contrast <- contrasts[[name]]
  results(dds, filterFun = ihw, contrast = contrast) |>
    as_tibble(rownames = "gene_id") |>
    mutate(comparison = name) |>
    arrange(pvalue)
}) |> bind_rows()

if(!sum2zero) {
write_tsv(resIHW, paste0(resDir,"Cancer_Epithelial_tumor_vs_normal_DESeq2_result.tsv"))
}else if (sum2zero) {
  write_tsv(resIHW, paste0(resDir,prefix, "_DESeq2_result.tsv"))
}