library("conflicted")
library("docopt")
library("BiocParallel")
library("DESeq2")
library("IHW")
library("ggplot2")
library("pcaExplorer")
library("topGO")
library("clusterProfiler")
library("ReactomePA")
library("writexl")
library("readr")
library("dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
library("EnhancedVolcano")
library("ggpubr")
library("tibble")
library("stringr")
library("ggrepel")
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)
library("enrichplot")
organism = "mouse"
gene_id_type = "ENSEMBL"
if (organism == "human") {
  anno_db = "org.Hs.eg.db"
  org_kegg = "hsa"
  org_reactome = "human"
  org_wp = "Homo sapiens"
} else if (organism == "mouse") {
  anno_db = "org.Mm.eg.db"
  org_kegg = "mmu"
  org_reactome = "mouse"
  org_wp = "Mus musculus"
} else {
  msg <- paste0("Organism not implemented: ", organism)
  stop(msg)
}
library(anno_db, character.only = TRUE)


require(data.table)

data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/intestine_CR_Ctrl_IHWallGenes.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/intestine_adlibNB_IHWallGenes.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/intestine_CR_NB_IHWallGenes.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/intestine_NB_IHWallGenes.tsv"))


data <- data_adlibNB
ensg_to_desc = AnnotationDbi::select(get(anno_db), data$gene_id %>% unique(), keytype = gene_id_type, columns = c("SYMBOL")) %>%
  distinct(across(!!gene_id_type), .keep_all = TRUE)

#data  <- data |>
#  left_join(ensg_to_desc, by = c("gene_id" = gene_id_type) ) |>
#  rename(gene_name = SYMBOL) 


##### Pathway enrichment analysis
hgnc_to_entrez = AnnotationDbi::select(get(anno_db), data %>% pull("gene_name") %>% unique(), keytype="SYMBOL", columns=c("ENTREZID"))

# full list with ENTREZIDs added
data_entrez = data %>%  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))
universe = data_entrez %>% pull("ENTREZID") %>% unique()


# list of significant genes with ENTREZIDs added
datasig_fc_entrez <- data %>%  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))
de_foldchanges <- datasig_fc_entrez$log2FoldChange
names(de_foldchanges) <- datasig_fc_entrez$ENTREZID



## ORA
ora_tests = list(
  "KEGG" = function(genes, universe) {
    enrichKEGG(
      gene         = genes,
      universe     = universe,
      organism     = org_kegg,
      pvalueCutoff = 0.05
    )
  },
  "Reactome" = function(genes, universe) {
    enrichPathway(
      gene = genes,
      organism = org_reactome,
      universe = universe,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  },
  "WikiPathway" = function(genes, universe) {
    enrichWP(
      gene = genes,
      universe     = universe,
      organism     = org_wp,
      pvalueCutoff = 0.05
    )
  },
  "GO_BP" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  },
  "GO_MF" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "MF",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  },
  "GO_CC" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = anno_db,
      ont = "CC",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  }
)


# Warmup GO database - work around https://github.com/YuLab-SMU/clusterProfiler/issues/207
._ = enrichGO(universe[1], OrgDb = get(anno_db), keyType = "ENTREZID", ont = "BP", universe = universe)

get_heatplot_dims <- function(p) {
  nr_gene <- length(unique(p$data$Gene))
  nr_cat <- length(unique(p$data$categoryID))
  
  hp_width = min(nr_gene * 0.25, 40)
  hp_height = min(nr_cat * 0.25, 40)
  
  return(c(hp_width, hp_height))
}


ora_name="GO_BP" 
message(paste0("Performing ", ora_name, " ORA-test..."))
results_dir="/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/figures"
prefix="facet_wrap"

test_fun = ora_tests[[ora_name]]
ora_res = test_fun(datasig_fc_entrez$ENTREZID, universe)
#######################
if (!is.null(ora_res)) {
  ora_res = setReadable(ora_res, OrgDb = anno_db, keyType="ENTREZID")
  res_tab = as_tibble(ora_res@result)
  if (min(res_tab$p.adjust) < 0.05 & length(unique(res_tab$geneID)) > 1) {
    # Change the orientation:Horizontal barplot plot
  }  else {
    message(paste0("Warning: No significant enrichment in ", ora_name, " ORA analysis. "))
  
  }
else {
    message(paste0("Warning: No gene can be mapped in ", ora_name, " ORA analysis. "))
}
}
############################
ora_res = setReadable(ora_res, OrgDb = anno_db, keyType="ENTREZID")
res_tab = as_tibble(ora_res@result)
write.csv(res_tab,paste0(results_dir,"/res_tab_adlibNB.csv"))

CR_Ctrl_res <- read.csv("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/figures/res_tab_CR_Ctrl.csv")
adlib_res <-  read.csv("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/figures/res_tab_CR_Ctrl.csv")
CR_Ctrl_res$condition = "CR_Ctrl"
adlib_res$condition = "adlib"
total <- rbind(CR_Ctrl_res, adlib_res)
conflicts_prefer(clusterProfiler::slice)
filtered_data <- total %>%
  group_by(condition) %>%
  slice(1:20)
p <- ggplot(filtered_data, aes(x = Description, y = Count, fill = GeneRatio)) + 
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(condition)) +  coord_flip()

p
p <- ggplot(filtered_data, aes(Count, Description)) + geom_point()
p + facet_grid(rows = vars(condition))
