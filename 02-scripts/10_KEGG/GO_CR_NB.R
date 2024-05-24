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

data<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/liver_NB_WB/NB_results_sig_df_ann.tsv"))

ensg_to_desc = AnnotationDbi::select(get(anno_db), data$ensembl_gene_id %>% unique(), keytype = gene_id_type, columns = c("SYMBOL")) %>%
  distinct(across(!!gene_id_type), .keep_all = TRUE)

data  <- data |>
  left_join(ensg_to_desc, by = c("ensembl_gene_id" = gene_id_type) ) |>
  rename(gene_name = SYMBOL) 


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

bplapply(names(ora_tests), function(ora_name) {
  message(paste0("Performing ", ora_name, " ORA-test..."))
  results_dir="/data/scratch/kvalem/projects/2024/duszka/CR_NB/out/liver_NB_WB/NB"
  prefix="NB"
  
  test_fun = ora_tests[[ora_name]]
  ora_res = test_fun(datasig_fc_entrez$ENTREZID, universe)
  
  if (!is.null(ora_res)) {
    ora_res = setReadable(ora_res, OrgDb = anno_db, keyType="ENTREZID")
    res_tab = as_tibble(ora_res@result)
    write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_ORA_", ora_name, ".tsv")))
    
    if (min(res_tab$p.adjust) < 0.05 & length(unique(res_tab$geneID)) > 1) {
      p = dotplot(ora_res, showCategory=30) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))
      
      save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_dotplot")), p, width = 15, height = 10)
      
      p <- cnetplot(ora_res,
                    categorySize="pvalue",
                    showCategory = min(5, length(ora_res@result$ID)),
                    foldChange=de_foldchanges,
                    vertex.label.font=6)
      
      save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_cnetplot")), p, width = 15, height = 12)
      
      p <- heatplot(ora_res, foldChange=de_foldchanges, showCategory=40) +
        scale_fill_gradient2(midpoint=0, low="blue4", mid="white", high="red4" )
      hp_dims <- get_heatplot_dims(p)
      
      save_plot(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_heatplot")), p, width = hp_dims[1], height = hp_dims[2])
      
    } else {
      message(paste0("Warning: No significant enrichment in ", ora_name, " ORA analysis. "))
    }
  } else {
    message(paste0("Warning: No gene can be mapped in ", ora_name, " ORA analysis. "))
  }
})

save_plot <- function(filename, p, width=NULL, height=NULL) {
  if (!is.null(width) && !is.null(height)) {
    ggsave(file.path(paste0(filename, ".png")), plot = p, width = width, height = height)
    ggsave(file.path(paste0(filename, ".svg")), plot = p, width = width, height = height)
  } else {
    ggsave(file.path(paste0(filename, ".png")), plot = p)
    ggsave(file.path(paste0(filename, ".svg")), plot = p)
  }
}

## GSEA
skip_gsea=FALSE

if(!skip_gsea) {
  results_dir="/data/scratch/kvalem/projects/2024/duszka/CR_NB/out/liver_NB_WB/NB"
  prefix="NB"
  # for GSEA use genes ranked by test statistic
  res_ihw_ranked = data_entrez %>%
    arrange(-stat) %>%
    select(ENTREZID, stat) %>%
    na.omit() %>%
    distinct(ENTREZID, .keep_all=TRUE)
  ranked_gene_list = res_ihw_ranked$stat
  names(ranked_gene_list) = res_ihw_ranked$ENTREZID
  
  gsea_tests = list(
    "KEGG"=function(ranked_gene_list) {
      gseKEGG(geneList = ranked_gene_list, organism = org_kegg, pvalueCutoff = 1)
    },
    "Reactome"=function(ranked_gene_list) {
      gsePathway(geneList = ranked_gene_list, organism = org_reactome, pvalueCutoff = 1)
    },
    "WikiPathway"=function(ranked_gene_list) {
      gseWP(geneList = ranked_gene_list, organism = org_wp, pvalueCutoff = 1)
    },
    "GO_BP"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
            keyType = "ENTREZID",
            OrgDb = anno_db,
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            minGSSize = 10)
    },
    "GO_MF"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
            keyType = "ENTREZID",
            OrgDb = anno_db,
            ont = "MF",
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            minGSSize = 10)
    },
    "GO_CC"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
            keyType = "ENTREZID",
            OrgDb = anno_db,
            ont = "CC",
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            minGSSize = 10)
    }
  )
  
  bplapply(names(gsea_tests), function(gsea_name) {
    message(paste0("Performing ", gsea_name, " GSEA-test..."))
    
    test_fun = gsea_tests[[gsea_name]]
    gsea_res = test_fun(ranked_gene_list)
    
    if (!is.null(gsea_res)) {
      gsea_res = setReadable(gsea_res, OrgDb = get(anno_db), keyType="ENTREZID")
      res_tab = gsea_res@result %>% as_tibble()
      
      write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, ".tsv")))
      if (min(res_tab$p.adjust) < 0.05) {
        p = dotplot(gsea_res, showCategory=15, split=".sign") + facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=30))
        
        save_plot(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_dotplot")), p, width = 20, height = 15)
        
        
        p <- cnetplot(gsea_res,
                      categorySize="pvalue",
                      showCategory = 5,
                      foldChange=de_foldchanges,
                      vertex.label.font=6)
        
        save_plot(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_cnetplot")), p, width = 15, height = 12)
        
        # GSEA generates to long gene lists so that the heatplot gets to overloaded
        # p <- heatplot(gsea_res, foldChange=de_foldchanges, showCategory=40) +
        #   scale_fill_gradient2(midpoint=0, low="blue4", mid="white", high="red4" )
        #
        # hp_dims <- get_heatplot_dims(p)
        #
        # ggsave(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_heatplot.png")), plot = p, width = hp_dims[1], height = hp_dims[2])
      } else {
        message(paste0("Warning: No significant enrichment in ", gsea_name, " GSEA analysis. "))
      }
    } else {
      message(paste0("Warning: No gene can be mapped in ", gsea_name, " GSEA analysis. "))
    }
  })
}




