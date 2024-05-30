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
 
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/tmp/intestine_CR_Ctrl_GSEA_KEGG.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/tmp/intestine_adlibNB_GSEA_KEGG.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/tmp/intestine_CR_NB_GSEA_KEGG.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/tmp/intestine_NB_GSEA_KEGG.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"

total <- rbind(data_CR_Ctrl, data_adlibNB,data_CR_NB,data_NB)
total$Description <- gsub("-.*", "", total$Description)

conflicts_prefer(clusterProfiler::slice)
filtered_data <- total %>%
  group_by(condition) %>%
  slice(1:5)
p <- ggplot(filtered_data, aes(x = Description, y = NES, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(condition)) +  
  labs(x = "KEGG Pathways", color = "black")  + # Change the y-axis label to "Pathways"
  coord_flip() 

# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/intestine_KEGG_pathways_facet_grid.png", plot = p, width = 8, height = 10)  # Adjust width and height as needed


###################################################################
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/intestine_CR_Ctrl_ORA_GO_BP.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/intestine_adlibNB_ORA_GO_BP.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/intestine_CR_NB_ORA_GO_BP.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/intestine_NB_ORA_GO_BP.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"

total <- rbind(data_CR_Ctrl, data_adlibNB,data_CR_NB,data_NB)
total$Description <- gsub("-.*", "", total$Description)

conflicts_prefer(clusterProfiler::slice)
filtered_data <- total %>%
  group_by(condition) %>%
  slice(1:5)
p <- ggplot(filtered_data, aes(x = Description, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(condition)) +  
  labs(x = "ORA GO BP", color = "black")  + # Change the y-axis label to "Pathways"
  coord_flip() 

# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/intestine_ORA_GO_BP_facet_grid.png", plot = p, width = 8, height = 10)  # Adjust width and height as needed

########################################################
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/intestine_CR_Ctrl_ORA_GO_CC.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/intestine_adlibNB_ORA_GO_CC.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/intestine_CR_NB_ORA_GO_CC.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/intestine_NB_ORA_GO_CC.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"

total <- rbind(data_CR_Ctrl, data_adlibNB,data_CR_NB,data_NB)
total$Description <- gsub("-.*", "", total$Description)

conflicts_prefer(clusterProfiler::slice)
filtered_data <- total %>%
  group_by(condition) %>%
  slice(1:5)
p <- ggplot(filtered_data, aes(x = Description, y = Count, fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(condition)) +  
  labs(x = "ORA GO CC", color = "black")  + # Change the y-axis label to "Pathways"
  coord_flip() 

# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/intestine_ORA_GO_CC_facet_grid.png", plot = p, width = 8, height = 10)  # Adjust width and height as needed

