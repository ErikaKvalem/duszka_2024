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
################### INTESTINE 
resDir = "/data/projects/2024/duszka/NB_project/out/30_other_plots_130624/intestine_NB_WB/"
############################################################################### 
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/GSEA/intestine_CR_Ctrl_GSEA_KEGG.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/GSEA/intestine_adlibNB_GSEA_KEGG.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/GSEA/intestine_CR_NB_GSEA_KEGG.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/GSEA/intestine_NB_GSEA_KEGG.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"


data_CR_Ctrl <- data_CR_Ctrl %>%
  filter(p.adjust < 0.05)
data_adlibNB <- data_adlibNB %>%
  filter(p.adjust < 0.05)
data_CR_NB <- data_CR_NB %>%
  filter(p.adjust < 0.05)
data_NB <- data_NB %>%
  filter(p.adjust < 0.05)


total <- rbind(head(data_CR_Ctrl),head(data_adlibNB),head(data_CR_NB),head(data_NB))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(NES,Description,fill = NES)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "KEGG Pathways", color = "black") # Change the y-axis label to "Pathways"
p2 <- ggplot(total, aes(NES,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "KEGG Pathways", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
p2 <- p2 + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)
print(p2)

# Save the plot as an image file (e.g., PNG)

ggsave(paste0(resDir,"intestine_KEGG_pathways_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

ggsave(paste0(resDir,"intestine_KEGG_pathways_padjust_facet_grid.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed


###################################################################
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/intestine_CR_Ctrl_ORA_GO_BP.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/intestine_adlibNB_ORA_GO_BP.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/intestine_CR_NB_ORA_GO_BP.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/intestine_NB_ORA_GO_BP.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"


#data_CR_Ctrl <- data_CR_Ctrl %>%
#  filter(p.adjust < 0.05)
#data_adlibNB <- data_adlibNB %>%
#  filter(p.adjust < 0.05)
#data_CR_NB <- data_CR_NB %>%
#  filter(p.adjust < 0.05)
#data_NB <- data_NB %>%
#  filter(p.adjust < 0.05)
#

total <- rbind(head(data_CR_Ctrl),head(data_adlibNB),head(data_CR_NB),head(data_NB))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO BP", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"intestine_ORA_GO_BP_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

########################################################
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/intestine_CR_Ctrl_ORA_GO_CC.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/intestine_adlibNB_ORA_GO_CC.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/intestine_CR_NB_ORA_GO_CC.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/intestine_NB_ORA_GO_CC.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"



total <- rbind(head(data_CR_Ctrl),head(data_adlibNB),head(data_CR_NB),head(data_NB))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO CC", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"intestine_ORA_GO_CC_facet_grid.png"), plot = p, width = 6, height = 10)  # Adjust width and height as needed

################################################### INFLAMATION 

data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_Ctrl/intestine_CR_Ctrl_ORA_GO_BP.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/adlibNB/intestine_adlibNB_ORA_GO_BP.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/CR_NB/intestine_CR_NB_ORA_GO_BP.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB/NB/intestine_NB_ORA_GO_BP.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"
total <- rbind(data_CR_Ctrl,data_adlibNB,data_CR_NB,data_NB)
# Function to convert fraction to decimal
fraction_to_decimal <- function(fraction) {
  parts <- strsplit(fraction, "/")[[1]]
  numerator <- as.numeric(parts[1])
  denominator <- as.numeric(parts[2])
  decimal <- numerator / denominator
  return(decimal)
}

# Apply the function to the column
total$GeneRatio <- sapply(total$GeneRatio, fraction_to_decimal)

subset_df <- total[grepl("inflammatory ", total$Description, ignore.case = TRUE), ]

# p plots counts and padjust 
p <- ggplot(subset_df, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO BP", color = "black") # Change the y-axis label to "Pathways"

#p2 plots gene ration , counts and padjust 
p2 <- ggplot(subset_df, aes(GeneRatio,Description)) + 
  geom_point(aes(size = Count, colour = p.adjust)) +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO BP")  + # Change the y-axis label to "Pathways"
  scale_colour_gradient(low = "blue", high = "red")



# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"intestine_ORA_GO_BP_inflammation_counts_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed
ggsave(paste0(resDir,"intestine_ORA_GO_BP_inflammation_generatio_facet_grid.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed

################### LIVER
resDir = "/data/projects/2024/duszka/NB_project/out/30_other_plots_130624/liver_NB_WB/"
############################################################################### 
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_Ctrl/GSEA/liver_CR_Ctrl_GSEA_KEGG.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/adlibNB/GSEA/liver_adlibNB_GSEA_KEGG.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_NB/GSEA/liver_CR_NB_GSEA_KEGG.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/NB/GSEA/liver_NB_GSEA_KEGG.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"

data_CR_Ctrl <- data_CR_Ctrl %>%
  filter(p.adjust < 0.05)
data_adlibNB <- data_adlibNB %>%
  filter(p.adjust < 0.05)
data_CR_NB <- data_CR_NB %>%
  filter(p.adjust < 0.05)
data_NB <- data_NB %>%
  filter(p.adjust < 0.05)


total <- rbind(head(data_CR_Ctrl),head(data_adlibNB),head(data_CR_NB),head(data_NB))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(NES,Description,fill = NES)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "KEGG Pathways", color = "black") # Change the y-axis label to "Pathways"


p2 <- ggplot(total, aes(NES,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "KEGG Pathways", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
p2 <- p2 + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)
print(p2)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"liver_KEGG_pathways_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

ggsave(paste0(resDir,"liver_KEGG_pathways_padjust_facet_grid.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed


###################################################################
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_Ctrl/liver_CR_Ctrl_ORA_GO_BP.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/adlibNB/liver_adlibNB_ORA_GO_BP.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_NB/liver_CR_NB_ORA_GO_BP.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/NB/liver_NB_ORA_GO_BP.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"

#data_CR_Ctrl <- data_CR_Ctrl %>%
#  filter(p.adjust < 0.05)
#data_adlibNB <- data_adlibNB %>%
#  filter(p.adjust < 0.05)
#data_CR_NB <- data_CR_NB %>%
#  filter(p.adjust < 0.05)
#data_NB <- data_NB %>%
#  filter(p.adjust < 0.05)

total <- rbind(head(data_CR_Ctrl),head(data_adlibNB),head(data_CR_NB),head(data_NB))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO BP", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"liver_ORA_GO_BP_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

########################################################
data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_Ctrl/liver_CR_Ctrl_ORA_GO_CC.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/adlibNB/liver_adlibNB_ORA_GO_CC.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_NB/liver_CR_NB_ORA_GO_CC.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/NB/liver_NB_ORA_GO_CC.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"


total <- rbind(head(data_CR_Ctrl),head(data_adlibNB),head(data_CR_NB),head(data_NB))
total$Description <- gsub("- M.*", "", total$Description)

p <- ggplot(total, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO CC", color = "black") # Change the y-axis label to "Pathways"


# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired

# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"liver_ORA_GO_CC_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed

################################################### INFLAMATION 

data_CR_Ctrl<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_Ctrl/liver_CR_Ctrl_ORA_GO_BP.tsv"))
data_adlibNB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/adlibNB/liver_adlibNB_ORA_GO_BP.tsv"))
data_CR_NB <- as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/CR_NB/liver_CR_NB_ORA_GO_BP.tsv"))
data_NB<-as.data.frame(fread("/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/NB/liver_NB_ORA_GO_BP.tsv"))
data_CR_Ctrl$condition = "CR_Ctrl"
data_adlibNB$condition = "adlib"
data_CR_NB$condition = "CR_NB"
data_NB$condition = "NB"



total <- rbind(data_CR_Ctrl,data_adlibNB,data_CR_NB,data_NB)
# Function to convert fraction to decimal
fraction_to_decimal <- function(fraction) {
  parts <- strsplit(fraction, "/")[[1]]
  numerator <- as.numeric(parts[1])
  denominator <- as.numeric(parts[2])
  decimal <- numerator / denominator
  return(decimal)
}

# Apply the function to the column
total$GeneRatio <- sapply(total$GeneRatio, fraction_to_decimal)

subset_df <- total[grepl("inflammatory ", total$Description, ignore.case = TRUE), ]

# p plots counts and padjust 
p <- ggplot(subset_df, aes(Count,Description,fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO BP", color = "black") # Change the y-axis label to "Pathways"

#p2 plots gene ration , counts and padjust 
p2 <- ggplot(subset_df, aes(GeneRatio,Description)) + 
  geom_point(aes(size = Count, colour = p.adjust)) +
  facet_grid(condition ~ .,scales = "free_x", space = "free") +  
  labs(y = "ORA GO BP")  + # Change the y-axis label to "Pathways"
  scale_colour_gradient(low = "blue", high = "red")



# Add continuous color scale
p <- p + scale_fill_gradient(low = "blue", high = "red")
p2 <- p2 + scale_fill_gradient(low = "blue", high = "red")
# Customize color scale for condition

# Add color to facet grid labels
p <- p + theme(strip.text.y = element_text(color = "black"))  # Adjust color as desired
# Print the plot
print(p)

# Save the plot as an image file (e.g., PNG)
ggsave(paste0(resDir,"liver_ORA_GO_BP_inflammation_counts_facet_grid.png"), plot = p, width = 8, height = 10)  # Adjust width and height as needed
ggsave(paste0(resDir,"liver_ORA_GO_BP_inflammation_generatio_facet_grid.png"), plot = p2, width = 8, height = 10)  # Adjust width and height as needed
