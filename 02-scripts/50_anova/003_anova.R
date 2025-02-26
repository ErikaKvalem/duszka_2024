# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(conflicted)
#conflicts_prefer(AnnotationDbi::select)
conflicts_prefer(clusterProfiler::select)
conflicts_prefer(clusterProfiler::filter)
#conflicts_prefer(stats::filter)
#conflicts_prefer(lmerTest::lmer)
conflicts_prefer(lme4::lmer)
# Load data
data <- read.csv("../tables/50_anova/SARA2.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Data preprocessing
data <- data[rowSums(is.na(data)) <= 5, ]  # Remove rows with excessive NAs
colnames(data) <- trimws(colnames(data))
data_sub <- select(data, -external_gene_name, -description)  # Remove non-numeric columns
rownames(data_sub) <- data_sub$gene_id  # Set gene_id as rownames
data_sub <- select(data_sub, -gene_id)  # Remove gene_id column
data_sub <- data_sub %>% mutate_all(as.numeric)  # Convert all values to numeric
data_sub <- rownames_to_column(data_sub, var = "gene_id")  # Recreate gene_id column

# Reshape data from wide to long format
data_long <- data_sub %>%
  pivot_longer(cols = -gene_id, names_to = "Sample", values_to = "Expression")

# Create group and time factors
data_long <- data_long %>%
  mutate(
    gene_id = as.character(gene_id),  
    Group = ifelse(str_starts(Sample, "S"), "S", "R"),
    Time = as.character(str_extract(Sample, "[0-9]+"))  # Extract numeric time points as character
  )

data_long$Group <- as.factor(data_long$Group)
data_long$Time <- factor(data_long$Time)

# Create a subset for time points 2 and 3
data_long_2_3 <- data_long %>%
  filter(Time %in% c("2", "3"))

# Initialize results storage
anova_results_all <- data.frame(gene_id = character(), P_value = numeric(), stringsAsFactors = FALSE)
anova_results_2_3 <- data.frame(gene_id = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene and perform two separate analyses
for (gene in unique(data_long$gene_id)) {
  gene_data_all <- filter(data_long, gene_id == gene)  # All time points
  gene_data_2_3 <- filter(data_long_2_3, gene_id == gene)  # Only time points 2 and 3
  
  # Ensure sufficient data before running models
  if (nrow(gene_data_all) > 2 && length(unique(gene_data_all$Time)) > 1 && length(unique(gene_data_all$Group)) > 1) {  
    model_all <- lm(Expression ~ Time * Group, data = gene_data_all)
    anova_res_all <- anova(model_all)
    p_val_all <- if ("Time:Group" %in% rownames(anova_res_all)) anova_res_all["Time:Group", "Pr(>F)"] else NA
    anova_results_all <- rbind(anova_results_all, data.frame(gene_id = gene, P_value = p_val_all))
  }
  
  if (nrow(gene_data_2_3) > 2 && length(unique(gene_data_2_3$Time)) > 1 && length(unique(gene_data_2_3$Group)) > 1) {  
    model_2_3 <- lm(Expression ~ Time * Group, data = gene_data_2_3)
    anova_res_2_3 <- anova(model_2_3)
    p_val_2_3 <- if ("Time:Group" %in% rownames(anova_res_2_3)) anova_res_2_3["Time:Group", "Pr(>F)"] else NA
    anova_results_2_3 <- rbind(anova_results_2_3, data.frame(gene_id = gene, P_value = p_val_2_3))
  }
}

# Adjust p-values using False Discovery Rate (FDR)
anova_results_all$P_adjusted <- p.adjust(na.omit(anova_results_all$P_value), method = "fdr")

adjusted_p_values <- rep(NA, nrow(anova_results_2_3)) 
valid_indices <- which(!is.na(anova_results_2_3$P_value)) 
adjusted_p_values[valid_indices] <- p.adjust(anova_results_2_3$P_value[valid_indices], method = "fdr")
anova_results_2_3$P_adjusted <- adjusted_p_values

gene_name_mapping <- data %>%
  select(gene_id, external_gene_name) %>%
  distinct()  # Ensure unique mapping


anova_results_2_3 <- anova_results_2_3 %>%
  left_join(gene_name_mapping, by = "gene_id") %>%
  rename(gene_name = external_gene_name)  # Rename the column

anova_results_all <- anova_results_all %>%
  left_join(gene_name_mapping, by = "gene_id") %>%
  rename(gene_name = external_gene_name)  # Rename the column

# Save results
write.csv(anova_results_all, "anova_results_all_timepoints_sara2.csv", row.names = FALSE)
write.csv(anova_results_2_3, "anova_results_timepoints_2_3_sara2.csv", row.names = FALSE)
