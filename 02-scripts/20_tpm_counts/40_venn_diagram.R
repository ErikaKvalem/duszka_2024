#install.packages("VennDiagram")
require(data.table)

#intestine_path = "/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/intestine_NB_WB"
#
#adlibNB_deg = read.table(file = paste0(intestine_path,"/adlibNB/intestine_adlibNB_IHWallGenes.tsv"), sep = '\t', header = TRUE)
#CR_Ctrl_deg =  read.table(file = paste0(intestine_path,"/CR_Ctrl/intestine_CR_Ctrl_IHWallGenes.tsv"), sep = '\t', header = TRUE)
#CR_NB_deg =  read.table(file = paste0(intestine_path,"/CR_NB/intestine_CR_NB_IHWallGenes.tsv"), sep = '\t', header = TRUE)
#NB_deg    = read.table(file = paste0(intestine_path,"/NB/intestine_NB_IHWallGenes.tsv"), sep = '\t', header = TRUE)


liver_path = "/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB"

adlibNB_deg = read.table(file = paste0(liver_path,"/adlibNB/liver_adlibNB_IHWallGenes.tsv"), sep = '\t', header = TRUE)
CR_Ctrl_deg =  read.table(file = paste0(liver_path,"/CR_Ctrl/liver_CR_Ctrl_IHWallGenes.tsv"), sep = '\t', header = TRUE)
CR_NB_deg =  read.table(file = paste0(liver_path,"/CR_NB/liver_CR_NB_IHWallGenes.tsv"), sep = '\t', header = TRUE)
NB_deg    = read.table(file = paste0(liver_path,"/NB/liver_NB_IHWallGenes.tsv"), sep = '\t', header = TRUE)



x <- list(
  adlibNB = adlibNB_deg$gene_id, 
  CR_Ctrl = CR_Ctrl_deg$gene_id, 
  CR_NB = CR_NB_deg$gene_id,
  NB = NB_deg$gene_id
)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename =NULL, ...)
  grid.draw(venn_object)
  # Save as PNG
  png(filename="/data/projects/2024/duszka/NB_project/out/20_tpm_counts_260524/liver_NB_WB/venn.png")
  grid.draw(venn_object)
  dev.off()
  
  
}

# Further customization
display_venn(
  x,
  category.names = c("adlib.NBvsadlib.Ctrl" , "CR.Ctrlvsadlib.Ctrl" , "CR.NBvsCR.Ctrl", "CR.NBvsadlib.NB"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1),
)



