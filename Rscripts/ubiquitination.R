#Load required libraries
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(tibble)
library(forcats)
library(readxl)

#Load pathway collections
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways_hallmark <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

msigdb_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways_go_bp <- split(msigdb_go_bp$gene_symbol, msigdb_go_bp$gs_name)

all_pathways <- list(Hallmark = pathways_hallmark, GO_BP = pathways_go_bp)

#Load the "all_substrates" sheet from your Excel file
substrate_file <- "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/ubiquitination_assessment.xlsx"  # Change if needed
df <- read_excel(substrate_file, sheet = "all_substrates")

#Clean and prepare gene list
df <- df %>%
  filter(!is.na(`Confidence Score`)) %>%
  rename(Gene = `Gene Symbol(Substrate)`) %>%
  mutate(Gene = toupper(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

#Create named score vector
ranks <- df$`Confidence Score`
names(ranks) <- df$Gene
ranks <- sort(ranks, decreasing = TRUE)

#Run fgsea
fgsea_combined_results <- list()

for (category_name in names(all_pathways)) {
  pathways <- all_pathways[[category_name]]
  
  fgsea_res <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 5, maxSize = 500)
  
  if (nrow(fgsea_res) > 0) {
    fgsea_res$Category <- category_name
    fgsea_combined_results[[category_name]] <- fgsea_res}}

#Combine all
fgsea_all <- bind_rows(fgsea_combined_results)

if (nrow(fgsea_all) == 0) stop("No significant fgsea results.")

#Clean and select top results
fgsea_filtered <- fgsea_all %>%
  filter(pval < 0.05) %>%
  mutate(pathway = sub("^HALLMARK_|^GOBP_", "", pathway))


top_fgsea <- fgsea_filtered %>%
  filter(NES > 0) %>%
  group_by(Category) %>%
  slice_min(order_by = pval, n = 15) %>%
  ungroup()

ggplot(top_fgsea, aes(x = NES, y = fct_reorder(pathway, NES), color = pval, size = size)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "P-Value") +
  scale_size_continuous(name = "Genes", labels = scales::number_format(accuracy = 1)) +
  scale_y_discrete(position = "right") +  # Move y-axis text to the right
  facet_grid(rows = vars(Category), scales = "free_y", space = "free", switch = "y") +
  labs(title = "Pathway Enrichment of RBCK1 Ubiquitination Targets",
       x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 12) +
  theme(legend.position = "left",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y.right = element_text(size = 11),  # Add this line
        axis.title.y = element_blank())
