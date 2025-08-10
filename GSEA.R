#Load required libraries
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(forcats)

###Cell subtype DEG GSEA
#Define your subtypes and DEG CSV files
cell_subtypes <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

deg_files <- list("Pre-Unciliated" = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Pre_Unciliated.csv",
                  "Unciliated"     = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Unciliated.csv",
                  "Ciliated"       = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Ciliated.csv",
                  "Secretory"      = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Secretory.csv",
                  "Pre-Ciliated"   = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Pre_Ciliated.csv",
                  "Proliferative"  = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Proliferative.csv")

#Load pathway collections (Hallmark, GO Biological Process, KEGG)
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways_hallmark <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

msigdb_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways_go_bp <- split(msigdb_go_bp$gene_symbol, msigdb_go_bp$gs_name)

#Combine pathway lists in a named list for iteration
all_pathways <- list(Hallmark = pathways_hallmark,
                     GO_BP = pathways_go_bp)

run_fgsea <- function(deg_file, subtype_name, pathways, category_name) {
  deg <- read_csv(deg_file)
  
  deg <- deg %>%
    filter(!is.na(stat)) %>%
    distinct(gene, .keep_all = TRUE)
  
  #Normalize gene names to uppercase
  deg$gene <- toupper(deg$gene)
  pathway_genes <- toupper(unlist(pathways))
  
  deg_filtered <- deg %>% filter(gene %in% pathway_genes)
  
  cat(sprintf("Subtype: %s | Category: %s | Genes in DEG: %d | Filtered genes in pathways: %d\n",
              subtype_name, category_name, nrow(deg), nrow(deg_filtered)))
  
  if (nrow(deg_filtered) == 0) {
    warning(paste("No DEG genes found in pathways for subtype:", subtype_name, "category:", category_name))
    return(tibble())}
  
  ranks <- deg_filtered$stat
  names(ranks) <- deg_filtered$gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  fgsea_res <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 5, maxSize = 500)
  
  fgsea_res_to_save <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = "; ")))
  
  #Save results
  write.csv(fgsea_res_to_save, paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/GSEA/fgsea_", subtype_name, "_", category_name, ".csv"), row.names = FALSE)
  
  #Plot top pathway enrichment if any
  if (nrow(fgsea_res) > 0) {
    top_path <- fgsea_res %>% arrange(padj) %>% slice(1)
    pdf(paste0("fgsea_", subtype_name, "_", category_name, "_plot.pdf"), width = 7, height = 5)
    print(plotEnrichment(pathways[[top_path$pathway]], ranks) + ggtitle(paste(subtype_name, "-", category_name, "-", top_path$pathway)))
    dev.off()}
  
  fgsea_res$Subtype <- subtype_name
  fgsea_res$Category <- category_name
  return(fgsea_res)}

#Run fgsea on each subtype and each pathway category
fgsea_results <- list()

for (category_name in names(all_pathways)) {
  for (subtype in names(deg_files)) {
    cat("Running fgsea for category:", category_name, "subtype:", subtype, "\n")
    res <- run_fgsea(deg_files[[subtype]], subtype, all_pathways[[category_name]], category_name)
    if (nrow(res) > 0) {
      fgsea_results[[paste(subtype, category_name, sep = "_")]] <- res}}}

#Combine all results
fgsea_all <- bind_rows(fgsea_results)

if (nrow(fgsea_all) == 0) {
  stop("No significant fgsea results to plot.")}

#Make sure Subtype is a factor with the desired order
fgsea_all$Subtype <- factor(fgsea_all$Subtype, levels = cell_subtypes)

#Filter again for pval < 0.05 before plotting (just to be safe)
fgsea_all <- fgsea_all %>% filter(pval < 0.05)

#Select top 5 pathways per subtype and category by p-value
top_fgsea <- fgsea_all %>%
  group_by(Subtype, Category) %>%
  slice_min(order_by = pval, n = 6) %>%
  ungroup()

#Remove GOBP_ and HALLMARK_
top_fgsea$pathway <- sub("^GOBP_|^HALLMARK_", "", top_fgsea$pathway)

#Plot top enriched pathways for each subtype and category combined into one plot
ggplot(top_fgsea, aes(x = NES, y = fct_reorder(pathway, NES), color = pval, size = size)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "P-Value") +
  scale_size_continuous(name = "Genes", labels = scales::number_format(accuracy = 1)) +
  scale_y_discrete(position = "right") +  # pathway names on the right
  facet_grid(rows = vars(Category), cols = vars(Subtype), 
             scales = "free_y", space = "free", switch = "y") +  # facet label left
  labs(title = "Top Infertile Enriched Pathways",
       x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 12) +
  theme(legend.position = "left",        # move legend to the left
        strip.placement = "outside",     # move strip text outside the plot panels
        strip.text.y.left = element_text(angle = 0),  # horizontal text
        axis.text.y = element_text(size = 11),
        axis.title.y = element_blank())






####Bulk DEG GSEA
#Load bulk DEG data
bulk_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/bulk_DEGs.csv")

#Load MSigDB pathways
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways_hallmark <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

msigdb_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways_go_bp <- split(msigdb_go_bp$gene_symbol, msigdb_go_bp$gs_name)

all_pathways <- list(Hallmark = pathways_hallmark, GO_BP = pathways_go_bp)

#Define a function to run fgsea
run_fgsea_bulk <- function(deg_df, pathways, category_name) {
  deg_df <- deg_df %>%
    filter(!is.na(stat)) %>%
    distinct(gene, .keep_all = TRUE)
  
  #Normalize gene names
  deg_df$gene <- toupper(deg_df$gene)
  pathway_genes <- toupper(unlist(pathways))
  
  deg_filtered <- deg_df %>% filter(gene %in% pathway_genes)
  
  cat(sprintf("Category: %s | Genes in DEG: %d | Filtered genes in pathways: %d\n",
              category_name, nrow(deg_df), nrow(deg_filtered)))
  
  if (nrow(deg_filtered) == 0) {
    warning(paste("No DEG genes found in pathways for category:", category_name))
    return(tibble())}
  
  ranks <- deg_filtered$stat
  names(ranks) <- deg_filtered$gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  fgsea_res <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 5, maxSize = 500)
  
  fgsea_res_to_save <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = "; ")))
  
  #Save results
  write.csv(fgsea_res_to_save,
            paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/GSEA/fgsea_bulk_", category_name, ".csv"),
            row.names = FALSE)
  
  #Plot top pathway enrichment
  if (nrow(fgsea_res) > 0) {
    top_path <- fgsea_res %>% arrange(padj) %>% slice(1)
    pdf(paste0("fgsea_bulk_", category_name, "_plot.pdf"), width = 7, height = 5)
    print(plotEnrichment(pathways[[top_path$pathway]], ranks) +
            ggtitle(paste("Bulk -", category_name, "-", top_path$pathway)))
    dev.off()}
  
  fgsea_res$Subtype <- "Bulk"
  fgsea_res$Category <- category_name
  return(fgsea_res)}

#Run GSEA on each pathway category
fgsea_results_bulk <- list()

for (category_name in names(all_pathways)) {
  cat("Running fgsea for bulk data - category:", category_name, "\n")
  res <- run_fgsea_bulk(bulk_DE, all_pathways[[category_name]], category_name)
  if (nrow(res) > 0) {
    fgsea_results_bulk[[category_name]] <- res}}

#Combine all results
fgsea_all <- bind_rows(fgsea_results_bulk)

if (nrow(fgsea_all) == 0) stop("No significant fgsea results to plot.")

#Filter for pval < 0.05
fgsea_all <- fgsea_all %>% filter(pval < 0.05)

#Select top 5–6 pathways per category
top_fgsea <- fgsea_all %>%
  group_by(Category) %>%
  slice_min(order_by = pval, n = 20) %>%
  ungroup()

# Clean pathway names
top_fgsea$pathway <- sub("^GOBP_|^HALLMARK_", "", top_fgsea$pathway)

# Plot
ggplot(top_fgsea, aes(x = NES, y = fct_reorder(pathway, NES), color = pval, size = size)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "P-Value") +
  scale_size_continuous(name = "Genes", labels = scales::number_format(accuracy = 1)) +
  scale_y_discrete(position = "right") +
  facet_grid(rows = vars(Category), cols = vars(Subtype), 
             scales = "free_y", space = "free", switch = "y") +
  labs(title = "Top Bulk-Enriched Pathways",
       x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.4),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 11),
        plot.title = element_text(size = 14, face ="bold"))

