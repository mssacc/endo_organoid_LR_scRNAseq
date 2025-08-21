### ====================================================================================================================   
#Endo Organoid Long-Read scRNA-seq Figures
### ====================================================================================================================   

#Load required libraries
library(Seurat)
library(ggplot2)
library(CellChat)
library(biomaRt)
library(clustree)
library(cowplot)
library(curl)
library(dittoSeq)
library(dorothea)
library(dplyr)
library(EnhancedVolcano)
library(fgsea)
library(forcats)
library(ggh4x)
library(ggraph)
library(ggrepel)
library(ggpattern)
library(grid)
library(gridExtra)
library(igraph)
library(IsoformSwitchAnalyzeR)
library(msigdbr)
library(patchwork)
library(pfamAnalyzeR)
library(pheatmap)
library(presto)
library(RColorBrewer)
library(readxl)
library(reshape2)
library(rlang)
library(scales)
library(stringr)
library(tibble)
library(tidygraph)
library(tidyr)
library(UpSetR)
library(VennDiagram)
library(viridis)


### ====================================================================================================================   
#Load seurat object RDS (RNA + iso assays of all 6 samples)
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")


### ====================================================================================================================  
###Figure 1 - General SC Data
#A - Method Outline (Created in BioRender)

#B - Endo organoid UMAP
DimPlot(all_samples_integrated, reduction = "umap", group.by = "cell_type") +
  ggtitle("Endometrial Epithelial Organoid UMAP\n (n = 8,077)") +
  theme(legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) 

#C - UMAP split by fertility  
fertile   <- subset(all_samples_integrated, subset = fertility == "Fertile")
infertile <- subset(all_samples_integrated, subset = fertility == "Infertile")
fer <- DimPlot(fertile, reduction = "umap", group.by = "fertility", cols = "palegreen3") +
          ggtitle("Fertile\n (n = 3,995)") +
          theme(legend.position = "none", 
                plot.title = element_text(size = 16, hjust = 0.5, face = "plain"),  
                axis.line = element_blank(),  
                axis.text = element_blank(), 
                axis.ticks = element_blank(),  
                axis.title = element_blank())  
inf <- DimPlot(fertile, reduction = "umap", group.by = "fertility", cols = "#F8766D") +
          ggtitle("Infertile\n (n = 4,082)") +
          theme(legend.position = "none", 
                plot.title = element_text(size = 16, hjust = 0.5, face = "plain"),
                axis.line = element_blank(), 
                axis.text = element_blank(),  
                axis.ticks = element_blank(),  
                axis.title = element_blank())  
fer + inf


#D - Combination dot plot of all marker genes
DotPlot(all_samples_integrated, group.by = "cell_type",
        features = c("MKI67","PCNA","PCLAF",       # Proliferative
                     "PAEP","DPP4","GADD45A",      # Secretory
                     "FOXJ1","TPPP3","PIFO",       # Ciliated
                     "CCNO","CDC20B","MUC12",      # Pre-Ciliated
                     "CXCL2","TIMP2","SDC4",       # Unciliated
                     "ABCG1","PTGS1","SULT1E1")) + # Pre-Unciliated
  scale_y_discrete(limits = c("Pre-Unciliated","Unciliated","Pre-Ciliated","Ciliated","Secretory","Proliferative")) +
  coord_flip()


#E - Endometrial mapping of SC subypes (made in BioRender)

#F - Glandular and Luminal Marker Gene Expression
#Luminal and glandular marker gene sets
luminal_genes   <- c("PTGS1","VTCN1","CLDN22","IL6","LEFTY1","LGR5","CRISP3","PAX2","SULT1E1")
glandular_genes <- c("FOXA2","MUC1","PAEP","CXCL14","DPP4","GPX3","SCGB2A2","ABCG1","IER3","MUC16","SLPI")

#Compute receptivity score using AddModuleScore
all_samples_integrated <- AddModuleScore(all_samples_integrated, 
                                         features = list(luminal_genes), name = "Luminal")
all_samples_integrated <- AddModuleScore(all_samples_integrated, 
                                         features = list(glandular_genes), name = "Glandular")

#Extract UMAP coordinates and receptivity scores
umap_data <- as.data.frame(Embeddings(all_samples_integrated, reduction = "umap"))
umap_data$Luminal <- all_samples_integrated$Luminal1  # Extract the score from metadata
umap_data$Glandular <- all_samples_integrated$Glandular1
umap_data$cell_type <- all_samples_integrated$cell_type

#Calculate common limits for color scale
score_range <- range(c(umap_data$Luminal, umap_data$Glandular), na.rm = TRUE)

#UMAP of glandular and luminal marker genes with consistent color scale
LE <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = Luminal)) +
          geom_point(size = 2, alpha = 0.8) +
          scale_color_viridis(option = "plasma", limits = score_range) +
          labs(color = "Expression", x = "UMAP 1", y = "UMAP 2", title = "Luminal Marker Gene Expression") +
          theme_minimal(base_size = 14) +
          theme(plot.title = element_text(hjust = 0.5, face="bold", size=20),
                legend.position = "none")
GE <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = Glandular)) +
          geom_point(size = 2, alpha = 0.8) +
          scale_color_viridis(option = "plasma", limits = score_range) +
          labs(color = "Expression", x = "UMAP 1", y = "UMAP 2", title = "Glandular Marker Gene Expression") +
          theme_minimal(base_size = 14) +
          theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))
LE + GE



### ====================================================================================================================   
###Figure 2 - Cell and gene level analysis
#A - Barplot of cell subtype percentages by fertility
#Desired cell type order
cell_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

#Prepare metadata
meta <- as.data.frame(all_samples_integrated@meta.data) %>%
  mutate(fertility = as.character(fertility),
         cell_type = factor(as.character(cell_type), levels = cell_order))

#Count cells and calculate percentages
counts <- meta %>%
  group_by(fertility, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup()

ggplot(counts, aes(x = fertility, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(aes(label = paste0(round(percentage * 100, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 5, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("#F8766D", "#ABA300", "#0CB702", "#00BFC4", "#849AFF", "#FF61CC"),
                    breaks = cell_order) +
  labs(title = "Cell Subtype Percentages", x = NULL, y = "Cell Percentage", fill = "Cell Subtypes") +
  theme_classic() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"))


#B - Differential cell-cell interactions
#Load CellChat object of each dataset and merge them together
cellchat.fer <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat/cellchat_fertile.rds")
cellchat.inf <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat/cellchat_infertile.rds")
object.list <- list(Fertile = cellchat.fer, Infertile = cellchat.inf)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
ptm = Sys.time()
execution.time = Sys.time() - ptm

#Define custom colours
colours.celltype <- c("Pre-Unciliated" = "#F8766D", "Unciliated"    = "#ABA300", 
                      "Ciliated"       = "#0CB702", "Secretory"     = "#00BFC4", 
                      "Pre-Ciliated"   = "#849AFF", "Proliferative" = "#FF61CC")
colours.fertility <- c("Fertile" = "palegreen3", "Infertile" = "#F8766D")

par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", 
                          color.use = colours.celltype, color.edge = c("#F8766D", "palegreen3"))


#C - Compare overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", 
               sources.use = NULL, targets.use = NULL, 
               stacked = TRUE, do.stat = TRUE, color.use = colours.fertility) +
        theme(legend.position = "none",  # Remove legend from gg1
              axis.text.x = element_text(size = 11),
              axis.text.y = element_text(size = 10),
              axis.title.x = element_text(size = 13))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", 
               sources.use = NULL, targets.use = NULL, 
               stacked = FALSE, do.stat = TRUE, color.use = colours.fertility) +
         theme(legend.text = element_text(size = 13),
               axis.text.x = element_text(size = 11),
               axis.text.y = element_text(size = 10),
               axis.title.x = element_text(size = 13),
               legend.position = "bottom")
gg1 + gg2


#D - Up and downregulated genes of each cell subtype
#Load your DEGs
pre_unciliated_df <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Pre_Unciliated.csv")
unciliated_df     <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Unciliated.csv")
ciliated_df       <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Ciliated.csv")
secretory_df      <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Secretory.csv")
pre_ciliated_df   <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Pre_Ciliated.csv")
proliferative_df  <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Proliferative.csv")

#Count up- and downregulated genes
count_regulation <- function(df) {
  up <- sum(df$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(df$log2FoldChange < 0, na.rm = TRUE)
  return(data.frame(Direction = c("Upregulated", "Downregulated"), Count = c(up, down)))}

#Combine all counts
deg_counts <- bind_rows(
  count_regulation(pre_unciliated_df) %>% mutate(CellSubtype = "Pre-Unciliated"),
  count_regulation(unciliated_df)     %>% mutate(CellSubtype = "Unciliated"),
  count_regulation(pre_ciliated_df)   %>% mutate(CellSubtype = "Pre-Ciliated"),
  count_regulation(ciliated_df)       %>% mutate(CellSubtype = "Ciliated"),
  count_regulation(secretory_df)      %>% mutate(CellSubtype = "Secretory"),
  count_regulation(proliferative_df)  %>% mutate(CellSubtype = "Proliferative"))

#Set plotting order
deg_counts$CellSubtype <- factor(deg_counts$CellSubtype, levels = c(
  "Pre-Unciliated", "Unciliated", "Pre-Ciliated", "Ciliated", "Secretory", "Proliferative"))

#Ensure "Upregulated" comes first in bar grouping
deg_counts$Direction <- factor(deg_counts$Direction, levels = c("Upregulated", "Downregulated"))

#Define the custom color scheme for each combination of subtype and gene type
fill_colors <- c("Pre-Unciliated.Upregulated"="#F8766D", "Pre-Unciliated.Downregulated"="#FBBEB5",
                 "Unciliated.Upregulated"    ="#ABA300", "Unciliated.Downregulated"    ="#E1DC85",
                 "Pre-Ciliated.Upregulated"  ="#849AFF", "Pre-Ciliated.Downregulated"  ="#CBD2FF",
                 "Ciliated.Upregulated"      ="#0CB702", "Ciliated.Downregulated"      ="#A8E7A3",
                 "Secretory.Upregulated"     ="#00BFC4", "Secretory.Downregulated"     ="#A1E6E8",
                 "Proliferative.Upregulated" ="#FF61CC", "Proliferative.Downregulated" ="#FFB3E9")

#Add a combined key for color mapping
deg_counts$fill_key <- interaction(deg_counts$CellSubtype, deg_counts$Direction)

ggplot(deg_counts, aes(x = CellSubtype, y = Count, fill = fill_key, pattern = Direction)) +
  geom_bar_pattern(stat = "identity", 
                   position = position_dodge(width = 0.7), width = 0.6, color = "grey",  # bar outline
                   pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.03,
                   pattern_colour = "grey", pattern_fill = "grey") +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.7), vjust = -0.3, size = 4) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = c("Upregulated" = "none", "Downregulated" = "stripe")) +
  theme_minimal(base_size = 14) +
  labs(title = "DEGs per Cell Subtype", y = "Number of DEGs", fill = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


#E - Upset plot of overlapping DEGs
#Load DEGs as character vectors
pre_unciliated_DEGs <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Pre_Unciliated.csv")$gene)
unciliated_DEGs     <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Unciliated.csv")$gene)
ciliated_DEGs       <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Ciliated.csv")$gene)
secretory_DEGs      <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Secretory.csv")$gene)
pre_ciliated_DEGs   <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Pre_Ciliated.csv")$gene)
proliferative_DEGs  <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_Proliferative.csv")$gene)

#Clean list with valid names (no hyphens)
subtype_lists <- list(Pre_Unciliated = pre_unciliated_DEGs,
                      Unciliated     = unciliated_DEGs,
                      Ciliated       = ciliated_DEGs,
                      Secretory      = secretory_DEGs,
                      Pre_Ciliated   = pre_ciliated_DEGs,
                      Proliferative  = proliferative_DEGs)
#Create full gene list (no bulk)
all_genes <- unique(unlist(subtype_lists))

#Create presence matrix
presence_matrix <- data.frame(Gene = all_genes,
                              Pre_Unciliated = all_genes %in% subtype_lists$Pre_Unciliated,
                              Unciliated     = all_genes %in% subtype_lists$Unciliated,
                              Ciliated       = all_genes %in% subtype_lists$Ciliated,
                              Secretory      = all_genes %in% subtype_lists$Secretory,
                              Pre_Ciliated   = all_genes %in% subtype_lists$Pre_Ciliated,
                              Proliferative  = all_genes %in% subtype_lists$Proliferative)
#Prepare matrix for UpSetR
rownames(presence_matrix) <- presence_matrix$Gene
presence_matrix <- presence_matrix[, -1]
numeric_df <- as.data.frame(lapply(presence_matrix, as.integer))
rownames(numeric_df) <- rownames(presence_matrix)

#Original subtype names with hyphens
original_names <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
colnames(numeric_df) <- original_names

#Upset Plot
upset(numeric_df, sets = original_names, nsets = length(original_names),
      nintersects = NA, order.by = "freq",
      mainbar.y.label = "DEG Cross-section", sets.x.label = "Total DEGs in Each Cell Subtype", 
      sets.bar.color = "#849AFF", main.bar.color = "#FF61CC", matrix.color = "#00BFC4",
      text.scale = c(2.5, 2, 1.5, 1.75, 1.75, 2), point.size = 4, line.size = 1)


#F - GSEA analysis
#Define path to GSEA CSVs
gsea_files <- list.files(
  path = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/GSEA/",
  pattern = "^fgsea_.*\\.csv$",
  full.names = TRUE) %>%
  .[!grepl("_bulk_", .) & !grepl("ubiquitination", ., ignore.case = TRUE)]

#Extract Subtype and Category from filename
parse_gsea_filename <- function(filename) {
  parts <- str_match(basename(filename), "fgsea_(.*?)_(.*?)\\.csv$")
  return(list(subtype = parts[,2], category = parts[,3]))}

#Load and combine GSEA files
fgsea_list <- lapply(gsea_files, function(file) {
  info <- parse_gsea_filename(file)
  df <- read_csv(file)
  if (nrow(df) == 0) return(NULL)
  df$Subtype <- info$subtype
  df$Category <- info$category
  return(df)})
fgsea_all <- bind_rows(fgsea_list)

#Define ordered factor levels for cell subtypes
cell_subtypes <- c("Pre-Unciliated","Unciliated","Pre-Ciliated","Ciliated","Secretory","Proliferative")
fgsea_all$Subtype <- factor(fgsea_all$Subtype, levels = cell_subtypes)

#Filter for significant pathways
fgsea_all <- fgsea_all %>% filter(pval < 0.05)

#Keep top 6 pathways per subtype and category
top_fgsea <- fgsea_all %>%
  group_by(Subtype, Category) %>%
  slice_min(order_by = pval, n = 7) %>%
  ungroup()

#Clean pathway names
top_fgsea$pathway <- top_fgsea$pathway %>%
  sub("^GOBP_|^HALLMARK_", "", .)

#Add direction and create ordered factor for plotting
top_fgsea <- top_fgsea %>%
  mutate(Direction = ifelse(NES >= 0, "Upregulated", "Downregulated"),
         Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) %>%
  group_by(Subtype, Category, Direction) %>%
  arrange(desc(NES)) %>%
  mutate(pathway_ordered = factor(pathway, levels = unique(pathway))) %>%
  ungroup()

#Flip the order of the y-axis so highest NES is at the top
top_fgsea$pathway_ordered <- fct_rev(top_fgsea$pathway_ordered)

ggplot(top_fgsea, aes(x = NES, y = pathway_ordered, color = pval, size = size)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "P-Value",
                       labels = scales::label_number(accuracy = 0.001),
                       guide = guide_colorbar(barwidth = 12, barheight = 0.6)) +
  scale_size_continuous(name = "Genes", labels = scales::number_format(accuracy = 1)) +
  scale_y_discrete(position = "right") +
  facet_grid2(rows = vars(Category), cols = vars(Subtype),
              scales = "free_y", space = "free", switch = "y",
              strip = strip_themed(background_x = elem_list_rect(
                fill = c("#F8766D", "#ABA300","#849AFF", "#0CB702", "#00BFC4", "#FF61CC"),
                color = "black"))) +
  labs(title = "Top Infertile Enriched Pathways per Cell Subtype - GSEA",
       x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.4),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 16, face ="bold"))



### ===================================================================================================================
###Figure 3 - Isoform classification, expression and usage
#A - Isoform classes (figure generated in BioRender)

#B - Isoform percentage by isoform classes
#Import SQANTI dataset
SQANTI <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/remove_unknownstrand_classification.txt", sep ='\t')

#Aggregate the expression data by cell type 
counts <- AggregateExpression(all_samples_integrated, assays = "iso", 
                              return.seurat = FALSE, group.by="fertility")
as.data.frame(counts) -> df
row.names(df) -> df$gene

#Split transcript ids into gene and transcript id
pseudobulk_data <- df %>%
  mutate(transcript_id = sub("-.*", "", gene),
         gene_id = sub("^[^-]+-", "", gene))

#Plot the structural category across fertility
merged_data <- merge(pseudobulk_data, SQANTI, 
                     by.x = "transcript_id", by.y = "isoform", all.x = TRUE)

#Pivot pseudobulk data to long format
long_data <- merged_data %>%
  pivot_longer(cols = starts_with("iso."),
               names_to = "fertility", values_to = "expression")

#Generate new df that we can use for plotting attributes defined by fertility
filtered_data <- long_data %>%
  filter(expression > 0) %>%
  distinct()

#Calculate unique counts of genes and isoforms per structural category and fertility
category_summary <- filtered_data %>%
  filter(!is.na(structural_category)) %>%
  group_by(fertility, structural_category, subcategory, coding) %>%
  summarise(num_genes = n_distinct(gene_id),
            num_isoforms = n_distinct(transcript_id),
            .groups = "drop")

#Calculate the total number of isoforms across all categories for each fertility level
category_summary_percent <- category_summary %>%
  group_by(fertility, structural_category) %>%
  summarise(num_isoforms = sum(num_isoforms)) %>%
  ungroup() %>%
  group_by(fertility) %>%
  mutate(percentage = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Define the labeller function for fertility levels
fertility_labeller <- c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")
desired_order <- c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog",
                   "antisense", "fusion", "genic", "genic_intron", "intergenic")
category_summary_percent$structural_category <- factor(category_summary_percent$structural_category, 
                                                       levels = desired_order)
levels(category_summary_percent$structural_category) <- c("FSM", "NIC", "NNC", "Antisense",
                                                          "Fusion","Genic", "Genic Intron", "Intergenic")

ggplot(category_summary_percent, aes(x = fertility, y = percentage, fill = structural_category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme_minimal() +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Isoform Structural Category", y = "Percentage of Isoforms", fill = "Structural Category") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))



#C - Up and downregulated TRANSCRIPTS of each cell subtype
#Load DETs
pre_unciliated_df <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Pre_Unciliated.csv")
unciliated_df     <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Unciliated.csv")
ciliated_df       <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Ciliated.csv")
secretory_df      <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Secretory.csv")
pre_ciliated_df   <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Pre_Ciliated.csv")
proliferative_df  <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Proliferative.csv")

#Count up- and downregulated genes
count_regulation <- function(df) {
  up <- sum(df$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(df$log2FoldChange < 0, na.rm = TRUE)
  return(data.frame(Direction = c("Upregulated", "Downregulated"), Count = c(up, down)))}

#Combine all counts
det_counts <- bind_rows(count_regulation(pre_unciliated_df) %>% mutate(CellSubtype = "Pre-Unciliated"),
                        count_regulation(unciliated_df)     %>% mutate(CellSubtype = "Unciliated"),
                        count_regulation(pre_ciliated_df)   %>% mutate(CellSubtype = "Pre-Ciliated"),
                        count_regulation(ciliated_df)       %>% mutate(CellSubtype = "Ciliated"),
                        count_regulation(secretory_df)      %>% mutate(CellSubtype = "Secretory"),
                        count_regulation(proliferative_df)  %>% mutate(CellSubtype = "Proliferative"))

#Set plotting order
det_counts$CellSubtype <- factor(det_counts$CellSubtype, levels = c(
  "Pre-Unciliated", "Unciliated", "Pre-Ciliated", "Ciliated", "Secretory", "Proliferative"))

#Ensure "Upregulated" comes first in bar grouping
det_counts$Direction <- factor(det_counts$Direction, levels = c("Upregulated", "Downregulated"))

#Define the custom color scheme for each combination of subtype and gene type
fill_colors <- c("Pre-Unciliated.Upregulated"="#F8766D", "Pre-Unciliated.Downregulated"="#FBBEB5",
                 "Unciliated.Upregulated"    ="#ABA300", "Unciliated.Downregulated"    ="#E1DC85",
                 "Pre-Ciliated.Upregulated"  ="#849AFF", "Pre-Ciliated.Downregulated"  ="#CBD2FF",
                 "Ciliated.Upregulated"      ="#0CB702", "Ciliated.Downregulated"      ="#A8E7A3",
                 "Secretory.Upregulated"     ="#00BFC4", "Secretory.Downregulated"     ="#A1E6E8",
                 "Proliferative.Upregulated" ="#FF61CC", "Proliferative.Downregulated" ="#FFB3E9")

#Add a combined key for color mapping
det_counts$fill_key <- interaction(det_counts$CellSubtype, det_counts$Direction)

ggplot(det_counts, aes(x = CellSubtype, y = Count, fill = fill_key, pattern = Direction)) +
  geom_bar_pattern(stat = "identity", 
                   position = position_dodge(width = 0.7), width = 0.6, color = "grey",  # bar outline
                   pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.03,
                   pattern_colour = "grey", pattern_fill = "grey") +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.7), vjust = -0.3, size = 4) +
  scale_fill_manual(values = fill_colors) +
  scale_pattern_manual(values = c("Upregulated" = "none", "Downregulated" = "stripe")) +
  theme_minimal(base_size = 14) +
  labs(title = "DEIs per Cell Subtype", y = "Number of DEIs", fill = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


#D - Upset plot of overlapping DETs
grid.newpage()
#Load DETs as character vectors
pre_unciliated_DETs <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Pre_Unciliated.csv")$transcript)
unciliated_DETs     <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Unciliated.csv")$transcript)
ciliated_DETs       <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Ciliated.csv")$transcript)
secretory_DETs      <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Secretory.csv")$transcript)
pre_ciliated_DETs   <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Pre_Ciliated.csv")$transcript)
proliferative_DETs  <- as.character(read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Proliferative.csv")$transcript)

#Clean list with valid names (no hyphens)
subtype_lists <- list(Pre_Unciliated = pre_unciliated_DETs,
                      Unciliated     = unciliated_DETs,
                      Ciliated       = ciliated_DETs,
                      Secretory      = secretory_DETs,
                      Pre_Ciliated   = pre_ciliated_DETs,
                      Proliferative  = proliferative_DETs)

#Create full gene list (no bulk)
all_transcripts <- unique(unlist(subtype_lists))

#Create presence matrix
presence_matrix <- data.frame(Transcript = all_transcripts,
                              Pre_Unciliated = all_transcripts %in% subtype_lists$Pre_Unciliated,
                              Unciliated     = all_transcripts %in% subtype_lists$Unciliated,
                              Ciliated       = all_transcripts %in% subtype_lists$Ciliated,
                              Secretory      = all_transcripts %in% subtype_lists$Secretory,
                              Pre_Ciliated   = all_transcripts %in% subtype_lists$Pre_Ciliated,
                              Proliferative  = all_transcripts %in% subtype_lists$Proliferative)

#Prepare matrix for UpSetR
rownames(presence_matrix) <- presence_matrix$Transcript
presence_matrix <- presence_matrix[, -1]
numeric_df <- as.data.frame(lapply(presence_matrix, as.integer))
rownames(numeric_df) <- rownames(presence_matrix)
original_names <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
colnames(numeric_df) <- original_names

upset(numeric_df, sets = original_names, nsets = length(original_names),
      nintersects = NA, order.by = "freq",
      mainbar.y.label = "DEI Cross-section", sets.x.label = "Total DEIs in Each Cell Subtype", 
      sets.bar.color = "#849AFF", main.bar.color = "#FF61CC", matrix.color = "#00BFC4",
      text.scale = c(2.5, 2.5, 1.5, 2, 2, 2), point.size = 4, line.size = 1)


#E - Venn diagram of DTU overlap
#Read in the combined dataset from Excel
combined_data <- read_excel("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DTU_final.xlsx")

#Clean transcript IDs (remove version numbers from ENST/Bambu/etc.)
combined_data$Transcript <- gsub("\\.\\d+$", "", combined_data$Transcript)

#Rename columns for consistency
combined_data <- combined_data %>%
  rename(cell_type     = `Cell Subtype`,
         transcript_id = Transcript,
         gene_name     = Gene,
         dataset       = Dataset)

#Define the desired cell type order
desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
combined_data$cell_type <- factor(combined_data$cell_type, levels = desired_order)

#Create sets of "cell_type|transcript_id" strings for each method
dtu_isa     <- combined_data %>%
  filter(dataset == "ISA") %>%
  mutate(pair = paste(cell_type, transcript_id, sep = "|")) %>%
  pull(pair) %>% unique()

dtu_isopod  <- combined_data %>%
  filter(dataset == "Isopod") %>%
  mutate(pair = paste(cell_type, transcript_id, sep = "|")) %>%
  pull(pair) %>% unique()

dtu_dturtle <- combined_data %>%
  filter(dataset == "DTUrtle") %>%
  mutate(pair = paste(cell_type, transcript_id, sep = "|")) %>%
  pull(pair) %>% unique()

#Create named list of these sets
dtu_lists <- list(DTUrtle = dtu_dturtle,
                  IsoPod = dtu_isopod,
                  IsoformSwitchAnalyzeR = dtu_isa)

#Create Venn diagram
venn.plot <- venn.diagram(x = dtu_lists, filename = NULL,
                          fill = c("#00BFC4","#849AFF","#FF61CC"),
                          alpha = 0.6, cex = 1.5,
                          fontfamily = "Arial", cat.fontfamily = "Arial",
                          cat.cex = 1.5,
                          cat.names = c("DTUrtle", "IsoPod", "IsoformSwitchAnalyzeR"),
                          cat.pos = c(-20, 20, 180),     # angles for each label
                          cat.dist = c(0.04, 0.04, 0.04))# distance of labels from the circle
grid.newpage()
grid.draw(venn.plot)


#F - DTU in 2 methods
#Helper function to remove version numbers from ENST IDs
strip_version <- function(id) sub("\\.\\d+$", "", id)

#Load DETs
pre_unciliated_df <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Unciliated.csv")
unciliated_df     <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Unciliated.csv")
ciliated_df       <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Ciliated.csv")
secretory_df      <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Secretory.csv")
pre_ciliated_df   <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Ciliated.csv")
proliferative_df  <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Proliferative.csv")

#Merge DETs
DE_results <- bind_rows(pre_unciliated_df %>% mutate(source = "Pre-Unciliated"),
                        unciliated_df     %>% mutate(source = "Unciliated"),
                        ciliated_df       %>% mutate(source = "Ciliated"),
                        secretory_df      %>% mutate(source = "Secretory"),
                        pre_ciliated_df   %>% mutate(source = "Pre-Ciliated"),
                        proliferative_df  %>% mutate(source = "Proliferative"))

#Split into transcript_id and gene_id, remove version number
DE_results_iso <- DE_results %>%
  separate(transcript, into = c("transcript_id", "gene_id"), sep = "-", extra = "merge") %>%
  mutate(transcript_id = strip_version(transcript_id))

#Load DTU data
file_path <- "/data/gpfs/projects/punim1901/SC_paper/DTU_comparison_new.xlsx"

ISA <- readxl::read_excel(file_path, sheet = "ISA") %>%
  select(cell_type, transcript_id, gene_name) %>%
  distinct() %>%
  mutate(dataset = "ISA", transcript_id = strip_version(transcript_id))

DTUrtle <- readxl::read_excel(file_path, sheet = "DTUrtle") %>%
  select(cell_type, transcript_id, gene_name) %>%
  distinct() %>%
  mutate(dataset = "DTUrtle", transcript_id = strip_version(transcript_id))

isopod <- readxl::read_excel(file_path, sheet = "IsoPod") %>%
  select(cell_type, transcript_id, gene_name) %>%
  distinct() %>%
  mutate(dataset = "isopod", transcript_id = strip_version(transcript_id))

#Create join keys
ISA <- ISA %>% mutate(join_key = paste(cell_type, transcript_id, sep = "_"))
DTUrtle <- DTUrtle %>% mutate(join_key = paste(cell_type, transcript_id, sep = "_"))
isopod <- isopod %>% mutate(join_key = paste(cell_type, transcript_id, sep = "_"))
DE_results_iso <- DE_results_iso %>%
  mutate(join_key = paste(source, transcript_id, sep = "_"))

#Combine DTU results
dtu_combined <- bind_rows(ISA, DTUrtle, isopod) %>%
  select(cell_type, transcript_id, dataset)

#Count number of datasets per (cell_type, transcript_id) pair
dtu_freq <- dtu_combined %>%
  distinct(cell_type, transcript_id, dataset) %>%
  group_by(cell_type, transcript_id) %>%
  summarise(n_methods = n(), .groups = "drop")

#Keep DTUs found in at least 2 datasets
DTU_DET <- dtu_freq %>%
  filter(n_methods >= 2)

#Filter DE_results_iso and Join Transcript Names
filtered_DE <- DE_results_iso %>%
  inner_join(DTU_DET, by = c("transcript_id", "source" = "cell_type"))

#Join to get full transcript names
filtered_DE <- filtered_DE %>%
  left_join(DE_results %>%
              mutate(transcript_id_match = strip_version(str_extract(transcript, "^[^-]+"))) %>%
              select(transcript_id_match, source, transcript),
            by = c("transcript_id" = "transcript_id_match", "source"))

#Prepare Heatmap Matrix
logFC_mat_df <- filtered_DE %>%
  group_by(transcript, source) %>%
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

logFC_wide <- logFC_mat_df %>%
  pivot_wider(names_from = source, values_from = log2FoldChange, values_fill = 0) %>%
  column_to_rownames(var = "transcript")

#Reorder columns
desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
logFC_wide <- logFC_wide[, intersect(desired_order, colnames(logFC_wide))]

#Define your custom transcript order (to appear at the top)
custom_order <- c("ENST00000383329.7-HLA-C", "ENST00000376228.9-HLA-C",
                  "BambuTx1337-PSORS1C3",    "BambuTx1338-PSORS1C3",
                  "ENST00000356286.9-RBCK1", "ENST00000475269.5-RBCK1",
                  "ENST00000422689.3-ZNF880")

#Extract gene names from transcript names
all_transcripts <- rownames(logFC_wide)
gene_names <- str_extract(all_transcripts, "(?<=-).+$")
avg_log2FC <- rowMeans(logFC_wide, na.rm = TRUE)

#Create a data frame of all transcripts
ordering_df <- data.frame(
  transcript = all_transcripts,
  gene = gene_names,
  avg_log2FC = avg_log2FC)

#Split into custom and remaining transcripts
custom_df <- ordering_df %>% filter(transcript %in% custom_order) %>%
  mutate(order = match(transcript, custom_order)) %>%
  arrange(order)

remaining_df <- ordering_df %>% filter(!transcript %in% custom_order)

#Order the remaining transcripts by gene (grouped) and avg expression (increasing)
remaining_df <- remaining_df %>%
  arrange(gene, avg_log2FC)

#Combine the two sets
final_order <- c(custom_df$transcript, remaining_df$transcript)

#Reorder the matrix
logFC_wide <- logFC_wide[final_order, ]

#Remove one iso that keeps coming up even though it's only in 1 DTU method
logFC_wide <- logFC_wide[!rownames(logFC_wide) %in% "ENST00000567998.5-SULT1A1", ]

my_breaks <- seq(-8, 8, length.out = 100)
legend_ticks <- seq(-8, 8, by = 2)

pheatmap::pheatmap(logFC_wide, show_rownames = TRUE,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   treeheight_col = 0, treeheight_row = 0,
                   angle_col = 0,
                   breaks = my_breaks, legend_breaks = legend_ticks,
                   main = "Isoform Expression of DTUs Found in ≥2 Methods")



### ====================================================================================================================
###Figure 4 - HLA-C
#Set colours  
cluster_colours <- c("#F8766D","#ABA300","#0CB702","#00BFC4","#849AFF","#FF61CC")  # Customize for cluster_names

#A - Isoform expression split by fertility
#Expression of a specific isoform split by fertility
split_plot <- FeaturePlot(all_samples_integrated, reduction = "umap",
                          features = "ENST00000383329.7-HLA-C", split.by = "fertility") &
                    theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 16))

#Create a dummy plot to extract legend
legend_plot <- FeaturePlot(all_samples_integrated, reduction = "umap", features = "ENST00000383329.7-HLA-C") +
  labs(color = "Log-normalized\nisoform expression")
legend_only <- get_legend(legend_plot)

#Combine
wrap_plots(split_plot, wrap_elements(legend_only), ncol = 2, widths = c(3, 0.6)) +
  plot_annotation(title = "ENST00000383329.7-HLA-C") &
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 16))

#Vln plot of isoform split by fertility
VlnPlot(all_samples_integrated, features = "ENST00000383329.7-HLA-C",
        split.by = "cell_type", group.by = "fertility", cols = cluster_colours, assay = "iso") +
  ggtitle("ENST00000383329.7-HLA-C") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold", size = 16),
        axis.text.x = element_text(angle = 0, hjust = 0.5))


#B - Gene expression split by fertility
#Expression of a gene split by fertility
split_plot <- FeaturePlot(all_samples_integrated, reduction = "umap",
                          features = "HLA-C", split.by = "fertility") &
                    theme(plot.title = element_text(hjust = 0.45, face = "bold", size = 16))

#Create a dummy plot to extract legend
legend_plot <- FeaturePlot(all_samples_integrated, reduction = "umap", features = "HLA-C") +
  labs(color = "Log-normalized\ngene expression")
legend_only <- get_legend(legend_plot)

#Combine
wrap_plots(split_plot, wrap_elements(legend_only), ncol = 2, widths = c(3, 0.6)) +
  plot_annotation(title = "HLA-C") &
  theme(plot.title = element_text(hjust = 0.45, face = "bold", size = 16))


#VlnPlot of gene split by fertility
VlnPlot(all_samples_integrated, features = "HLA-C",
        split.by="cell_type", group.by="fertility", cols=cluster_colours) +
  ggtitle("HLA-C") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 0, hjust = 0.5))


#C - Isoform expression percentages
#HLA-C isoforms
gene_isoforms <- c("iso_ENST00000376228.9-HLA-C", "iso_ENST00000383329.7-HLA-C", 
                   "iso_ENST00000487245.5-HLA-C", "iso_ENST00000376237.8-HLA-C", 
                   "iso_ENST00000470363.5-HLA-C", "iso_ENST00000466892.5-HLA-C", 
                   "iso_ENST00000484378.1-HLA-C", "iso_ENST00000495835.1-HLA-C", 
                   "iso_ENST00000640219.1-HLA-C")

#Extract expression levels of each isoform along with fertility metadata
isoform_counts <- FetchData(all_samples_integrated, vars = c(gene_isoforms, "fertility"))

#Sum across cells separately for Fertile and Infertile samples
isoform_sums <- isoform_counts %>%
  group_by(fertility) %>%
  summarise(across(all_of(gene_isoforms), sum, na.rm = TRUE))

#Compute isoform percentages separately for each fertility group
isoform_percentages <- isoform_sums %>%
  pivot_longer(cols = -fertility, names_to = "isoform", values_to = "count") %>%
  group_by(fertility) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

#Remove "iso_" prefix from isoform names in both the data and the levels list
isoform_percentages$isoform <- gsub("^iso_", "", isoform_percentages$isoform)
gene_isoforms_clean <- gsub("^iso_", "", gene_isoforms)

#Make sure isoform is a factor with correct levels
isoform_percentages$isoform <- factor(isoform_percentages$isoform, levels = gene_isoforms_clean)

#Recalculate average isoform percentages across both fertility groups
isoform_means <- isoform_percentages %>%
  group_by(isoform) %>%
  summarise(mean_pct = mean(percentage), .groups = "drop")

#Define which isoforms to group as "Other"
low_abundance_isoforms <- isoform_means %>%
  filter(mean_pct < 0.5) %>%
  pull(isoform)

#Replace those isoforms with "Other" and re-aggregate
isoform_percentages_grouped <- isoform_percentages %>%
  mutate(isoform = ifelse(isoform %in% low_abundance_isoforms, "Other Isoforms", as.character(isoform))) %>%
  group_by(fertility, isoform) %>%
  summarise(percentage = sum(percentage), .groups = "drop")

#Ensure factor levels reflect only what's in the final dataset (important for fill color mapping)
isoform_percentages_grouped$isoform <- factor(isoform_percentages_grouped$isoform,
                                              levels = unique(c(
                                                setdiff(gene_isoforms_clean, low_abundance_isoforms),
                                                "Other Isoforms")))

ggplot(isoform_percentages_grouped, aes(x = "", y = percentage, fill = isoform)) +
  geom_bar(stat = "identity", width = 1, colour = "black", size = 0.1) +
  coord_polar(theta = "y", start = pi / 6) +
  facet_wrap(~fertility) +
  theme_void() +
  labs(title = "HLA-C Isoform Expression", fill = "HLA-C Isoforms") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 2)),
        strip.text = element_text(size = 15, margin = margin(b = 0, t = 10)),  #Reduce vertical space in facet strip
        strip.placement = "outside",
        strip.background = element_blank(),  # Removes grey background padding
        plot.margin = margin(t = 5, r = 10, b = 10, l = 10),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11))


#D - Isoform switch analysis
#Load cell subtype IsoformSwitchAnalyzeR rds
pre_unciliated_DTU <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR/SwitchListAnalyzed_Pre-Unciliated.rds")
unciliated_DTU <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR/SwitchListAnalyzed_Unciliated.rds")
ciliated_DTU <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR/SwitchListAnalyzed_Ciliated.rds")
secretory_DTU <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR/SwitchListAnalyzed_Secretory.rds")
pre_ciliated_DTU <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR/SwitchListAnalyzed_Pre-Ciliated.rds")
proliferative_DTU <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/IsoformSwitchAnalyzeR/SwitchListAnalyzed_Proliferative.rds")

#Custom bar colors
custom_colors <- c("Fertile" = "palegreen3", "Infertile" = "#F8766D")

#Generate isoform expression plots
p1 <- switchPlotIsoExp(unciliated_DTU, gene = "HLA-C",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE, alphas = c(0.05, 0.05),
                       localTheme = theme_bw(base_size = 13))
p2 <- switchPlotIsoExp(ciliated_DTU, gene = "HLA-C",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE, alphas = c(0.05, 0.05),
                       localTheme = theme_bw(base_size = 13))
p3 <- switchPlotIsoExp(secretory_DTU, gene = "HLA-C",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE, alphas = c(0.05, 0.05),
                       localTheme = theme_bw(base_size = 13))
p4 <- switchPlotIsoExp(proliferative_DTU, gene = "HLA-C",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE, alphas = c(0.05, 0.05),
                       localTheme = theme_bw(base_size = 13))

#Generate isoform usage plots
p5 <- switchPlotIsoUsage(unciliated_DTU, gene = "HLA-C",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13))
p6 <- switchPlotIsoUsage(ciliated_DTU, gene = "HLA-C",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13))
p7 <- switchPlotIsoUsage(secretory_DTU, gene = "HLA-C",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13))
p8 <- switchPlotIsoUsage(proliferative_DTU, gene = "HLA-C",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13))

#Combine into a list
plot_list <- list(p1, p2, p3, p4, p5, p6, p7, p8)

#Apply styling and color changes
plots_cleaned <- lapply(seq_along(plot_list), function(i) {
  p <- plot_list[[i]]
  base_plot <- p +
    scale_fill_manual(values = custom_colors) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(t = 5, r = 30, b = 5, l = 5))
  #Remove x-axis text/ticks for top row only
  if (i <= 4) {
    base_plot + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  } else {
    base_plot}})

#Reassign cleaned plots
p1 <- plots_cleaned[[1]]; p2 <- plots_cleaned[[2]]; p3 <- plots_cleaned[[3]]
p4 <- plots_cleaned[[4]]; p5 <- plots_cleaned[[5]]; p6 <- plots_cleaned[[6]]
p7 <- plots_cleaned[[7]]; p8 <- plots_cleaned[[8]]

#Combine into rows and final layout
row1 <- p1 | p2 | p3 | p4 
row2 <- p5 | p6 | p7 | p8
(row1 / row2) + plot_layout(guides = "collect") & theme(legend.position = "right")

#E - HLA-C IsoVis
#F - HLA-C AlphaFold



### ====================================================================================================================
###Figure 5 - Transcription factors
data(dorothea_hs, package = "dorothea")
net <- dorothea_hs                   
net = dorothea_hs %>%
  filter(confidence %in% c("A","B","C"))
net <- net[,c('tf','target','mor','confidence')]
colnames(net) <- c('source','target','mor','confidence')

#Load in DE analysis and add a cell_type column
pre_unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Pre_Unciliated.csv") %>%
  mutate(cell_type = "Pre-Unciliated")
unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Unciliated.csv") %>%
  mutate(cell_type = "Unciliated")
ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Ciliated.csv") %>%
  mutate(cell_type = "Ciliated")
secretory_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Secretory.csv") %>%
  mutate(cell_type = "Secretory")
pre_ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Pre_Ciliated.csv") %>%
  mutate(cell_type = "Pre-Ciliated")
proliferative_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Proliferative.csv") %>%
  mutate(cell_type = "Proliferative")

#Combine all into one data frame
combined_DE <- bind_rows(pre_unciliated_DE, unciliated_DE, ciliated_DE,
                         secretory_DE, pre_ciliated_DE, proliferative_DE)
results_2 <- combined_DE
results_2 <- results_2 %>%
  mutate(sig = ifelse(abs(log2FoldChange) >= 0.585 & padj < 0.05, TRUE, FALSE))

#Get TFs that are DE and in the DoRothEA regulon
DE_TFs <- results_2 %>%
  filter(gene %in% net$source, sig) %>%
  select(gene, cell_type, log2FoldChange, padj)

#Join regulon with DE TF info
regulon_with_DE_TFs <- net %>%
  inner_join(DE_TFs, by = c("source" = "gene")) %>%
  rename(TF = source, TF_log2FC = log2FoldChange, TF_padj = padj, TF_cell_type = cell_type)

#Now join to find if the target gene is also a DEG in the same cell type
TF_target_DEG_pairs <- regulon_with_DE_TFs %>%
  inner_join(results_2 %>% filter(sig) %>%
               select(gene, log2FoldChange, padj, cell_type),
             by = c("target" = "gene", "TF_cell_type" = "cell_type")) %>%
  rename(target_log2FC = log2FoldChange, target_padj = padj, cell_type = TF_cell_type)

#Heatmap of TFs
target_TFs <- c("POU2F2", "EGR1", "POU5F1", "BCL11A", "MEIS2", "FOXM1", "SOX9", "TFAP2A")
tf_matrix <- combined_DE %>%
  filter(gene %in% target_TFs) %>%
  select(gene, cell_type, log2FoldChange) %>%
  pivot_wider(names_from = cell_type, values_from = log2FoldChange, values_fill = 0) %>%
  column_to_rownames("gene")

desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
tf_matrix <- tf_matrix[, intersect(desired_order, colnames(tf_matrix))]
breaks_continuous <- seq(-2, 2, by = 0.1)
legend_ticks <- seq(-2, 2, by = 1)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaks_continuous) - 1)

pheatmap(as.matrix(tf_matrix),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Gene Log2FC of Differentially Expressed TFs",
         fontsize_row = 10, fontsize_col = 11, angle_col = 0,
         treeheight_row = 0, treeheight_col = 0,
         breaks = breaks_continuous, legend_breaks = legend_ticks, color = my_colors)


#B - TF Networks
#Determine global log2FC limits
global_log2fc_vals <- TF_target_DEG_pairs %>%
  select(TF_log2FC, target_log2FC) %>%
  unlist() %>%
  as.numeric()
log2fc_limit <- max(abs(global_log2fc_vals), na.rm = TRUE)
log2fc_breaks <- seq(-log2fc_limit, log2fc_limit, by = 1)

#Network plot function with fixed scale
plot_tf_network <- function(cell_type_name) {
  network_df <- TF_target_DEG_pairs %>%
    filter(cell_type == cell_type_name) %>%
    select(TF, target, mor, TF_log2FC, target_log2FC)
  
  if (nrow(network_df) == 0) return(NULL)
  tf_nodes <- network_df %>%
    select(name = TF, log2FC = TF_log2FC) %>%
    mutate(node_type = "TF")
  target_nodes <- network_df %>%
    select(name = target, log2FC = target_log2FC) %>%
    mutate(node_type = "Target")
  node_attrs <- bind_rows(tf_nodes, target_nodes) %>%
    group_by(name) %>%
    slice(1) %>%
    ungroup()
  edges <- network_df %>%
    mutate(mor = factor(mor, levels = c("-1", "1"))) %>%
    select(from = TF, to = target, mor)
  graph <- tbl_graph(nodes = node_attrs, edges = edges, directed = TRUE)
  
  ggraph(graph, layout = "fr") +
    geom_edge_link(aes(color = factor(mor)),
                   arrow = arrow(length = unit(3, 'mm')),
                   end_cap = circle(3, 'mm'), width = 1.2) +
    geom_node_point(aes(fill = log2FC), shape = 21, size = 6, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4,
                   box.padding = 0.2, point.padding = 0.1, max.overlaps = 15) +
    scale_fill_gradientn(colours = my_colors, limits = c(-2.5, 2.5), breaks = legend_ticks,
                         labels = scales::number_format(accuracy = 1), name = "DEG log2FC") +
    scale_edge_color_manual(values = c("-1" = "#F8766D", "1" = "palegreen3"),
                            labels = c("-1" = "Repression", "1" = "Activation"),
                            name = "Regulation") +
    theme_graph() +
    labs(title = paste("TF-target network in\n", cell_type_name, "cells")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_text(face = "bold"),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "bottom")}

p1 <- plot_tf_network("Pre-Unciliated") + theme(plot.margin = margin(5, 5, 5, 5))
p2 <- plot_tf_network("Unciliated") + theme(plot.margin = margin(5, 5, 5, 5))
p3 <- plot_tf_network("Ciliated") + theme(plot.margin = margin(5, 5, 5, 5))
p4 <- plot_tf_network("Secretory") + theme(plot.margin = margin(5, 5, 5, 5))

(p1 + p2 + p3 + p4) + plot_layout(guides = "collect", nrow = 1) & 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12, hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5))


#C + D - TF iso
#Load in DE analysis and add a cell_type column
pre_unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Unciliated.csv") %>%
  mutate(cell_type = "Pre-Unciliated")
unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Unciliated.csv") %>%
  mutate(cell_type = "Unciliated")
ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Ciliated.csv") %>%
  mutate(cell_type = "Ciliated")
secretory_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Secretory.csv") %>%
  mutate(cell_type = "Secretory")
pre_ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Ciliated.csv") %>%
  mutate(cell_type = "Pre-Ciliated")
proliferative_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Proliferative.csv") %>%
  mutate(cell_type = "Proliferative")

#Combine all into one data frame
combined_DE <- bind_rows(pre_unciliated_DE, unciliated_DE, ciliated_DE,
                         secretory_DE, pre_ciliated_DE, proliferative_DE)
combined_DE <- combined_DE %>%
  tidyr::separate(transcript, into = c("isoform", "gene"), sep = "-", extra = "merge", remove = FALSE) %>%
  mutate(sig = ifelse(abs(log2FoldChange) >= 0.585 & padj < 0.05, TRUE, FALSE))

#Print out DEIs that are also TFs
results <- combined_DE
results$sig <- ifelse(abs(results$log2FoldChange) >=0.585 & results$padj<0.05,TRUE,FALSE)
DEGTFs <- results %>%
  filter(gene %in% net$source & sig)

#Pull out DEIs that are targets
DEG_targets <- results %>%
  filter(gene %in% net$target & sig)

#Create combined dataframe of common TF and TARGET DEIs per cell subtype
results_2 <- combined_DE
results_2 <- results_2 %>%
  mutate(sig = ifelse(abs(log2FoldChange) >= 0.585 & padj < 0.05, TRUE, FALSE))

#Get TFs that are DE and in the DoRothEA regulon
DE_TFs <- results_2 %>%
  filter(gene %in% net$source, sig) %>%
  select(gene, cell_type, log2FoldChange, padj)

#Join regulon with DE TF info
regulon_with_DE_TFs <- net %>%
  inner_join(DE_TFs, by = c("source" = "gene")) %>%
  rename(TF = source, TF_log2FC = log2FoldChange, TF_padj = padj, TF_cell_type = cell_type)

#Now join to find if the target gene is also a DEI in the same cell type
TF_target_DEG_pairs <- regulon_with_DE_TFs %>%
  inner_join(results_2 %>% filter(sig) %>%
               select(gene, log2FoldChange, padj, cell_type),
             by = c("target" = "gene", "TF_cell_type" = "cell_type")) %>%
  rename(target_log2FC = log2FoldChange,
         target_padj = padj,
         cell_type = TF_cell_type)


#C - Heatmap of iso TFs
#Create a wide matrix of log2FC per TF per cell type
desired_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

#Set cell_type as factor with all levels
DEGTFs$cell_type <- factor(DEGTFs$cell_type, levels = desired_order)

tf_matrix <- DEGTFs %>%
  mutate(iso_gene = paste(isoform, gene, sep = "-")) %>%
  select(iso_gene, cell_type, log2FoldChange) %>%
  pivot_wider(names_from = cell_type, values_from = log2FoldChange, values_fill = 0) %>%
  column_to_rownames("iso_gene")

#Add missing columns with zeros
missing_cols <- setdiff(desired_order, colnames(tf_matrix))
for (col in missing_cols) {
  tf_matrix[[col]] <- 0}

#Reorder columns
tf_matrix <- tf_matrix[, desired_order]

breaks_continuous <- seq(-3, 3, by = 0.1)
legend_ticks <- seq(-3, 3, by = 1)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaks_continuous) - 1)

pheatmap(as.matrix(tf_matrix),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "Isoform Log2FC of Differentially Expressed TFs",
         fontsize_row = 10, fontsize_col = 11, angle_col = 0,
         treeheight_row = 0, treeheight_col = 0,
         breaks = breaks_continuous, legend_breaks = legend_ticks, color = my_colors)


#D - Network of iso TFs
#Define colors
my_colors <- c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#ffffbf", 
               "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026")
legend_ticks <- c(-2.5, -1.25, 0, 1.25, 2.5)

#Filter DEGs and build isoform-gene names
results_2 <- results_2 %>%
  mutate(sig = ifelse(abs(log2FoldChange) >= 0.585 & padj < 0.05, TRUE, FALSE),
         iso_gene = paste(isoform, gene, sep = "-"))

#Get DE TFs
DE_TFs <- results_2 %>%
  filter(gene %in% net$source, sig) %>%
  select(gene, isoform, cell_type, log2FoldChange, padj) %>%
  mutate(iso_gene = paste(isoform, gene, sep = "-"))

#Join with regulon network
regulon_with_DE_TFs <- net %>%
  inner_join(DE_TFs, by = c("source" = "gene")) %>%
  rename(TF = source, TF_log2FC = log2FoldChange, TF_padj = padj, TF_cell_type = cell_type, TF_iso_gene = iso_gene)

#Get DE targets
target_DEGs <- results_2 %>%
  filter(sig) %>%
  select(gene, isoform, log2FoldChange, padj, cell_type) %>%
  mutate(iso_gene = paste(isoform, gene, sep = "-"))

#Final joined table: TF → target (both at isoform-gene level)
TF_target_DEG_pairs <- regulon_with_DE_TFs %>%
  inner_join(target_DEGs, by = c("target" = "gene", "TF_cell_type" = "cell_type")) %>%
  rename(target_log2FC = log2FoldChange,
         target_padj = padj,
         target_iso_gene = iso_gene,
         cell_type = TF_cell_type)

#Function to plot network per cell type
plot_tf_network <- function(cell_type_name) {
  network_df <- TF_target_DEG_pairs %>%
    filter(cell_type == cell_type_name) %>%
    select(TF_iso_gene, target_iso_gene, mor, TF_log2FC, target_log2FC)
  
  if (nrow(network_df) == 0) return(NULL)
  
  tf_nodes <- network_df %>%
    select(name = TF_iso_gene, log2FC = TF_log2FC) %>%
    mutate(node_type = "TF")
  
  target_nodes <- network_df %>%
    select(name = target_iso_gene, log2FC = target_log2FC) %>%
    mutate(node_type = "Target")
  
  node_attrs <- bind_rows(tf_nodes, target_nodes) %>%
    group_by(name) %>%
    slice(1) %>%
    ungroup()
  
  edges <- network_df %>%
    mutate(mor = factor(mor, levels = c("-1", "1"))) %>%
    select(from = TF_iso_gene, to = target_iso_gene, mor)
  
  graph <- tbl_graph(nodes = node_attrs, edges = edges, directed = TRUE)
  
  ggraph(graph, layout = "fr") +
    geom_edge_link(aes(color = factor(mor)),
                   arrow = arrow(length = unit(3, 'mm')),
                   end_cap = circle(3, 'mm'), width = 1.2) +
    geom_node_point(aes(fill = log2FC), shape = 21, size = 6, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    scale_fill_gradientn(colours = my_colors, limits = c(-2.5, 2.5), breaks = legend_ticks,
                         labels = scales::number_format(accuracy = 1), name = "log2FC") +
    scale_edge_color_manual(values = c("-1" = "#F8766D", "1" = "palegreen3"),
                            labels = c("-1" = "Repression", "1" = "Activation"),
                            name = "Regulation", drop = FALSE) +
    theme_graph() +
    labs(title = paste("TF-target network in\n", cell_type_name, "cells")) +
    theme(plot.title = element_text(hjust = 0.5))}

#Create plots
cell_types <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
plots <- lapply(cell_types, plot_tf_network)

#Replace NULLs with blank plots for consistent layout
plots <- lapply(plots, function(p) if (is.null(p)) ggplot() + theme_void() else p)

#Combine plots with patchwork
wrap_plots(plots, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, hjust = 0.5))  # smaller title size


#E - ZNF880 ISA Plotting
#Custom colors for Fertile and Infertile
custom_colors <- c("Fertile" = "palegreen3", "Infertile" = "#F8766D")

#Generate isoform expression plot
IE <- switchPlotIsoExp(ciliated_DTU, gene = "ZNF880",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE, alphas = c(0.05, 0.05),
                       localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors) +
  theme(plot.margin = margin(10, 10, 10, 10))

#Generate isoform usage plot
IU <- switchPlotIsoUsage(ciliated_DTU, gene = "ZNF880",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors) +
  theme(plot.margin = margin(10, 10, 10, 10))

#Combine plots side by side with shared legend on the right
(IE | IU) + plot_layout(guides = "collect") & theme(legend.position = "right")

#F - ZNF880 IsoVis



### ====================================================================================================================
###Figure 6 - - RBCK1
#A - RBCK1 ISA Plotting
#Define custom colors
custom_colors <- c("Fertile" = "palegreen3", "Infertile" = "#F8766D")

#Generate plots
p1 <- switchPlotGeneExp(pre_unciliated_DTU, gene = "RBCK1",
                        condition1 = "Fertile", condition2 = "Infertile",
                        addErrorbars = FALSE,
                        localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors)
p2 <- switchPlotIsoExp(pre_unciliated_DTU, gene = "RBCK1",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE,
                       localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors)
p3 <- switchPlotIsoUsage(pre_unciliated_DTU, gene = "RBCK1",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors)
p4 <- switchPlotGeneExp(secretory_DTU, gene = "RBCK1",
                        condition1 = "Fertile", condition2 = "Infertile",
                        addErrorbars = FALSE, 
                        localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors)
p5 <- switchPlotIsoExp(secretory_DTU, gene = "RBCK1",
                       condition1 = "Fertile", condition2 = "Infertile",
                       addErrorbars = FALSE, 
                       localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors)
p6 <- switchPlotIsoUsage(secretory_DTU, gene = "RBCK1",
                         condition1 = "Fertile", condition2 = "Infertile",
                         localTheme = theme_bw(base_size = 13)) +
  scale_fill_manual(values = custom_colors)

#Extract plot data
df_gene <- p1$data

#Rebuild the gene expression plot manually with full control
p1_fixed <- ggplot(df_gene, aes(x = Condition, y = gene_expression, fill = Condition)) +
              geom_bar(stat = "identity", position = "dodge") +
              scale_fill_manual(values = custom_colors) +
              theme_bw(base_size = 13) +
              theme(legend.position = "none",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_blank(),
                    plot.margin = margin(5, 50, 5, 5))

df_gene_p4 <- p4$data

p4_fixed <- ggplot(df_gene_p4, aes(x = Condition, y = gene_expression, fill = Condition)) +
              geom_bar(stat = "identity", position = "dodge") +
              scale_fill_manual(values = custom_colors) +
              theme_bw(base_size = 13) +
              theme(legend.position = "none",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    plot.title = element_blank(),
                    plot.margin = margin(5, 50, 5, 5))

plots_cleaned <- lapply(list(p4_fixed, p5, p6, p1_fixed, p2, p3), function(p) {
  p + theme(legend.position = "none",
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(5, 30, 5, 5))})

plots_cleaned_mod <- lapply(plots_cleaned, function(p) {
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))})

wrap_plots(plots_cleaned_mod, nrow = 1, widths = c(1.4,1.5,1.5, 1.4,2.2,2.2))  + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

#B - RBCK1 Isovis
#C - RBCK1 AlphaFold

#D - Ubiquitination GSEA
#Load pathway collections
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways_hallmark <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)

msigdb_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
pathways_go_bp <- split(msigdb_go_bp$gene_symbol, msigdb_go_bp$gs_name)

all_pathways <- list(Hallmark = pathways_hallmark, GO_BP = pathways_go_bp)

#Load the "all_substrates" sheet from Excel file
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
  theme(legend.position = "right",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y.right = element_text(size = 11),  # Add this line
        axis.title.y = element_blank())

                
#E - NOTCH signalling (generated in BioRender)

                

### ====================================================================================================================
### ====================================================================================================================
###Supplementary Figures
### ====================================================================================================================
### ====================================================================================================================
###Supplementary Figure 1
#A -  Clustree plot
all_samples_clustree <- FindClusters(all_samples_integrated, reduction = "umap", resolution = seq(0.05, 0.6, 0.05), dims = 1:10)
clustree <- clustree(all_samples_clustree, prefix = "RNA_snn_res.", node_size_range = c(5, 15))
clustree + theme(axis.text = element_text(size = 0), axis.title = element_text(size = 0),    
                 legend.text = element_text(size = 12),
                 legend.title = element_text(size = 13)) +
  labs(edge_alpha = "Proportion", colour = "Resolution", size = "Node Size")

#B - Epithelial cell marker genes
FeaturePlot(all_samples_integrated, reduction="umap", features=c("EPCAM","CDH1","KRT8"), ncol=3)

#C - Stem cell marker genes
FeaturePlot(all_samples_integrated, reduction="umap", features=c("SOX9","AXIN2","CDH2","FUT4"), ncol=4)



### ====================================================================================================================
###Supplementary Figure 2
#A - David Analysis of clusters
#Import all cell DAVID analysis datasets
file_path <- "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/david_analysis/cluster_pathways.xlsx"

#Get the names of all sheets
sheet_names <- excel_sheets(file_path)
sheets_list <- lapply(sheet_names, function(sheet) {
  read_excel(path = file_path, sheet = sheet)})
names(sheets_list) <- sheet_names

#Define 6 upregulated cell subtypes only
subtypes <- list("Pre-Unciliated" = sheets_list[["pre_unciliated_up"]],
                 "Unciliated"     = sheets_list[["unciliated_up"]],
                 "Ciliated"       = sheets_list[["ciliated_up"]],
                 "Secretory"      = sheets_list[["secretory_up"]],
                 "Pre-Ciliated"   = sheets_list[["pre_ciliated_up"]],
                 "Proliferative"  = sheets_list[["proliferative_up"]])

#Loop over subtypes and plot
for (i in seq_along(subtypes)) {cell_type <- names(subtypes)[i]
up_df <- subtypes[[i]]
if (nrow(up_df) == 0) next  # Skip if empty
up_df$Direction <- "Up"

combined_df <- up_df %>%
  filter(Pvalue < 0.05)

combined_df$log10P <- -log10(combined_df$Pvalue)

combined_df <- combined_df %>%
  group_by(Term, Direction) %>%
  slice(1) %>%
  ungroup()

#Order terms
combined_df <- combined_df %>%
  arrange(Category, desc(log10P))
combined_df$Term <- factor(combined_df$Term, levels = rev(unique(combined_df$Term)))

#Define the desired order for categories
category_levels <- c("Biological Process", "Cellular Component", "Molecular Function", "KEGG Pathway")

#Recode and set as factor
combined_df$Category <- recode(combined_df$Category,
                          "GOTERM_BP_DIRECT" = "Biological Process",
                          "GOTERM_CC_DIRECT" = "Cellular Component",
                          "GOTERM_MF_DIRECT" = "Molecular Function",
                          "KEGG_PATHWAY"     = "KEGG Pathway")
combined_df$Category <- factor(combined_df$Category, levels = category_levels)
custom_colors <- c("Biological Process" = "#849AFF", "Cellular Component" = "#CBD2FF",
                   "Molecular Function" = "#00C094", "KEGG Pathway"     = "#53B400")

p <- ggplot(combined_df, aes(x = log10P, y = Term, fill = Category)) +
  geom_col(width = 0.6) +
  geom_point(aes(size = Count), shape = 21, color = "black", stroke = 0.3) +
  scale_size_continuous(range = c(5, 10)) +
  geom_text(aes(label = Count), size = 3, color = "black", vjust = 0.5, fontface = "bold") +
  scale_x_continuous(name = expression(-log[10](P[adj]))) +
  scale_fill_manual(values = custom_colors) +
  labs(title = paste(cell_type, "Cluster"),
       y = NULL, fill = "Category", size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold"))

#Save plot
assign(paste0("plot_", i), p)}

#Combine and display all 6 plots with shared legend
legend <- get_legend(plot_2 + theme(legend.position = "right"))
combined_plot <- plot_grid(plot_1 + theme(legend.position = "none"),
                           plot_2 + theme(legend.position = "none"),
                           plot_3 + theme(legend.position = "none"),
                           plot_4 + theme(legend.position = "none"),
                           plot_5 + theme(legend.position = "none"),
                           plot_6 + theme(legend.position = "none"), ncol = 3)
plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(0.85, 0.15))


#B - Cell Cycle Phase
seurat_obj <- all_samples_integrated
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes,
                               set.ident = TRUE)
colours <- c("G1" = "#849AFF", "S" = "#0CB702", "G2M" = "#FF61CC")
DimPlot(seurat_obj, group.by = "Phase", reduction = "umap", cols = colours) + ggtitle("Cell Cycle Phase")


#C - Cell subtype percentages by sample
cell_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

#Prepare metadata
meta <- as.data.frame(all_samples_integrated@meta.data) %>%
  mutate(sample = as.character(sample),
         cell_type = factor(as.character(cell_type), levels = cell_order))

#Count cells and calculate percentages
counts <- meta %>%
  group_by(sample, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup()

ggplot(counts, aes(x = sample, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_text(aes(label = paste0(round(percentage * 100, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 5, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("#F8766D", "#ABA300", "#0CB702", "#00BFC4", "#849AFF", "#FF61CC"),
                    breaks = cell_order) +
  labs(title = "Cell Subtype Percentages by Sample", x = NULL, y = "Cell Percentage", fill = "Cell Subtypes") +
  theme_classic() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"))
                


### ====================================================================================================================
###Supplementary Figure 3 - Cell-cell signalling
colours.celltype <- c("#F8766D","#ABA300","#0CB702","#00BFC4","#849AFF","#FF61CC")
                
#A - Cell populations with significant changes in sending or receiving signals
#Compute centrality for both objects
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]])
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]])

#Calculate weights
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))  # control dot size

p1 <- netAnalysis_signalingRole_scatter(object.list[[1]],
                                        title = names(object.list)[1],
                                        color.use = colours.celltype,
                                        weight.MinMax = weight.MinMax, label.size = 0) +
  coord_cartesian(xlim = c(3, 7), ylim = c(2, 8)) +
  theme_bw() + theme(text = element_text(size = 16),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 15),
                     plot.title = element_text(size = 18, face = "bold"),
                     panel.grid.minor = element_line(color = "grey80", size = 0.3),
                     legend.position = "none")
p2 <- netAnalysis_signalingRole_scatter(object.list[[2]],
                                        title = names(object.list)[2],
                                        color.use = colours.celltype,
                                        weight.MinMax = weight.MinMax, label.size = 0) +
  coord_cartesian(xlim = c(3, 7), ylim = c(2, 8)) +
  theme_bw() + theme(text = element_text(size = 16),
                     axis.text = element_text(size = 14),
                     axis.title = element_text(size = 15),
                     plot.title = element_text(size = 18, face = "bold"),
                     panel.grid.minor = element_line(color = "grey80", size = 0.3))
patchwork::wrap_plots(p1, p2) 


#B - Identify the signaling changes of specific cell populations
gg0 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Pre-Unciliated",
                                            color.use = c("black", "palegreen3","#F8766D"))
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Unciliated",
                                            color.use = c("black", "palegreen3","#F8766D"))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Ciliated",
                                            color.use = c("black", "palegreen3","#F8766D"))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Secretory",
                                            color.use = c("black", "palegreen3","#F8766D"))
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Pre-Ciliated",
                                            color.use = c("black", "palegreen3","#F8766D"))
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Proliferative",
                                            color.use = c("black", "palegreen3","#F8766D"))

gg0 <- gg0 + theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title = element_text(size = 11),
                   axis.text = element_text(size = 10)) + labs(title = "Signalling Changes in Pre-Unciliated")
gg1 <- gg1 + theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title = element_text(size = 11),
                   axis.text = element_text(size = 10)) + labs(title = "Signalling Changes in Unciliated")
gg2 <- gg2 + theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title = element_text(size = 11),
                   axis.text = element_text(size = 10)) + labs(title = "Signalling Changes in Ciliated")
gg3 <- gg3 + theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title = element_text(size = 11),
                   axis.text = element_text(size = 10)) + labs(title = "Signalling Changes in Secretory")
gg4 <- gg4 + theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title = element_text(size = 11),
                   axis.text = element_text(size = 10)) + labs(title = "Signalling Changes in Pre-Ciliated")
gg5 <- gg5 + theme(legend.position = "none",plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.title = element_text(size = 11),
                   axis.text = element_text(size = 10)) + labs(title = "Signalling Changes in Proliferative")

#Visualizing differential outgoing and incoming signaling changes from fertile to infertile
patchwork::wrap_plots(list(gg0, gg1, gg2, gg3, gg4, gg5), ncol = 3) +
  theme(legend.position = "none")


#C - SPP1 chord plot
group.names <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
color.use <- c("Pre-Unciliated" = "#F8766D", "Unciliated"   = "#ABA300", "Ciliated"     = "#0CB702",
               "Secretory"      = "#00BFC4", "Pre-Ciliated" = "#849AFF","Proliferative" = "#FF61CC")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = "SPP1")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = "SPP1", 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste("SPP1", names(object.list)[i]),
                      color.use = color.use)}


#D - SPP1 LR graph (Fer vs Inf)
#Redo merge with shortened suffixes for bubble plot
#Load CellChat object of each dataset and merge them together
cellchat.fer.short <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_fertil_short.rds")
cellchat.inf.short <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_infertile_short.rds")
object.list <- list(Fertile = cellchat.fer.short, Infertile = cellchat.inf.short)
cellchat_short <- mergeCellChat(object.list, add.names = c("F", "I"))
ptm = Sys.time()
execution.time = Sys.time() - ptm

#Define colours
colours.fertility <- c("F" = "palegreen3", "I" = "#F8766D")

netVisual_bubble(cellchat_short,
                 sources.use = c(1,2,4,6), targets.use = c(1:6),
                 signaling = c("SPP1"), comparison = c(1, 2),
                 angle.x = 45, font.size = 8,
                 color.text = colours.fertility,
                 title.name = "SPP1 Signalling Pathways L-R Interactions between Cell Subtypes Divided by Fertility", 
                 font.size.title = 20) +
  theme(axis.text.x = element_text(size = 11.5),  # Increase x-axis text size
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),  # Increase legend text size
        legend.title = element_text(size = 12, face = "bold"), # Increase legend title size
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom",
        plot.title = element_text(margin = margin(b = 30), hjust = 1.2, face = "bold")) +
  guides(color = guide_colorbar(
                  title = "Communication Prob", barwidth = 9, barheight = 0.7, title.vjust = 1.5),
         size = guide_legend(
                  title = "P-value", override.aes = list(shape = 21, fill = "black", color = "black"),
                  label.position = "right", title.vjust = 1.5))



### ====================================================================================================================
###Supplementary Figure 4 - Gene Comparison
#Aggregate the expression data by cell type 
counts <- AggregateExpression(all_samples_integrated, assays = "iso", 
                              return.seurat = FALSE, group.by="fertility")
as.data.frame(counts) -> df
row.names(df) -> df$gene

#Split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "-",  extra = "merge") 

#Select cell type columns and gather into long format
long_data <- pseudobulk_data %>%
  pivot_longer(cols = starts_with("iso."),
               names_to = "fertility", values_to = "expression")

#Filter for non-zero expression to identify genes and isoforms expressed in each cell type
filtered_data <- long_data %>%
  filter(expression > 0) %>%
  dplyr::select(fertility, transcript_id, gene_id) %>%
  distinct()

#Calculate unique gene and isoform counts for each cell type
cell_type_summary <- filtered_data %>%
  group_by(fertility) %>%
  summarise(num_genes = n_distinct(gene_id),
            num_isoforms = n_distinct(transcript_id)) %>%
  pivot_longer(cols = c(num_genes, num_isoforms), 
               names_to = "Category", values_to = "Count")

#A - Plot gene and isoform count
ggplot(cell_type_summary, aes(x = fertility, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  scale_fill_manual(values = c("num_genes" = "#849AFF", "num_isoforms" = "#CBD2FF"),
                    labels = c("num_genes" = "Genes", "num_isoforms" = "Isoforms")) +
  labs(title = "Total Number of\n Genes and Isoforms", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#B - Novel genes
#Filter novel genes
novel_genes <- pseudobulk_data %>%
  filter(grepl("^BambuGene\\d+$", gene_id))

#Determine detection status in each condition
novel_status <- novel_genes %>%
  mutate(detected_fertile = iso.Fertile > 0,
         detected_infertile = iso.Infertile > 0) %>%
  select(gene_id, detected_fertile, detected_infertile) %>%
  as_tibble()

#Classify genes by detection status
gene_classification <- novel_status %>%
  mutate(status = case_when(detected_fertile & detected_infertile ~ "Common",
                            detected_fertile & !detected_infertile ~ "Unique to Fertile",
                            !detected_fertile & detected_infertile ~ "Unique to Infertile",
                            TRUE ~ NA_character_)) %>%
                  filter(!is.na(status))

#Reshape to long format for plotting: one row per (gene, fertility condition)
plot_data <- gene_classification %>%
  pivot_longer(
    cols = c(detected_fertile, detected_infertile),
    names_to = "fertility",
    values_to = "detected") %>%
  filter(detected) %>%
  mutate(fertility = ifelse(fertility == "detected_fertile", "Fertile", "Infertile")) %>%
  mutate(category = case_when(status == "Common" ~ "Common",
                              status == "Unique to Fertile" & fertility == "Fertile" ~ "Unique to Fertile",
                              status == "Unique to Infertile" & fertility == "Infertile" ~ "Unique to Infertile",
                              TRUE ~ NA_character_)) %>%
                      filter(!is.na(category))

#Count genes per fertility and category
counts <- plot_data %>%
  group_by(fertility, category) %>%
  summarise(n = n_distinct(gene_id), .groups = "drop") %>%
  mutate(category = factor(category, levels = c("Common", "Unique to Fertile", "Unique to Infertile")))

#Calculate percentages within each fertility group
counts <- counts %>%
  group_by(fertility) %>%
  mutate(percent = n / sum(n) * 100,
         ypos = cumsum(n) - 0.5 * n) %>%  # for label positioning in the middle of the bar segment
  ungroup()

#Reverse factor levels to put unique on top in the stacked bar
counts$category <- factor(counts$category, levels = rev(levels(counts$category)))

#Plot novel genes
ggplot(counts, aes(x = fertility, y = n, fill = category)) +
  geom_bar(stat = "identity", width = 0.5) +   #narrowed bars 
  geom_text(aes(y = ypos, label = paste0(round(percent, 1), "%")), color = "black", size = 4) +
  scale_fill_manual(values = c("Unique to Fertile" = "palegreen3", "Unique to Infertile" = "#F8766D",
                               "Common" = "#CBD2FF")) +
  labs(title = "Novel Genes",  x = "Fertility Condition", y = "Number of Novel Genes", fill = "Category") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#C - Bulk DEG volcano plot
#Load in bulk DE analysis
bulk_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/bulk_DE_analysis.csv")

#Define colours for volcano plot
keyvals <- ifelse(bulk_DE$log2FoldChange < -0.585 & bulk_DE$padj < 0.05, '#F8766D',  
                  ifelse(bulk_DE$log2FoldChange > 0.585 & bulk_DE$padj < 0.05,  'palegreen3',  
                         '#717171'))  
keyvals[is.na(keyvals)] <- '#717171'
names(keyvals)[keyvals ==  'palegreen3'] <- 'Up in Infertile'
names(keyvals)[keyvals ==  '#F8766D'] <-    'Down in Infertile'
names(keyvals)[keyvals ==  '#717171'] <-    'Not Significant'

p <- EnhancedVolcano(bulk_DE, lab = bulk_DE$gene,
                     x = 'log2FoldChange', y = 'padj',
                     pCutoff = 0.05, FCcutoff = 0.585,
                     pointSize = 5, labSize = 4, axisLabSize = 14, 
                     xlim = c(-7.5, 8.5),
                     title = "Global DEGs - Infertile vs Fertile", titleLabSize = 20,
                     subtitle = NULL, colCustom = keyvals,
                     drawConnectors = TRUE, max.overlaps = 20, arrowheads = FALSE,
                     legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5)
p + theme(plot.margin = margin(10, 10, 10, 10), 
          legend.spacing = unit(0.1, 'cm'),     
          legend.margin = margin(t = -10))      


#D - Cell subtype DEG volcano plots
#Load in DE analysis
pre_unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Pre_Unciliated.csv")
unciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Unciliated.csv")
ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Ciliated.csv")
secretory_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Secretory.csv")
pre_ciliated_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Pre_Ciliated.csv")
proliferative_DE <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DE_analysis_Proliferative.csv")

#Define volcano plot colours
generate_keyvals <- function(df, lfc_col = "log2FoldChange", padj_col = "padj", 
                             lfc_thresh = 0.585, padj_thresh = 0.05,
                             col_up = "palegreen3", col_down = "#F8766D", col_ns = "#717171") {
  keyvals <- ifelse(df[[lfc_col]] < -lfc_thresh & df[[padj_col]] < padj_thresh, col_down,
                    ifelse(df[[lfc_col]] >  lfc_thresh & df[[padj_col]] < padj_thresh, col_up,
                           col_ns))
  keyvals[is.na(keyvals)] <- col_ns
  names(keyvals)[keyvals == col_up] <- "Up in Infertile"
  names(keyvals)[keyvals == col_down] <- "Down in Infertile"
  names(keyvals)[keyvals == col_ns] <- "Not Significant"
  return(keyvals)}
keyvals_pu <- generate_keyvals(pre_unciliated_DE)
keyvals_u  <- generate_keyvals(unciliated_DE)
keyvals_c  <- generate_keyvals(ciliated_DE)
keyvals_s  <- generate_keyvals(secretory_DE)
keyvals_pc <- generate_keyvals(pre_ciliated_DE)
keyvals_p  <- generate_keyvals(proliferative_DE)

#Volcano plots per cell subtype
pu <- EnhancedVolcano(pre_unciliated_DE, lab = pre_unciliated_DE$gene,
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585, pointSize = 4, labSize = 4, axisLabSize = 14,  
                      title = "Pre-Unciliated DEGs", titleLabSize = 18,
                      subtitle = NULL,  colCustom = keyvals_pu,
                      legendPosition = 'none', caption = NULL,
                      drawConnectors = TRUE, arrowheads = FALSE,
                      max.overlaps = 10) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank())
u <- EnhancedVolcano(unciliated_DE, lab = unciliated_DE$gene,
                     x = 'log2FoldChange', y = 'padj',
                     pCutoff = 0.05, FCcutoff = 0.585, pointSize = 4, labSize = 4, axisLabSize = 14,  
                     title = "Unciliated DEGs", titleLabSize = 18,
                     subtitle = NULL,  colCustom = keyvals_u,
                     legendPosition = 'none', caption = NULL,
                     drawConnectors = TRUE, arrowheads = FALSE,
                     max.overlaps = 15) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank())
c <- EnhancedVolcano(ciliated_DE, lab = ciliated_DE$gene,
                     x = 'log2FoldChange', y = 'padj', 
                     pCutoff = 0.05, FCcutoff = 0.585, pointSize = 4, labSize = 4, axisLabSize = 14,  
                     title = "Ciliated DEGs", titleLabSize = 18,
                     subtitle = NULL, colCustom = keyvals_c,
                     legendPosition = 'none', caption = NULL,
                     drawConnectors = TRUE, arrowheads = FALSE,
                     max.overlaps = 15) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank())
s <- EnhancedVolcano(secretory_DE, lab = secretory_DE$gene,
                     x = 'log2FoldChange', y = 'padj',
                     pCutoff = 0.05, FCcutoff = 0.585, pointSize = 4, labSize = 4, axisLabSize = 14,  
                     title = "Secretory DEGs", titleLabSize = 18,
                     subtitle = NULL, colCustom = keyvals_s,
                     legendPosition = 'none', caption = NULL,
                     drawConnectors = TRUE, arrowheads = FALSE,
                     max.overlaps = 15) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank())
pc <- EnhancedVolcano(pre_ciliated_DE, lab = pre_ciliated_DE$gene,
                      x = 'log2FoldChange', y = 'padj',
                      pCutoff = 0.05, FCcutoff = 0.585, pointSize = 4, labSize = 4, axisLabSize = 14, 
                      title = "Pre-Ciliated DEGs", titleLabSize = 18,
                      subtitle = NULL, colCustom = keyvals_pc,
                      legendPosition = 'none', caption = NULL,
                      drawConnectors = TRUE, arrowheads = FALSE,
                      max.overlaps = 15) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank())
p <- EnhancedVolcano(proliferative_DE, lab = proliferative_DE$gene,
                     x = 'log2FoldChange', y = 'padj',
                     pCutoff = 0.05, FCcutoff = 0.585, pointSize = 4, labSize = 4, axisLabSize = 14,  
                     title = "Proliferative DEGs", titleLabSize = 18,
                     subtitle = NULL, colCustom = keyvals_p,
                     legendPosition = 'none', caption = NULL,
                     drawConnectors = TRUE, arrowheads = FALSE,
                     max.overlaps = 15) +
                theme(plot.title = element_text(hjust = 0.5),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank())
pu + u + c + s + pc + p


#E - Top 10 DEGs per cluster
#Define gene lists per cell subtype
gene_lists <- list(
  "Pre-Unciliated" = c("AP000462.1","AL365273.2","CC2D2B","RERE-AS1","AP000462.3",
                       "AC092865.3","AL136309.4","AP003096.1","EIF1AXP2","MIR3658"),
  "Unciliated" = c("AL365273.2","AP000462.1","AL365273.1","RERE-AS1","AP000462.3",
                   "AC010542.5","AC005082.2","EIF1AXP2","AP003096.1","MIR3658"),
  "Ciliated" = c("AP000462.1","RF01871","AP000462.3","AL365273.2","AL356475.1",
                 "AC124784.1","AC010542.5","TRPM2","AP003096.1","MIR3658"),
  "Secretory" = c("AL109955.1","AL365273.2","AP000462.1","AP002748.4","IL1RL1",
                  "TRPM2-AS","EIF1AXP2","PRLR","AP003096.1","MIR3658"),
  "Pre-Ciliated" = c("AP000462.1","AL365273.2","AP002748.4","CC2D2B","RERE-AS1",
                     "AL133500.1","AP002954.1","AP003096.1","AC012379.1","MIR3658"),
  "Proliferative" = c("AL365273.2","AL109955.1","AP002748.4","AP000462.1","CC2D2B",
                      "CHST11","TAGLN","AC010542.5","AP001992.1","AP003096.1"))
#Store heatmaps
heatmap_grobs <- list()

#Loop through cell subtypes and generate heatmaps
for (cell_subtype in names(gene_lists)) {
  gene_list <- gene_lists[[cell_subtype]]
  
  #Subset cells of this subtype
  subtype_cells <- subset(all_samples_integrated, subset = cell_type == cell_subtype)
  
  #Filter for genes present in dataset
  gene_list_filtered <- gene_list[gene_list %in% rownames(subtype_cells)]
  if (length(gene_list_filtered) == 0) {
    message(paste("No genes found for", cell_subtype, "in dataset."))
    next}
  
  #Average expression
  avg_matrix <- AggregateExpression(subtype_cells, features = gene_list_filtered, 
                                    assays = "RNA", group.by = "fertility", return.seurat = FALSE)$RNA
  if (!"Fertile" %in% colnames(avg_matrix)) {
    message(paste("Skipping", cell_subtype, "- 'Fertile' group not found."))
    next}
  
  #log1p transform
  avg_matrix_log <- log1p(avg_matrix)
  
  #Create heatmap matrix
  heatmap_matrix <- as.matrix(avg_matrix_log[gene_list_filtered, , drop = FALSE])
  
  #Sort by fertile column
  gene_order <- order(heatmap_matrix[, "Fertile"], decreasing = TRUE)
  heatmap_matrix_ordered <- heatmap_matrix[gene_order, , drop = FALSE]
  
  #Save heatmap as grob
  p <- pheatmap(heatmap_matrix_ordered, scale = "none",
                cluster_rows = FALSE, cluster_columns = FALSE, 
                treeheight_col = 0, main = cell_subtype,
                fontsize_row = 8, fontsize_col = 10, angle_col = 0,
                silent = TRUE)
  
  #Modify title font size
  grob <- p[[4]]
  title_index <- which(sapply(grob$children, function(x) x$name) == "title")
  if (length(title_index) == 1) {
    grob$children[[title_index]] <- editGrob(
      grob$children[[title_index]],
      gp = gpar(fontsize = 15, fontface = "bold"))}
  heatmap_grobs[[cell_subtype]] <- grob}

#Plot all heatmaps together in a grid 
grid.arrange(grobs = heatmap_grobs, ncol = 3)



### ====================================================================================================================
###Supplementary Figure 5 - DEGs 2
#A - DEGs common to all cell subtypes
# Combine DE results (you already have this)
DE_results <- bind_rows(pre_unciliated_DE  %>% mutate(source = "Pre-Unciliated"),
                        unciliated_DE      %>% mutate(source = "Unciliated"),
                        ciliated_DE        %>% mutate(source = "Ciliated"),
                        secretory_DE       %>% mutate(source = "Secretory"),
                        pre_ciliated_DE    %>% mutate(source = "Pre-Ciliated"),
                        proliferative_DE   %>% mutate(source = "Proliferative"))

common_genes <- c("ABHD15-AS1","AC021733.3","AC026202.2","AC079447.1","AL109955.1","AL365273.2",
                  "AP000462.1","AP002748.4","AP003096.1","CC2D2B","KRT8","RPS10-NUDT3","SLC19A1","UPK1B")

#Filter for genes of interest and select relevant columns
filtered_DE <- DE_results %>%
  filter(gene %in% common_genes) %>%
  select(gene, source, log2FoldChange, padj)

#Pivot log2FoldChange wide for heatmap expression data
expression_data <- filtered_DE %>%
  select(gene, source, log2FoldChange) %>%
  pivot_wider(names_from = gene, values_from = log2FoldChange) %>%
  column_to_rownames("source")

#Reorder subtypes (rows)
sample_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
expression_data <- expression_data[sample_order, ]

#Order genes by increasing mean log2FC
gene_means <- colMeans(expression_data, na.rm = TRUE)
ordered_genes <- names(sort(gene_means))
expression_data <- expression_data[, ordered_genes]

#Pivot padj wide in same format for significance matrix
sig_data <- filtered_DE %>%
  select(gene, source, padj) %>%
  pivot_wider(names_from = gene, values_from = padj) %>%
  column_to_rownames("source")

#Reorder rows and columns in sig_data to match expression_data exactly
sig_data <- sig_data[sample_order, ordered_genes]

#Map padj to significance symbols
sig_mat <- matrix("", nrow = nrow(sig_data), ncol = ncol(sig_data),
                  dimnames = dimnames(sig_data))

sig_mat[sig_data < 0.0001] <- "****"
sig_mat[sig_data >= 0.0001 & sig_data < 0.001] <- "***"
sig_mat[sig_data >= 0.001 & sig_data < 0.01] <- "**"
sig_mat[sig_data >= 0.01 & sig_data < 0.05] <- "*"

#Plot heatmap with significance symbols
pheatmap(expression_data,
         cluster_rows = FALSE, cluster_cols = FALSE,
         angle_col = 45, cellwidth = 50,
         fontsize_row = 10, fontsize_col = 10,
         main = "DEGs in All Cell Subtypes",
         legend_breaks = seq(-8, 8, by = 2),
         display_numbers = sig_mat, number_color = "black", number_format = "%s")


#B - Global DEG GSEA analysis
#Define path to GSEA CSVs
gsea_files <- list.files(
  path = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/GSEA/",
  pattern = "^fgsea_.*\\.csv$",
  full.names = TRUE) %>%
  .[grepl("_bulk_", .)]

#Extract Subtype and Category from filename
parse_gsea_filename <- function(filename) {
  parts <- str_match(basename(filename), "fgsea_(.*?)_(.*?)\\.csv$")
  return(list(subtype = parts[,2], category = parts[,3]))}

#Load and combine GSEA files
fgsea_list <- lapply(gsea_files, function(file) {
  info <- parse_gsea_filename(file)
  df <- read_csv(file)
  if (nrow(df) == 0) return(NULL)
  df$Subtype <- info$subtype
  df$Category <- info$category
  return(df)})
fgsea_all <- bind_rows(fgsea_list)

#Filter for significant pathways
fgsea_all <- fgsea_all %>% filter(pval < 0.05)

#Keep top pathways per subtype and category
top_fgsea <- fgsea_all %>%
  group_by(Subtype, Category) %>%
  slice_min(order_by = pval, n = 50) %>%
  ungroup()

#Clean pathway names
top_fgsea$pathway <- top_fgsea$pathway %>%
  sub("^GOBP_|^HALLMARK_", "", .)

#Add direction and create ordered factor for plotting
top_fgsea <- top_fgsea %>%
  mutate(Direction = ifelse(NES >= 0, "Upregulated", "Downregulated"),
         Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) %>%
  group_by(Subtype, Category, Direction) %>%
  arrange(desc(NES)) %>%
  mutate(pathway_ordered = factor(pathway, levels = unique(pathway))) %>%
  ungroup()

#Flip the order of the y-axis so highest NES is at the top
top_fgsea$pathway_ordered <- fct_rev(top_fgsea$pathway_ordered)

#Rename factor levels for your single panel
top_fgsea$Subtype  <- factor(top_fgsea$Subtype, levels = unique(top_fgsea$Subtype), labels = "Infertile")

ggplot(top_fgsea, aes(x = NES, y = pathway_ordered, color = pval, size = size)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", name = "P-Value",
                       labels = scales::label_number(accuracy = 0.001),
                       guide = guide_colorbar(barwidth = 0.6, barheight = 12)) +  # switched width & height for vertical legend
  scale_size_continuous(name = "Genes",
                        labels = scales::number_format(accuracy = 1)) +
  scale_y_discrete(position = "right") +
  facet_grid2(rows = vars(Category), cols = vars(Subtype),
              scales = "free_y", space = "free", switch = "y",
              strip = strip_themed(background_x = elem_list_rect(
                fill = c("#F8766D"), color = "black"))) +
  labs(title = "Top Globally Enriched Infertile Pathways - GSEA",
       x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_bw(base_size = 12) +
  theme(legend.position = "left",
        legend.direction = "vertical",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "grey80", size = 0.4),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"))


#C - Apoptosis bubble plot
netVisual_bubble(cellchat_short, sources.use = c(1:6), targets.use = c(1:6),
                 signaling = c("IGFBP","TRAIL"), comparison = c(1, 2),
                 angle.x = 45, font.size = 8,
                 color.text = colours.fertility,
                 title.name = "Apoptosis Signalling Pathways L-R Interactions Divided by Fertility", 
                 font.size.title = 20) +
          theme(axis.text.x = element_text(size = 10),  
                axis.text.y = element_text(size = 10),
                legend.text = element_text(size = 10),  
                legend.title = element_text(size = 10), 
                legend.key.size = unit(0.5, "cm"),
                legend.position = "right",
                plot.title = element_text(size = 14, hjust = 0.5))


#D - IGFBP signalling LR Violin plot
plots <- plotGeneExpression(cellchat, signaling = "IGFBP", split.by = "datasets",
                            color.use = colours.fertility, type = "violin")
plots_rotated <- lapply(seq_along(plots), function(i) {
  if (i == length(plots)) {
    # Bottom plot: show rotated x-axis labels
    plots[[i]] + theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
                       axis.title.x = element_text()) } else {
                         # Other plots: hide x-axis text and ticks
                         plots[[i]] + theme(axis.text.x = element_blank(),
                                            axis.ticks.x = element_blank(),
                                            axis.title.x = element_blank())}})
wrap_plots(plots_rotated, ncol = 1)



### ====================================================================================================================
###Supplementary Figure 6 - DETs
#A - Coding isoforms
ggplot(category_summary %>% filter(!is.na(coding)),  # Remove NA values
       aes(x = fertility, y = num_isoforms, fill = coding)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = c("coding" = "#00BFC4", "non_coding" = "#A1E6E8"),
                    labels = c("coding" = "Coding", "non_coding" = "Non-Coding")) +
  labs(title = "Proportion of Coding Isoforms", 
       y = "Percentage of Isoforms", fill = "Coding Status") +
  scale_y_continuous(labels = scales::percent) + 
  scale_x_discrete(labels = c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_flip()


#B - Number of isoforms per gene
#Aggregate isoform expression by cell_type and fertility
counts <- AggregateExpression(all_samples_integrated, 
                              assays = "iso", return.seurat = FALSE,
                              group.by = "cell_type")

as.data.frame(counts) -> df
row.names(df) -> df$gene

#Split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "-",  extra = "merge") 

#Count the number of isoforms per gene
isoform_count_per_gene <- pseudobulk_data %>%
  group_by(gene_id) %>%
  summarise(n_isoforms = n_distinct(transcript_id))

#Count isforms per category 
isoform_count_per_gene <- isoform_count_per_gene %>%
  mutate(isoform_category = case_when(
      n_isoforms == 1 ~ "1",
      n_isoforms >= 2 & n_isoforms <= 3 ~ "2-3",
      n_isoforms >= 4 & n_isoforms <= 5 ~ "4-5",
      n_isoforms >= 6 & n_isoforms <= 9 ~ "6-9",
      n_isoforms >= 10 ~ "10+")) %>%
  dplyr::ungroup() %>%
  as_tibble()

#Calculate the percentage of genes in each bin
isoform_count_summary <- isoform_count_per_gene %>%
  dplyr::count(isoform_category) %>%
  mutate(percent = (n / sum(n)) * 100)

#Ensure consistent order for x-axis
isoform_count_summary$isoform_category <- factor(
  isoform_count_summary$isoform_category,
  levels = c("1", "2-3", "4-5", "6-9", "10+"))

ggplot(isoform_count_summary, aes(x = isoform_category, y = percent)) +
  geom_bar(stat = "identity", fill = "#849AFF", color = "black") +  
  labs(title = "Number of Isoforms per Gene", x = "Isoforms per Gene", y = "Genes Percentage") +
  theme_minimal() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))


#C - Transcript biotypes of all isoforms
#Aggregate the expression data by fertility
counts <- AggregateExpression(all_samples_integrated, assays = "iso", 
                              return.seurat = FALSE, group.by="fertility")

as.data.frame(counts) -> df
row.names(df) -> df$gene

#Split transcript ids into gene and transcript id
pseudobulk_data <- df %>% separate(gene, into = c("transcript_id", "gene_id"), sep = "-",  extra = "merge") 
pseudobulk_data$transcript_id <- sub("\\..*", "", pseudobulk_data$transcript_id)

#Select fertility columns and gather into long format
long_data <- pseudobulk_data %>%
  pivot_longer(cols = starts_with("iso."), names_to = "fertility", values_to = "expression")

#Filter for non-zero expression
filtered_data <- long_data %>%
  filter(expression > 0) %>%
  dplyr::select(fertility, transcript_id, gene_id) %>%
  distinct()

#Retrieve gene biotype from Ensembl
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://jan2020.archive.ensembl.org")
biotype_info_transcript <- getBM(attributes = c("ensembl_transcript_id","transcript_biotype","transcript_length"),
                                 filters = 'ensembl_transcript_id',
                                 values = filtered_data$transcript_id, mart = mart)

#Merge BioMart data into filtered expression data
merged_biomart_data <- merge(filtered_data, biotype_info_transcript, 
                             by.x = "transcript_id", by.y = "ensembl_transcript_id", all.x = TRUE)


#Reassign transcript_biotype to "BambuTx" if transcript_id starts with "BambuTx"
merged_biomart_data$transcript_id <- trimws(as.character(merged_biomart_data$transcript_id))
merged_biomart_data <- merged_biomart_data %>%
  mutate(transcript_biotype = ifelse(grepl("^BambuTx", transcript_id), "BambuTx", transcript_biotype))

#Prepare the summary data
category_summary <- merged_biomart_data %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = n(), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Group transcript_biotypes with proportion < 0.5% into "Other"
category_summary <- category_summary %>%
  mutate(transcript_biotype = ifelse(proportion < 0.005, "Other", as.character(transcript_biotype))) %>%
  group_by(fertility, transcript_biotype) %>%
  summarise(num_isoforms = sum(num_isoforms), .groups = "drop") %>%
  group_by(fertility) %>%
  mutate(proportion = num_isoforms / sum(num_isoforms)) %>%
  ungroup()

#Ensure consistent factor order for fertility and transcript_biotype
category_summary <- category_summary %>%
  mutate(fertility = factor(fertility, levels = unique(fertility)),
         transcript_biotype = factor(transcript_biotype, levels = c(
           "protein_coding", "lncRNA", "retained_intron", "processed_transcript",
           "nonsense_mediated_decay", "BambuTx", "processed_pseudogene", "TEC", "Other")))

#Create pie chart of transcript subtypes
ggplot(category_summary, aes(x = "", y = proportion, fill = transcript_biotype)) +
  geom_bar(stat = "identity", width = 1, colour="black", size=0.2) +  
  coord_polar(theta = "y", start = pi / 6) +
  facet_wrap(~fertility, labeller = as_labeller(c("iso.Fertile" = "Fertile", "iso.Infertile" = "Infertile"))) +
  scale_fill_discrete(
    labels = c("protein_coding" = "Protein Coding", "lncRNA" = "lncRNA", 
               "retained_intron" = "Retained Intron", "processed_transcript" = "Processed Transcript",
               "nonsense_mediated_decay" = "NMD","BambuTx" = "BambuTx",
               "processed_pseudogene" = "Processed Pseudogene", "TEC" = "TEC", "Other" = "Other")) +
  labs(title = "Transcript Biotypes by Fertility", fill = "Transcript Biotypes") +
  theme_void() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.margin = margin(10, 100, 10, 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold", margin = margin(b = 0)),
        panel.grid = element_blank())


#C - Cell subtype DEI volcano plots
#Load in DETs
pre_unciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Unciliated.csv")
unciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Unciliated.csv")
ciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Ciliated.csv")
secretory_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Secretory.csv")
pre_ciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Ciliated.csv")
proliferative_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Proliferative.csv")

#Define volcano plot colours
generate_keyvals <- function(df, lfc_col = "log2FoldChange", padj_col = "padj", 
                             lfc_thresh = 0.585, padj_thresh = 0.05,
                             col_up = "palegreen3", col_down = "#F8766D", col_ns = "#717171") {
  keyvals <- ifelse(df[[lfc_col]] < -lfc_thresh & df[[padj_col]] < padj_thresh, col_down,
                    ifelse(df[[lfc_col]] >  lfc_thresh & df[[padj_col]] < padj_thresh, col_up,
                           col_ns))
  keyvals[is.na(keyvals)] <- col_ns
  names(keyvals)[keyvals == col_up] <- "Up in Infertile"
  names(keyvals)[keyvals == col_down] <- "Down in Infertile"
  names(keyvals)[keyvals == col_ns] <- "Not Significant"
  return(keyvals)}
keyvals_pu <- generate_keyvals(pre_unciliated_DE_iso)
keyvals_u  <- generate_keyvals(unciliated_DE_iso)
keyvals_c  <- generate_keyvals(ciliated_DE_iso)
keyvals_s  <- generate_keyvals(secretory_DE_iso)
keyvals_pc <- generate_keyvals(pre_ciliated_DE_iso)
keyvals_p  <- generate_keyvals(proliferative_DE_iso)

#Volcano plots per cell subtype
pu_iso <- EnhancedVolcano(pre_unciliated_DE_iso, lab = pre_unciliated_DE_iso$transcript,
                          x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 0.585,
                          pointSize = 4, labSize = 4, axisLabSize = 14, 
                          title = "Pre-Unciliated DEIs", titleLabSize = 18,
                          subtitle = NULL,  colCustom = keyvals_pu,
                          legendPosition = 'none', caption = NULL,
                          drawConnectors = TRUE, arrowheads = FALSE,
                          max.overlaps = 6) +
                    theme(plot.title = element_text(hjust = 0.5),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
u_iso <- EnhancedVolcano(unciliated_DE_iso, lab = unciliated_DE_iso$transcript,
                         x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 0.585,
                         pointSize = 4, labSize = 4, axisLabSize = 14,   
                         title = "Unciliated DEIs", titleLabSize = 18,
                         subtitle = NULL,  colCustom = keyvals_u,
                         legendPosition = 'none', caption = NULL,
                         drawConnectors = TRUE, arrowheads = FALSE,
                         max.overlaps = 6) +
                   theme(plot.title = element_text(hjust = 0.5),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
c_iso <- EnhancedVolcano(ciliated_DE_iso, lab = ciliated_DE_iso$transcript,
                         x = 'log2FoldChange', y = 'padj',
                         pCutoff = 0.05, FCcutoff = 0.585,
                         pointSize = 4, labSize = 4, axisLabSize = 14,  
                         title = "Ciliated DEIs", titleLabSize = 18,
                         subtitle = NULL, colCustom = keyvals_c,
                         legendPosition = 'none', caption = NULL,
                         drawConnectors = TRUE, arrowheads = FALSE,
                         max.overlaps = 6) +
                  theme(plot.title = element_text(hjust = 0.5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
s_iso <- EnhancedVolcano(secretory_DE_iso, lab = secretory_DE_iso$transcript,
                         x = 'log2FoldChange', y = 'padj',
                         pCutoff = 0.05, FCcutoff = 0.585,
                         pointSize = 4, labSize = 4, axisLabSize = 14,  
                         title = "Secretory DEIs", titleLabSize = 18,
                         subtitle = NULL, colCustom = keyvals_s,
                         legendPosition = 'none', caption = NULL,
                         drawConnectors = TRUE, arrowheads = FALSE,
                         max.overlaps = 6) +
                  theme(plot.title = element_text(hjust = 0.5),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
pc_iso <- EnhancedVolcano(pre_ciliated_DE_iso, lab = pre_ciliated_DE_iso$transcript,
                          x = 'log2FoldChange', y = 'padj',
                          pCutoff = 0.05, FCcutoff = 0.585,
                          pointSize = 4, labSize = 4, axisLabSize = 14, 
                          title = "Pre-Ciliated DEIs", titleLabSize = 18,
                          subtitle = NULL, colCustom = keyvals_pc,
                          legendPosition = 'none', caption = NULL,
                          drawConnectors = TRUE, arrowheads = FALSE,
                          max.overlaps = 6) +
                    theme(plot.title = element_text(hjust = 0.5),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
p_iso <- EnhancedVolcano(proliferative_DE_iso, lab = proliferative_DE_iso$transcript,
                         x = 'log2FoldChange', y = 'padj',
                         pCutoff = 0.05, FCcutoff = 0.585,
                         pointSize = 4, labSize = 4, axisLabSize = 14, 
                         title = "Proliferative DEIs", titleLabSize = 18,
                         subtitle = NULL, colCustom = keyvals_p,
                         legendPosition = 'none', caption = NULL,
                         drawConnectors = TRUE, arrowheads = FALSE,
                         max.overlaps = 6) +
                    theme(plot.title = element_text(hjust = 0.5),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
pu_iso + u_iso + c_iso + s_iso + pc_iso + p_iso


#D - DETs common to all cell subtypes
#Load in DETs
bulk_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/bulk_DE_analysis.csv")
pre_unciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Unciliated.csv")
unciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Unciliated.csv")
ciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Ciliated.csv")
secretory_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Secretory.csv")
pre_ciliated_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Pre_Ciliated.csv")
proliferative_DE_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DE_analysis_Proliferative.csv")

#Combine DE results
DE_results_iso <- bind_rows(pre_unciliated_DE_iso  %>% mutate(source = "Pre-Unciliated"),
                            unciliated_DE_iso      %>% mutate(source = "Unciliated"),
                            ciliated_DE_iso        %>% mutate(source = "Ciliated"),
                            secretory_DE_iso       %>% mutate(source = "Secretory"),
                            pre_ciliated_DE_iso    %>% mutate(source = "Pre-Ciliated"),
                            proliferative_DE_iso   %>% mutate(source = "Proliferative"))
#DEGs in all clusters
common_iso <- c("ENST00000383329.7-HLA-C","BambuTx1337-PSORS1C3","ENST00000264234.8-UPK1B")

#Filter for genes of interest
filtered_DE_iso <- DE_results_iso %>%
  filter(transcript %in% common_iso) %>%
  dplyr::select(transcript, source, log2FoldChange)

#Reshape to wide format
expression_data <- filtered_DE_iso %>%
  pivot_wider(names_from = transcript, values_from = log2FoldChange) %>%
  column_to_rownames("source")

#Reorder rows and columns
sample_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", 
                  "Secretory", "Pre-Ciliated", "Proliferative")
expression_data <- expression_data[sample_order,]
expression_data_t <- t(expression_data)

#Calculate row distance and clustering manually
row_dist <- dist(expression_data_t)      
row_clust <- hclust(row_dist)           

#Reorder rows by clustering order
expression_data_t_ordered <- expression_data_t[row_clust$order, ]

#Fine-grained continuous color breaks for smooth gradient
breaks_continuous <- seq(-8, 8, by = 0.1)
legend_ticks <- seq(-8, 8, by = 2)

#Smooth red-yellow-blue palette (similar to pheatmap default)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaks_continuous) - 1)

pheatmap(expression_data_t_ordered,
         cluster_rows = FALSE, cluster_cols = FALSE,
         angle_col = 0, fontsize_row = 10,  fontsize_col = 10,
         main = "DEIs in All Cell Subtypes",
         breaks = breaks_continuous, legend_breaks = legend_ticks, color = my_colors)



###====================================================================================================================
###Supplementary Figure 7 - DET 2
#A - DEG-DEI overlap venn diagrams
#Function to extract gene symbol from transcript ID (everything after first hyphen)
extract_gene_from_transcript <- function(transcript) {
  sub("^[^-]*-", "", transcript)}

#Cell subtypes and nice names
subtypes <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")
display_names <- gsub("-", " ", subtypes)

#File paths for DEG and DEI csvs
deg_paths <- paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DEGs/DEGs_", gsub("-", "_", subtypes), ".csv")
det_paths <- paste0("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_", gsub("-", "_", subtypes), ".csv")

#Colors for venns: one for DEGs, one for DEIs
venn_colors_DEGs <- c("#F8766D", "#ABA300", "#0CB702", "#00BFC4", "#849AFF", "#FF61CC")
venn_colors_DEIs <- c("#FBBEB5","#E1DC85","#A8E7A3","#A1E6E8","#CBD2FF","#FFB3E9")

#Lists to hold all genes and transcripts for overall venn
all_deg_genes <- list()
all_dei_transcripts <- list()
all_dei_genes <- list()

#Store venn grobs for layout
venn_plots <- list()

for (i in seq_along(subtypes)) {
  #Load DEG gene data
  deg_df <- read.csv(deg_paths[i])
  deg_genes <- unique(deg_df$gene)
  
  #Load DEI transcript data
  det_df <- read.csv(det_paths[i])
  # Extract gene from transcript ID
  det_df$gene <- extract_gene_from_transcript(det_df$transcript)
  
  #Filter significant DEIs if needed
  det_df_filtered <- det_df[det_df$padj < 0.05 & !is.na(det_df$padj), ]
  dei_transcripts <- unique(det_df_filtered$transcript)
  dei_genes <- det_df_filtered$gene
  
  #Calculate counts for venn
  num_deg <- length(deg_genes)
  num_dei <- length(dei_transcripts)
  overlap <- sum(dei_genes %in% deg_genes)
  
  #Save for overall venn
  all_deg_genes[[subtypes[i]]] <- deg_genes
  all_dei_transcripts[[subtypes[i]]] <- dei_transcripts
  all_dei_genes[[subtypes[i]]] <- dei_genes
  
  #Create Venn diagram grob
  venn_grob <- draw.pairwise.venn(
    area1 = num_deg,  area2 = num_dei,
    cross.area = overlap,
    category = c("DEGs", "DEIs"),
    fill = c(venn_colors_DEGs[i], venn_colors_DEIs[i]),
    alpha = 0.5,  cex = 1.8,
    cat.cex = 1.3, cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05),
    fontfamily = "Arial", cat.fontfamily = "Arial",
    scaled = FALSE, margin = 0.1)
  
  #Combine Venn and label
  label_grob <- textGrob(display_names[i], y = unit(0.95, "npc"),
                         gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Arial"))
  
  venn_plots[[i]] <- grobTree(gList(venn_grob, label_grob))}

grid.newpage()
grid.arrange(grobs = venn_plots, ncol = 6,
             top = textGrob("", gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "Arial")))


#B - DEI biotype classification
#Load DETs
pre_unciliated_df <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Pre_Unciliated.csv")
unciliated_df     <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Unciliated.csv")
ciliated_df       <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Ciliated.csv")
secretory_df      <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Secretory.csv")
pre_ciliated_df   <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Pre_Ciliated.csv")
proliferative_df  <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/DETs/DETs_Proliferative.csv")

#Merge DET data frames
det_list <- list('Pre-Unciliated' = pre_unciliated_df, 'Unciliated' = unciliated_df,
                 'Ciliated' = ciliated_df,             'Secretory' = secretory_df,
                 'Pre-Ciliated' = pre_ciliated_df,     'Proliferative' = proliferative_df)

#Add cell_type column and merge into one table
all_dets <- bind_rows(lapply(names(det_list), function(cell) {det_list[[cell]] %>%
    mutate(cell_type = cell)}))

#Split and rename columns
all_dets <- all_dets %>%
  mutate(transcript_id = sub("-.*", "", transcript), gene_id = sub("^[^\\-]+-", "", transcript)) %>%
  dplyr::select(-transcript)

#SQANTI
sqanti <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/remove_unknownstrand_classification.txt", sep ='\t') %>%
  rename(transcript_id = isoform)
all_dets_annot <- all_dets %>%
  left_join(sqanti %>% dplyr::select(transcript_id, structural_category), by = "transcript_id")

#BIOMART
all_dets_annot$transcript_id_no_version <- sub("\\.[0-9]+$", "", all_dets_annot$transcript_id)
biomart_annot <- getBM(attributes = c("ensembl_transcript_id", 'transcript_biotype'),
                       filters = 'ensembl_transcript_id',
                       values = all_dets_annot$transcript_id_no_version, mart = mart)

#Merge BioMart and define Bambu transcripts as novel
biomart_annot <- biomart_annot %>%
  rename(transcript_id_no_version = ensembl_transcript_id)
det_annot <- all_dets_annot %>%
  left_join(biomart_annot, by = "transcript_id_no_version") %>%
  mutate(transcript_biotype = ifelse(grepl("^Bambu", transcript_id), "BambuTx", transcript_biotype))

#Define your desired cell subtype order
cell_order <- c("Proliferative","Secretory","Ciliated","Pre-Ciliated","Unciliated","Pre-Unciliated")

#Plot BioMart biotypes as percentages
biotype_counts <- det_annot %>%
  dplyr::count(cell_type, transcript_biotype) %>%
  group_by(cell_type) %>%
  mutate(percent = n / sum(n) * 100) %>%
  mutate(cell_type = factor(cell_type, levels = cell_order)) %>%
  ungroup()
biotype_counts$transcript_biotype <- factor(biotype_counts$transcript_biotype,
                                            levels = c("protein_coding", "lncRNA", "processed_transcript", "retained_intron",
                                                       "nonsense_mediated_decay", "processed_pseudogene", "TEC", "BambuTx"),
                                            labels = c("Protein Coding", "lncRNA", "Processed Transcript", "Retained Intron",
                                                       "NMD","Processed Pseudogene", "TEC", "BambuTx"))

ggplot(biotype_counts, aes(x = cell_type, y = percent, fill = transcript_biotype)) +
  geom_col(position = "stack") +  coord_flip() +
  labs(title = "DEI Transcript Biotype by Cell Subtype", y = "Percentage of DEIs", fill = "Transcript Biotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#Biotype summary
biotype_summary <- biotype_counts %>%
  group_by(transcript_biotype) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  arrange(desc(total_n)) %>%
  mutate(percent = 100 * total_n / sum(total_n), cumulative_percent = cumsum(percent))

biotype_summary$transcript_biotype <- factor(biotype_summary$transcript_biotype,
                                            labels = c("Protein Coding", "lncRNA", "Processed Transcript", "Retained Intron",
                                                       "NMD","Processed Pseudogene", "TEC", "BambuTx"))

ggplot(biotype_summary, aes(x = "All Cells", y = percent, fill = transcript_biotype)) +
  geom_col(position = "stack") +
  coord_flip() +
  labs(title = "DEI Transcript Biotypes - All Cells", y = "Percentage of DEIs", fill = "Transcript Biotype") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 15, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#C - Bar plot of DTU genes and transcripts comparison
#Summarize DET vs non-DET counts per dataset/cell type
counts <- combined_data %>%
  group_by(cell_type, dataset, DET) %>%
  summarise(transcript_count = n_distinct(transcript_id), .groups = "drop") %>%
  complete(cell_type, dataset, DET, fill = list(transcript_count = 0))

#Gene and transcript count summary
summary_df <- combined_data %>%
  group_by(cell_type, dataset) %>%
  summarise(gene_count = n_distinct(gene_name),
            transcript_count = n_distinct(transcript_id), .groups = "drop") %>%
  pivot_longer(cols = c(gene_count, transcript_count), names_to = "type", values_to = "count") %>%
  mutate(type = recode(type, gene_count = "Genes", transcript_count = "Transcripts"),
         dataset_type = paste(dataset, type, sep = "_"))

summary_df$cell_type <- factor(summary_df$cell_type, levels = desired_order)

#Custom fill colors
custom_colors <- c("ISA_Genes"           = "#F8766D",
                   "ISA_Transcripts"     = "#ABA300",
                   "DTUrtle_Genes"       = "#0CB702",
                   "DTUrtle_Transcripts" = "#00BFC4",
                   "Isopod_Genes"        = "#849AFF",
                   "Isopod_Transcripts"  = "#FF61CC")

#D - Plot gene and transcript bar plot
ggplot(summary_df, aes(x = cell_type, y = count, fill = dataset_type)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = count), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "DTU by Cell Subtype - Analysis Comparison", y = "Count", fill = "Dataset") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 0),
        legend.position = "bottom")



###====================================================================================================================
###Supplementary Figure 8 - HLA-C
#A - Isoform division of HLA-C 
features <- rownames(all_samples_integrated@assays$iso@features)
gene <- "HLA-C"
expression_matrix <- GetAssayData(all_samples_integrated, assay = "iso", slot = "data")

#Filter features containing the gene name
matching_features <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"), rownames(expression_matrix), value = TRUE)

#Subset the expression matrix to include only the matching features
subset_expression <- expression_matrix[matching_features, , drop = FALSE]

#Calculate the total expression for each matching feature
total_expression <- Matrix::rowSums(subset_expression)

#Rank features by average expression
top_features <- names(sort(total_expression, decreasing = TRUE))

#Plot the top features in descending order of their average expression
plots <- FeaturePlot(all_samples_integrated, features = head(top_features, 12), reduction = "umap",
                     order = TRUE, pt.size = 1)

#Adjust title size for each plot
plots <- lapply(plots, function(plot) {
  plot + theme(plot.title = element_text(size = 14))})
CombinePlots(plots = plots, ncol = 4)    


#B - Check for patient bias - HLA-C isoform
VlnPlot(all_samples_integrated, features = "ENST00000383329.7-HLA-C",
        split.by="cell_type", group.by="sample",
        cols=cluster_colours) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ggtitle("ENST00000383329.9-HLA-C Isoform Expression by Sample")


                              
### ====================================================================================================================
###Supplementary Figure 9 - HLA-C Mapping
#A + B - Generated using IGV

                              
### ====================================================================================================================
###Supplementary Figure 10 - RBCK1
#A - Isoform division of RBCK1
features <- rownames(all_samples_integrated@assays$iso@features)
gene <- "RBCK1"
expression_matrix <- GetAssayData(all_samples_integrated, assay = "iso", slot = "data")

#Filter features containing the gene name
matching_features <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"), rownames(expression_matrix), value = TRUE)

#Subset the expression matrix to include only the matching features
subset_expression <- expression_matrix[matching_features, , drop = FALSE]

#Calculate the total expression for each matching feature
total_expression <- Matrix::rowSums(subset_expression)

#Rank features by average expression
top_features <- names(sort(total_expression, decreasing = TRUE))

#Plot the top features in descending order of their average expression
plots <- FeaturePlot(all_samples_integrated, features = head(top_features, 12), reduction = "umap",
                     order = TRUE, pt.size = 1)

#Adjust title size for each plot
plots <- lapply(plots, function(plot) {
  plot + theme(plot.title = element_text(size = 14))})
CombinePlots(plots = plots, ncol = 4)


#B - AKR1C1
genes <- c("AKR1C1")
sample_order <- c("Pre-Unciliated", "Unciliated", "Ciliated", "Secretory", "Pre-Ciliated", "Proliferative")

filtered_DE <- DE_results %>%
  filter(gene %in% genes)
expression_data <- filtered_DE %>%
  dplyr::select(gene, source, log2FoldChange) %>%
  spread(key = gene, value = log2FoldChange) %>%
  column_to_rownames(var = "source")
expression_data <- expression_data[sample_order, , drop = FALSE]

#Extract gene name if transcript has combined label
DE_results_iso <- DE_results_iso %>%
  mutate(gene = sub(".*-", "", transcript))
filtered_iso_DE <- DE_results_iso %>%
  filter(gene %in% genes) %>%
  mutate(isoform_label = transcript)

#Recalculate isoform order
isoform_order <- filtered_iso_DE %>%
  dplyr::select(isoform_label, gene) %>%
  distinct() %>%
  arrange(gene, isoform_label) %>%
  pull(isoform_label)

#Reshape for heatmap
expression_data_iso <- filtered_iso_DE %>%
  dplyr::select(source, isoform_label, log2FoldChange) %>%
  pivot_wider(names_from = isoform_label, values_from = log2FoldChange) %>%
  column_to_rownames("source")
expression_data_iso <- expression_data_iso[sample_order, isoform_order, drop = FALSE]
expression_data_iso <- expression_data_iso %>%
  select(-'ENST00000477661.1-AKR1C1')
combined_expression <- cbind(expression_data, expression_data_iso)

#Heatmap settings
breaksList <- seq(-2, 2, length.out = 100)
heat_colors <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(length(breaksList)))
pheatmap(combined_expression,
         cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 0,
         main = "Gene and Isoform-Level Differential Expression of AKR1C1",
         fontsize_row = 11, fontsize_col = 10,
         color = heat_colors, breaks = breaksList)


#C - NOTCH L-R CellChat
netVisual_bubble(cellchat_short,
                 sources.use = c(1:6), targets.use = c(1:6),
                 signaling = c("NOTCH"), comparison = c(1, 2),
                 angle.x = 45, font.size = 8,
                 color.text = colours.fertility,
                 title.name = "NOTCH Signalling Pathways L-R Interactions between Cell Subtypes Divided by Fertility", 
                 font.size.title = 20) +
          theme(axis.text.x = element_text(size = 11.5),  # Increase x-axis text size
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 12),  # Increase legend text size
                legend.title = element_text(size = 12, face = "bold"), # Increase legend title size
                legend.key.size = unit(1, "cm"),
                legend.position = "bottom",
                plot.title = element_text(margin = margin(b = 30), hjust = 0.5, face = "bold")) +
          guides(color = guide_colorbar(title = "Communication Prob", barwidth = 9, barheight = 0.7,
                                        title.vjust = 1.5),
                 size = guide_legend(title = "P-value", 
                                     override.aes = list(shape = 21, fill ="black", color ="black"),
                                     label.position = "right", title.vjust = 1.5))
