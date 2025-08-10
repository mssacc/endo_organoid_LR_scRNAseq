#Merge multisamples, run unintegrated analysis, run integrated analysis

#Set library
.libPaths("/home/mssacc/R_libs_4.4")

#Load required libraries
library(SeuratObject)
library(Seurat)
library(harmony)
library(devtools)
library(presto)
library(dplyr)
library(ggplot2)
library(ggraph)
library(clustree)
library(dittoSeq)

#Load gene and isoform objects
E170 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_gene_and_isoform_seurat.rds")
E191 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_gene_and_isoform_seurat.rds")
E226 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_gene_and_isoform_seurat.rds")
E231 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_gene_and_isoform_seurat.rds")
E333 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_gene_and_isoform_seurat.rds")
E435 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_gene_and_isoform_seurat.rds")

#Set seed
set.seed(123)

#Merge datasets
all_samples <- merge(E170, y=c(E191, E226, E231, E333, E435), 
                           add.cell.ids=c("E170", "E191", "E226", "E231", "E333","E435"), 
                           project="all_samples")

#Rename orig.ident
all_samples@meta.data$orig.ident[all_samples@meta.data$orig.ident == "E170_genes"] <- "E170"
all_samples@meta.data$orig.ident[all_samples@meta.data$orig.ident == "E191_genes"] <- "E191"
all_samples@meta.data$orig.ident[all_samples@meta.data$orig.ident == "E226_genes"] <- "E226"
all_samples@meta.data$orig.ident[all_samples@meta.data$orig.ident == "E231_genes"] <- "E231"
all_samples@meta.data$orig.ident[all_samples@meta.data$orig.ident == "E333_genes"] <- "E333"
all_samples@meta.data$orig.ident[all_samples@meta.data$orig.ident == "E435_genes"] <- "E435"

#Create a column for fertility
all_samples$fertility=all_samples$orig.ident
all_samples$fertility[all_samples$fertility %in% c("E170", "E226", "E231")] <- "Infertile"
all_samples$fertility[all_samples$fertility %in% c("E191", "E333", "E435")] <- "Fertile"
#Create a column for batch
all_samples$batch=all_samples$orig.ident
all_samples$batch[all_samples$batch %in% c("E191", "E231")] <- "Batch-1"
all_samples$batch[all_samples$batch %in% c("E333", "E170")] <- "Batch-2"
all_samples$batch[all_samples$batch %in% c("E435", "E226")] <- "Batch-3"
#Create a column for endo ID
all_samples$endo.ID=all_samples$orig.ident
all_samples$endo.ID[all_samples$endo.ID %in% "E170"] <- "E170"
all_samples$endo.ID[all_samples$endo.ID %in% "E191"] <- "E191"
all_samples$endo.ID[all_samples$endo.ID %in% "E226"] <- "E226"
all_samples$endo.ID[all_samples$endo.ID %in% "E231"] <- "E231"
all_samples$endo.ID[all_samples$endo.ID %in% "E333"] <- "E333"
all_samples$endo.ID[all_samples$endo.ID %in% "E435"] <- "E435"
#Create a column for fertility + endo ID
all_samples$fertility_ID=all_samples$orig.ident
all_samples$fertility_ID[all_samples$fertility_ID %in% "E170"] <- "Infertile_E170"
all_samples$fertility_ID[all_samples$fertility_ID %in% "E191"] <- "Fertile_E191"
all_samples$fertility_ID[all_samples$fertility_ID %in% "E226"] <- "Infertile_E226"
all_samples$fertility_ID[all_samples$fertility_ID %in% "E231"] <- "Infertile_E231"
all_samples$fertility_ID[all_samples$fertility_ID %in% "E333"] <- "Fertile_E333"
all_samples$fertility_ID[all_samples$fertility_ID %in% "E435"] <- "Fertile_E435"
#Create a column for sample number
all_samples$sample=all_samples$orig.ident
all_samples$sample[all_samples$sample %in% "E170"] <- "Infertile 2"
all_samples$sample[all_samples$sample %in% "E191"] <- "Fertile 1"
all_samples$sample[all_samples$sample %in% "E226"] <- "Infertile 3"
all_samples$sample[all_samples$sample %in% "E231"] <- "Infertile 1"
all_samples$sample[all_samples$sample %in% "E333"] <- "Fertile 2"
all_samples$sample[all_samples$sample %in% "E435"] <- "Fertile 3"

#Confirm combining all data worked
table(all_samples$orig.ident)
table(all_samples$fertility)
table(all_samples$batch)

#Perform analysis without integration
all_samples <- NormalizeData(all_samples)
all_samples <- FindVariableFeatures(all_samples)
all_samples <- ScaleData(all_samples, features= rownames(all_samples))
all_samples <- RunPCA(all_samples)
ElbowPlot(all_samples)
all_samples <- FindNeighbors(all_samples, dims=1:10, reduction="pca")
all_samples <- FindClusters(all_samples, resolution=0.6, cluster.name="unintegrated_cluster")
all_samples <- RunUMAP(all_samples, dims=1:10, reduction="pca", reduction.name="umap.unintegrated", seed.use=123)

#QC filtering together - look at data distribution
VlnPlot(all_samples, group.by="endo.ID", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(all_samples, group.by="endo.ID", features = c("percent.mt"))

FeaturePlot(all_samples, features = c("percent.mt","nFeature_RNA"))
DimPlot(all_samples, reduction = "umap.unintegrated", group.by = c("seurat_clusters", "orig.ident", "fertility"))
DimPlot(all_samples, reduction = "umap.unintegrated")

#Perform integration
all_samples_integrated <- IntegrateLayers(object=all_samples,
                                          method=HarmonyIntegration,
                                          orig.reduction="pca", new.reduction="integrated",
                                          verbose=FALSE)
all_samples_integrated <- FindNeighbors(all_samples_integrated, dims=1:10, reduction="integrated")
all_samples_integrated <- FindClusters(all_samples_integrated, resolution=0.25, cluster.name="integrated_cluster")
all_samples_integrated <- RunUMAP(all_samples_integrated, reduction="integrated", dims=1:10, seed.use=123)
pca <- RunPCA(all_samples_integrated)

#Visualize data
DimPlot(all_samples_integrated, reduction="umap")
DimPlot(all_samples_integrated, reduction="umap", group.by=c("seurat_clusters", "orig.ident", "fertility"))
DimPlot(all_samples_integrated, reduction="umap", split.by="fertility")
DimPlot(all_samples_integrated, reduction="umap", split.by="batch")
DimPlot(all_samples_integrated, reduction="umap", split.by="endo.ID")


#Convert to character to allow renaming
all_samples_integrated$cell_type <- as.character(all_samples_integrated$seurat_clusters)
#Rename clusters
all_samples_integrated$cell_type[all_samples_integrated$cell_type == "0"] <- "Pre-Unciliated"
all_samples_integrated$cell_type[all_samples_integrated$cell_type == "1"] <- "Unciliated"
all_samples_integrated$cell_type[all_samples_integrated$cell_type == "2"] <- "Ciliated"
all_samples_integrated$cell_type[all_samples_integrated$cell_type == "3"] <- "Secretory"
all_samples_integrated$cell_type[all_samples_integrated$cell_type == "4"] <- "Pre-Ciliated"
all_samples_integrated$cell_type[all_samples_integrated$cell_type == "5"] <- "Proliferative"
#Explicitly set factor levels for 'cell_type'
all_samples_integrated$cell_type <- factor(all_samples_integrated$cell_type,
                                        levels = c("Pre-Unciliated", "Unciliated", "Ciliated","Secretory", "Pre-Ciliated", "Proliferative"))
#Plot cell subtype percentages
plot <- dittoBarPlot(all_samples_integrated,
                "cell_type",
                group.by = "fertility",
                retain.factor.levels = TRUE,  # Retain the custom factor levels
                main = "Cell Subtype Percentages", sub = NULL, xlab = NULL,
                x.labels.rotate = FALSE,
                color.panel = c("#F8766D", "#ABA300", "#0CB702", "#00BFC4", "#849AFF", "#FF61CC"))
plot <- plot + theme(plot.title = element_text(size = 18, face = "bold"),
                     axis.text = element_text(size = 12),  # Increase tick mark label size
                     axis.title = element_text(size = 14),
                     legend.text = element_text(size = 12)) # Increase axis title size
plot

#Reset clusters idents
all_samples_integrated <- RenameIdents(all_samples_integrated,
                                       "0" = "Pre-Unciliated",
                                       "1" = "Unciliated",
                                       "2" = "Ciliated",
                                       "3" = "Secretory",
                                       "4" = "Pre-Ciliated",
                                       "5" = "Proliferative")

#Regenerate graphs with cluster names
plot1 <- DimPlot(all_samples_integrated, reduction="umap")+
  theme(legend.text = element_text(size = 18)) 
plot1
DimPlot(all_samples_integrated, reduction="umap", group.by="cell_type")
DimPlot(all_samples_integrated, reduction="umap", group.by="cell_type", split.by="fertility")

#Rejoin gene datasets after integration
all_samples_integrated <- JoinLayers(all_samples_integrated)
all_samples_integrated <- JoinLayers(all_samples_integrated, assay="iso")

#Save RDS
saveRDS(all_samples_integrated, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")
