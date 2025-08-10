#Rscript - Gene QC (doublets, empty droplets, ambient RNA, low cell quality)

#Load required libraries
library(DropletUtils)
library(SeuratObject)
library(Seurat)
library(BiocParallel)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(devtools)
library(DoubletFinder)
library(patchwork)
library(cowplot)

#Define some factors
lower = 100
fdr_threshold = 0.001
output_path = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/"

#Load gene data
E170_genes <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_gene_count_symbols.csv", header=T, row.names=1)
E191_genes <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_gene_count_symbols.csv", header=T, row.names=1)
E226_genes <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_gene_count_symbols.csv", header=T, row.names=1)
E231_genes <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_gene_count_symbols.csv", header=T, row.names=1)
E333_genes <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_gene_count_symbols.csv", header=T, row.names=1)
E435_genes <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_gene_count_symbols.csv", header=T, row.names=1)


###Empty droplet filtering
#Load empty droplet data
E170_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_empty_droplets.csv", header=T, row.names = 1)
E191_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_empty_droplets.csv", header=T, row.names = 1)
E226_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_empty_droplets.csv", header=T, row.names = 1)
E231_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_empty_droplets.csv", header=T, row.names = 1)
E333_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_empty_droplets.csv", header=T, row.names = 1)
E435_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_empty_droplets.csv", header=T, row.names = 1)

#Combine the gene and empty droplet dataframes by row names
#E170
  E170_combined_df <- merge(E170_genes, E170_emptydrops, by = "row.names", all = TRUE)
  rownames(E170_combined_df) <- E170_combined_df[, 1]
  E170_combined_df[, 1] <- NULL
  E170_combined_df[is.na(E170_combined_df)] <- 0
#E191
  E191_combined_df <- merge(E191_genes, E191_emptydrops, by = "row.names", all = TRUE)
  rownames(E191_combined_df) <- E191_combined_df[, 1]
  E191_combined_df[, 1] <- NULL
  E191_combined_df[is.na(E191_combined_df)] <- 0
#E226
  E226_combined_df <- merge(E226_genes, E226_emptydrops, by = "row.names", all = TRUE)
  rownames(E226_combined_df) <- E226_combined_df[, 1]
  E226_combined_df[, 1] <- NULL
  E226_combined_df[is.na(E226_combined_df)] <- 0
#E231
  E231_combined_df <- merge(E231_genes, E231_emptydrops, by = "row.names", all = TRUE)
  rownames(E231_combined_df) <- E231_combined_df[, 1]
  E231_combined_df[, 1] <- NULL
  E231_combined_df[is.na(E231_combined_df)] <- 0
#E333
  E333_combined_df <- merge(E333_genes, E333_emptydrops, by = "row.names", all = TRUE)
  rownames(E333_combined_df) <- E333_combined_df[, 1]
  E333_combined_df[, 1] <- NULL
  E333_combined_df[is.na(E333_combined_df)] <- 0
#E435
  E435_combined_df <- merge(E435_genes, E435_emptydrops, by = "row.names", all = TRUE)
  rownames(E435_combined_df) <- E435_combined_df[, 1]
  E435_combined_df[, 1] <- NULL
  E435_combined_df[is.na(E435_combined_df)] <- 0

#Make seurat objects
E170_gene_object <- CreateSeuratObject(counts=E170_genes, project="E170")
E191_gene_object <- CreateSeuratObject(counts=E191_genes, project="E191")
E226_gene_object <- CreateSeuratObject(counts=E226_genes, project="E226")
E231_gene_object <- CreateSeuratObject(counts=E231_genes, project="E231")
E333_gene_object <- CreateSeuratObject(counts=E333_genes, project="E333")
E435_gene_object <- CreateSeuratObject(counts=E435_genes, project="E435")

##Perform standard pre-processing for each individual sample
#E170 Processing
    E170_gene_object <- NormalizeData(E170_gene_object, normalization.method="LogNormalize", scale.factor=10000)
    E170_gene_object <- FindVariableFeatures(E170_gene_object, selection.method="vst", nfeatures=2000)
    E170_gene_object <- ScaleData(E170_gene_object, features=rownames(E170_gene_object))
    E170_gene_object <- RunPCA(E170_gene_object, features=VariableFeatures(object=E170_gene_object))
    ElbowPlot(E170_gene_object)
    E170_gene_object <- FindNeighbors(E170_gene_object, dims=1:10)
    E170_gene_object <- FindClusters(E170_gene_object, resolution=0.2)
    E170_gene_object <- RunUMAP(E170_gene_object, dims=1:10)
    DimPlot(E170_gene_object, reduction = "umap")
#E191 Processing
    E191_gene_object <- NormalizeData(E191_gene_object, normalization.method="LogNormalize", scale.factor=10000)
    E191_gene_object <- FindVariableFeatures(E191_gene_object, selection.method="vst", nfeatures=2000)
    E191_gene_object <- ScaleData(E191_gene_object, features=rownames(E191_gene_object))
    E191_gene_object <- RunPCA(E191_gene_object, features=VariableFeatures(object=E191_gene_object))
    ElbowPlot(E191_gene_object)
    E191_gene_object <- FindNeighbors(E191_gene_object, dims=1:10)
    E191_gene_object <- FindClusters(E191_gene_object, resolution=0.2)
    E191_gene_object <- RunUMAP(E191_gene_object, dims=1:10)
    DimPlot(E191_gene_object, reduction = "umap")    
#E226 Processing
    E226_gene_object <- NormalizeData(E226_gene_object, normalization.method="LogNormalize", scale.factor=10000)
    E226_gene_object <- FindVariableFeatures(E226_gene_object, selection.method="vst", nfeatures=2000)
    E226_gene_object <- ScaleData(E226_gene_object, features=rownames(E226_gene_object))
    E226_gene_object <- RunPCA(E226_gene_object, features=VariableFeatures(object=E226_gene_object))
    ElbowPlot(E226_gene_object)
    E226_gene_object <- FindNeighbors(E226_gene_object, dims=1:10)
    E226_gene_object <- FindClusters(E226_gene_object, resolution=0.2)
    E226_gene_object <- RunUMAP(E226_gene_object, dims=1:10)
    DimPlot(E226_gene_object, reduction = "umap")    
#E231 Processing
    E231_gene_object <- NormalizeData(E231_gene_object, normalization.method="LogNormalize", scale.factor=10000)
    E231_gene_object <- FindVariableFeatures(E231_gene_object, selection.method="vst", nfeatures=2000)
    E231_gene_object <- ScaleData(E231_gene_object, features= rownames(E231_gene_object))
    E231_gene_object <- RunPCA(E231_gene_object, features=VariableFeatures(object=E231_gene_object))
    ElbowPlot(E231_gene_object)
    E231_gene_object <- FindNeighbors(E231_gene_object, dims=1:10)
    E231_gene_object <- FindClusters(E231_gene_object, resolution=0.2)
    E231_gene_object <- RunUMAP(E231_gene_object, dims=1:10)
    DimPlot(E231_gene_object, reduction = "umap")    
#E333 Processing
    E333_gene_object <- NormalizeData(E333_gene_object, normalization.method="LogNormalize", scale.factor=10000)
    E333_gene_object <- FindVariableFeatures(E333_gene_object, selection.method="vst", nfeatures=2000)
    E333_gene_object <- ScaleData(E333_gene_object, features= rownames(E333_gene_object))
    E333_gene_object <- RunPCA(E333_gene_object, features=VariableFeatures(object=E333_gene_object))
    ElbowPlot(E333_gene_object)
    E333_gene_object <- FindNeighbors(E333_gene_object, dims=1:10)
    E333_gene_object <- FindClusters(E333_gene_object, resolution=0.2)
    E333_gene_object <- RunUMAP(E333_gene_object, dims=1:10)
    DimPlot(E333_gene_object, reduction = "umap")    
#E435 Processing
    E435_gene_object <- NormalizeData(E435_gene_object, normalization.method="LogNormalize", scale.factor=10000)
    E435_gene_object <- FindVariableFeatures(E435_gene_object, selection.method="vst", nfeatures=2000)
    E435_gene_object <- ScaleData(E435_gene_object, features= rownames(E435_gene_object))
    E435_gene_object <- RunPCA(E435_gene_object, features=VariableFeatures(object=E435_gene_object))
    ElbowPlot(E435_gene_object)
    E435_gene_object <- FindNeighbors(E435_gene_object, dims=1:10)
    E435_gene_object <- FindClusters(E435_gene_object, resolution=0.2)
    E435_gene_object <- RunUMAP(E435_gene_object, dims=1:10)
    DimPlot(E435_gene_object, reduction = "umap")

#E170 Empty Drops Analysis
  #Define function to make dgCMatrix from combined counts
    makedgcmatrix <- function(count.matrix) {
      E170_gene_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell")
      list(E170_gene_object[["RNA"]]$counts)}
    E170_combined_df[] <- lapply(E170_combined_df, function(x) as.numeric(as.character(x)))
    outs.ddcmatrix <- makedgcmatrix(E170_combined_df)[[1]]
    br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
    e.out <- emptyDrops(outs.ddcmatrix, lower=100, niters=10000, test.ambient=TRUE, BPPARAM=SerialParam())
    is.cell <- e.out$FDR < fdr_threshold
    #Create a dataframe with FDR of TRUE cells
      is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
      is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
      is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
    #Function for retrieving the Seurat cells and cluster in dataframe
      overlap_true_cell <- function(E170_gene_object) {
        seurat_cluster.df <- as.data.frame(E170_gene_object$seurat_clusters)
        seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
        seurat_cluster.df}
    #Obtain cluster dataframe from Seurat object
      overlap_CR <- overlap_true_cell(E170_gene_object)
    #Check overlaps between Seurat object and true cells
      summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
    #Function to add metadata to Seurat object
      True.cells <- function(e.out) {
        cells <- as.data.frame(e.out@rownames)
        fdr <- as.data.frame(e.out$FDR)
        T.F.cells <- cbind(cells, fdr)
        T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
        setnames(T.F.cells, c('FDR'))
        T.F.cells %>%
          mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))}
      cells_CR <- True.cells(e.out)
      E170_gene_object <- AddMetaData(E170_gene_object, metadata=cells_CR, col.name='is.cell')
  #Plot Empty drops on Gene UMAP
    #Create a ggplot object
    rankplot <- ggplot(br.out, aes(x=rank, y=total)) +
      geom_point() + scale_x_log10() + scale_y_log10() +
      labs(x="Rank", y="Total") +
      geom_line(aes(y = fitted), color="red", linetype="solid") +
      geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed") +
      geom_hline(yintercept=metadata(br.out)$inflection, color="forestgreen", linetype="dashed") +
      theme(legend.position="bottomleft") +
      guides(colour=guide_legend(override.aes=list(linetype=c("dashed", "dashed")))) +
      annotate("text", x=Inf, y=metadata(br.out)$knee, label="knee", color="dodgerblue", vjust=-1, hjust=1) +
      annotate("text", x=Inf, y=metadata(br.out)$inflection, label="inflection", color="forestgreen", vjust=-1, hjust=1)
  #Summary table -> may want to add a bunch of other summary metrics 
    #Extract counts with checks for NULL
    cell_counts <- as.data.frame(table(E170_gene_object@meta.data$is.cell))
    count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
                                      cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
    count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
                                       cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  #Create the summary table
    summary_table <- data.frame(
      Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
      Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops))
    summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  #Create the combined plot
    plot1 <- grid.arrange(
      rankplot,
      DimPlot(E170_gene_object, reduction = "umap", group.by = 'is.cell') + 
        labs(color = "is.cell", title = 'Seurat Object') + 
        theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
      FeaturePlot(E170_gene_object, features = "nCount_RNA") + theme(plot.background = element_rect(fill = "white")),
      FeaturePlot(E170_gene_object, features = "nFeature_RNA") + theme(plot.background = element_rect(fill = "white")),
      summary_grob, ncol = 2, top = textGrob('Empty drops vs real cells'))
  #Subset the Seurat object to remove cells marked as empty drops
    E170_gene_object_rm_empty <- subset(E170_gene_object, subset = is.cell == 'Cells')


#E191 Empty Drops Analysis
  E170_gene_object <- E191_gene_object
  E170_combined_df <- E191_combined_df
  #Define function to make dgCMatrix from combined counts
  makedgcmatrix <- function(count.matrix) {
    E170_gene_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell")
    list(E170_gene_object[["RNA"]]$counts)}
  E170_combined_df[] <- lapply(E170_combined_df, function(x) as.numeric(as.character(x)))
  outs.ddcmatrix <- makedgcmatrix(E170_combined_df)[[1]]
  br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
  e.out <- emptyDrops(outs.ddcmatrix, lower=100, niters=10000, test.ambient=TRUE, BPPARAM=SerialParam())
  is.cell <- e.out$FDR < fdr_threshold
  #Create a dataframe with FDR of TRUE cells
  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
  #Function for retrieving the Seurat cells and cluster in dataframe
  overlap_true_cell <- function(E170_gene_object) {
    seurat_cluster.df <- as.data.frame(E170_gene_object$seurat_clusters)
    seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
    seurat_cluster.df}
  #Obtain cluster dataframe from Seurat object
  overlap_CR <- overlap_true_cell(E170_gene_object)
  #Check overlaps between Seurat object and true cells
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
  #Function to add metadata to Seurat object
  True.cells <- function(e.out) {
    cells <- as.data.frame(e.out@rownames)
    fdr <- as.data.frame(e.out$FDR)
    T.F.cells <- cbind(cells, fdr)
    T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
    setnames(T.F.cells, c('FDR'))
    T.F.cells %>%
      mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))}
  cells_CR <- True.cells(e.out)
  E170_gene_object <- AddMetaData(E170_gene_object, metadata=cells_CR, col.name='is.cell')
  #Plot Empty drops on Gene UMAP
  #Create a ggplot object
  rankplot <- ggplot(br.out, aes(x=rank, y=total)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x="Rank", y="Total") +
    geom_line(aes(y = fitted), color="red", linetype="solid") +
    geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed") +
    geom_hline(yintercept=metadata(br.out)$inflection, color="forestgreen", linetype="dashed") +
    theme(legend.position="bottomleft") +
    guides(colour=guide_legend(override.aes=list(linetype=c("dashed", "dashed")))) +
    annotate("text", x=Inf, y=metadata(br.out)$knee, label="knee", color="dodgerblue", vjust=-1, hjust=1) +
    annotate("text", x=Inf, y=metadata(br.out)$inflection, label="inflection", color="forestgreen", vjust=-1, hjust=1)
  #Summary table -> may want to add a bunch of other summary metrics 
  #Extract counts with checks for NULL
  cell_counts <- as.data.frame(table(E170_gene_object@meta.data$is.cell))
  count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
                             cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
  count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
                              cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  #Create the summary table
  summary_table <- data.frame(
    Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops))
  summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  #Create the combined plot
  plot1 <- grid.arrange(
    rankplot,
    DimPlot(E170_gene_object, reduction = "umap", group.by = 'is.cell') + 
      labs(color = "is.cell", title = 'Seurat Object') + 
      theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nCount_RNA") + theme(plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nFeature_RNA") + theme(plot.background = element_rect(fill = "white")),
    summary_grob, ncol = 2, top = textGrob('Empty drops vs real cells'))
  #Subset the Seurat object to remove cells marked as empty drops
  E191_gene_object_rm_empty <- subset(E170_gene_object, subset = is.cell == 'Cells')

#E226 Empty Drops Analysis
  E170_gene_object <- E226_gene_object
  E170_combined_df <- E226_combined_df
  #Define function to make dgCMatrix from combined counts
  makedgcmatrix <- function(count.matrix) {
    E170_gene_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell")
    list(E170_gene_object[["RNA"]]$counts)}
  E170_combined_df[] <- lapply(E170_combined_df, function(x) as.numeric(as.character(x)))
  outs.ddcmatrix <- makedgcmatrix(E170_combined_df)[[1]]
  br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
  e.out <- emptyDrops(outs.ddcmatrix, lower=100, niters=10000, test.ambient=TRUE, BPPARAM=SerialParam())
  is.cell <- e.out$FDR < fdr_threshold
  #Create a dataframe with FDR of TRUE cells
  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
  #Function for retrieving the Seurat cells and cluster in dataframe
  overlap_true_cell <- function(E170_gene_object) {
    seurat_cluster.df <- as.data.frame(E170_gene_object$seurat_clusters)
    seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
    seurat_cluster.df}
  #Obtain cluster dataframe from Seurat object
  overlap_CR <- overlap_true_cell(E170_gene_object)
  #Check overlaps between Seurat object and true cells
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
  #Function to add metadata to Seurat object
  True.cells <- function(e.out) {
    cells <- as.data.frame(e.out@rownames)
    fdr <- as.data.frame(e.out$FDR)
    T.F.cells <- cbind(cells, fdr)
    T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
    setnames(T.F.cells, c('FDR'))
    T.F.cells %>%
      mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))}
  cells_CR <- True.cells(e.out)
  E170_gene_object <- AddMetaData(E170_gene_object, metadata=cells_CR, col.name='is.cell')
  #Plot Empty drops on Gene UMAP
  #Create a ggplot object
  rankplot <- ggplot(br.out, aes(x=rank, y=total)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x="Rank", y="Total") +
    geom_line(aes(y = fitted), color="red", linetype="solid") +
    geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed") +
    geom_hline(yintercept=metadata(br.out)$inflection, color="forestgreen", linetype="dashed") +
    theme(legend.position="bottomleft") +
    guides(colour=guide_legend(override.aes=list(linetype=c("dashed", "dashed")))) +
    annotate("text", x=Inf, y=metadata(br.out)$knee, label="knee", color="dodgerblue", vjust=-1, hjust=1) +
    annotate("text", x=Inf, y=metadata(br.out)$inflection, label="inflection", color="forestgreen", vjust=-1, hjust=1)
  #Summary table -> may want to add a bunch of other summary metrics 
  #Extract counts with checks for NULL
  cell_counts <- as.data.frame(table(E170_gene_object@meta.data$is.cell))
  count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
                             cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
  count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
                              cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  #Create the summary table
  summary_table <- data.frame(
    Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops))
  summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  #Create the combined plot
  plot1 <- grid.arrange(
    rankplot,
    DimPlot(E170_gene_object, reduction = "umap", group.by = 'is.cell') + 
      labs(color = "is.cell", title = 'Seurat Object') + 
      theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nCount_RNA") + theme(plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nFeature_RNA") + theme(plot.background = element_rect(fill = "white")),
    summary_grob, ncol = 2, top = textGrob('Empty drops vs real cells'))
  #Subset the Seurat object to remove cells marked as empty drops
  E226_gene_object_rm_empty <- subset(E170_gene_object, subset = is.cell == 'Cells')

#E231 Empty Drops Analysis
  E170_gene_object <- E231_gene_object
  E170_combined_df <- E231_combined_df
  #Define function to make dgCMatrix from combined counts
  makedgcmatrix <- function(count.matrix) {
    E170_gene_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell")
    list(E170_gene_object[["RNA"]]$counts)}
  E170_combined_df[] <- lapply(E170_combined_df, function(x) as.numeric(as.character(x)))
  outs.ddcmatrix <- makedgcmatrix(E170_combined_df)[[1]]
  br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
  e.out <- emptyDrops(outs.ddcmatrix, lower=100, niters=10000, test.ambient=TRUE, BPPARAM=SerialParam())
  is.cell <- e.out$FDR < fdr_threshold
  #Create a dataframe with FDR of TRUE cells
  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
  #Function for retrieving the Seurat cells and cluster in dataframe
  overlap_true_cell <- function(E170_gene_object) {
    seurat_cluster.df <- as.data.frame(E170_gene_object$seurat_clusters)
    seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
    seurat_cluster.df}
  #Obtain cluster dataframe from Seurat object
  overlap_CR <- overlap_true_cell(E170_gene_object)
  #Check overlaps between Seurat object and true cells
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
  #Function to add metadata to Seurat object
  True.cells <- function(e.out) {
    cells <- as.data.frame(e.out@rownames)
    fdr <- as.data.frame(e.out$FDR)
    T.F.cells <- cbind(cells, fdr)
    T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
    setnames(T.F.cells, c('FDR'))
    T.F.cells %>%
      mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))}
  cells_CR <- True.cells(e.out)
  E170_gene_object <- AddMetaData(E170_gene_object, metadata=cells_CR, col.name='is.cell')
  #Plot Empty drops on Gene UMAP
  #Create a ggplot object
  rankplot <- ggplot(br.out, aes(x=rank, y=total)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x="Rank", y="Total") +
    geom_line(aes(y = fitted), color="red", linetype="solid") +
    geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed") +
    geom_hline(yintercept=metadata(br.out)$inflection, color="forestgreen", linetype="dashed") +
    theme(legend.position="bottomleft") +
    guides(colour=guide_legend(override.aes=list(linetype=c("dashed", "dashed")))) +
    annotate("text", x=Inf, y=metadata(br.out)$knee, label="knee", color="dodgerblue", vjust=-1, hjust=1) +
    annotate("text", x=Inf, y=metadata(br.out)$inflection, label="inflection", color="forestgreen", vjust=-1, hjust=1)
  #Summary table -> may want to add a bunch of other summary metrics 
  #Extract counts with checks for NULL
  cell_counts <- as.data.frame(table(E170_gene_object@meta.data$is.cell))
  count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
                             cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
  count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
                              cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  #Create the summary table
  summary_table <- data.frame(
    Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops))
  summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  #Create the combined plot
  plot1 <- grid.arrange(
    rankplot,
    DimPlot(E170_gene_object, reduction = "umap", group.by = 'is.cell') + 
      labs(color = "is.cell", title = 'Seurat Object') + 
      theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nCount_RNA") + theme(plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nFeature_RNA") + theme(plot.background = element_rect(fill = "white")),
    summary_grob, ncol = 2, top = textGrob('Empty drops vs real cells'))
  #Subset the Seurat object to remove cells marked as empty drops
  E231_gene_object_rm_empty <- subset(E170_gene_object, subset = is.cell == 'Cells')

  
#E333 Empty Drops Analysis
  E170_gene_object <- E333_gene_object
  E170_combined_df <- E333_combined_df
  #Define function to make dgCMatrix from combined counts
  makedgcmatrix <- function(count.matrix) {
    E170_gene_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell")
    list(E170_gene_object[["RNA"]]$counts)}
  E170_combined_df[] <- lapply(E170_combined_df, function(x) as.numeric(as.character(x)))
  outs.ddcmatrix <- makedgcmatrix(E170_combined_df)[[1]]
  br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
  e.out <- emptyDrops(outs.ddcmatrix, lower=100, niters=10000, test.ambient=TRUE, BPPARAM=SerialParam())
  is.cell <- e.out$FDR < fdr_threshold
  #Create a dataframe with FDR of TRUE cells
  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
  #Function for retrieving the Seurat cells and cluster in dataframe
  overlap_true_cell <- function(E170_gene_object) {
    seurat_cluster.df <- as.data.frame(E170_gene_object$seurat_clusters)
    seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
    seurat_cluster.df}
  #Obtain cluster dataframe from Seurat object
  overlap_CR <- overlap_true_cell(E170_gene_object)
  #Check overlaps between Seurat object and true cells
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
  #Function to add metadata to Seurat object
  True.cells <- function(e.out) {
    cells <- as.data.frame(e.out@rownames)
    fdr <- as.data.frame(e.out$FDR)
    T.F.cells <- cbind(cells, fdr)
    T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
    setnames(T.F.cells, c('FDR'))
    T.F.cells %>%
      mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))}
  cells_CR <- True.cells(e.out)
  E170_gene_object <- AddMetaData(E170_gene_object, metadata=cells_CR, col.name='is.cell')
  #Plot Empty drops on Gene UMAP
  #Create a ggplot object
  rankplot <- ggplot(br.out, aes(x=rank, y=total)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x="Rank", y="Total") +
    geom_line(aes(y = fitted), color="red", linetype="solid") +
    geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed") +
    geom_hline(yintercept=metadata(br.out)$inflection, color="forestgreen", linetype="dashed") +
    theme(legend.position="bottomleft") +
    guides(colour=guide_legend(override.aes=list(linetype=c("dashed", "dashed")))) +
    annotate("text", x=Inf, y=metadata(br.out)$knee, label="knee", color="dodgerblue", vjust=-1, hjust=1) +
    annotate("text", x=Inf, y=metadata(br.out)$inflection, label="inflection", color="forestgreen", vjust=-1, hjust=1)
  #Summary table -> may want to add a bunch of other summary metrics 
  #Extract counts with checks for NULL
  cell_counts <- as.data.frame(table(E170_gene_object@meta.data$is.cell))
  count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
                             cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
  count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
                              cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  #Create the summary table
  summary_table <- data.frame(
    Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops))
  summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  #Create the combined plot
  plot1 <- grid.arrange(
    rankplot,
    DimPlot(E170_gene_object, reduction = "umap", group.by = 'is.cell') + 
      labs(color = "is.cell", title = 'Seurat Object') + 
      theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nCount_RNA") + theme(plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nFeature_RNA") + theme(plot.background = element_rect(fill = "white")),
    summary_grob, ncol = 2, top = textGrob('Empty drops vs real cells'))
  #Subset the Seurat object to remove cells marked as empty drops
  E333_gene_object_rm_empty <- subset(E170_gene_object, subset = is.cell == 'Cells')

#E435 Empty Drops Analysis
  E170_gene_object <- E435_gene_object
  E170_combined_df <- E435_combined_df
  #Define function to make dgCMatrix from combined counts
  makedgcmatrix <- function(count.matrix) {
    E170_gene_object <- CreateSeuratObject(counts = count.matrix, project = "singlecell")
    list(E170_gene_object[["RNA"]]$counts)}
  E170_combined_df[] <- lapply(E170_combined_df, function(x) as.numeric(as.character(x)))
  outs.ddcmatrix <- makedgcmatrix(E170_combined_df)[[1]]
  br.out <- DropletUtils::barcodeRanks(outs.ddcmatrix)
  e.out <- emptyDrops(outs.ddcmatrix, lower=100, niters=10000, test.ambient=TRUE, BPPARAM=SerialParam())
  is.cell <- e.out$FDR < fdr_threshold
  #Create a dataframe with FDR of TRUE cells
  is.true.cell_CR <- as.data.frame(e.out@listData[["FDR"]], e.out@rownames)
  is.true.cell_CR <- is.true.cell_CR %>% filter(is.true.cell_CR$`e.out@listData[["FDR"]]` <= fdr_threshold)
  is.true.cell_CR <- tibble::rownames_to_column(is.true.cell_CR, "cell_id")
  #Function for retrieving the Seurat cells and cluster in dataframe
  overlap_true_cell <- function(E170_gene_object) {
    seurat_cluster.df <- as.data.frame(E170_gene_object$seurat_clusters)
    seurat_cluster.df <- tibble::rownames_to_column(seurat_cluster.df, "cell_id")
    seurat_cluster.df}
  #Obtain cluster dataframe from Seurat object
  overlap_CR <- overlap_true_cell(E170_gene_object)
  #Check overlaps between Seurat object and true cells
  summary(overlap_CR$cell_id %in% is.true.cell_CR$cell_id)
  #Function to add metadata to Seurat object
  True.cells <- function(e.out) {
    cells <- as.data.frame(e.out@rownames)
    fdr <- as.data.frame(e.out$FDR)
    T.F.cells <- cbind(cells, fdr)
    T.F.cells <- data.frame(T.F.cells[,-1], row.names = T.F.cells[,1])
    setnames(T.F.cells, c('FDR'))
    T.F.cells %>%
      mutate(FDR = case_when(FDR < fdr_threshold ~ "Cells", FDR > fdr_threshold ~ "Empty_drops"))}
  cells_CR <- True.cells(e.out)
  E170_gene_object <- AddMetaData(E170_gene_object, metadata=cells_CR, col.name='is.cell')
  #Plot Empty drops on Gene UMAP
  #Create a ggplot object
  rankplot <- ggplot(br.out, aes(x=rank, y=total)) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x="Rank", y="Total") +
    geom_line(aes(y = fitted), color="red", linetype="solid") +
    geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed") +
    geom_hline(yintercept=metadata(br.out)$inflection, color="forestgreen", linetype="dashed") +
    theme(legend.position="bottomleft") +
    guides(colour=guide_legend(override.aes=list(linetype=c("dashed", "dashed")))) +
    annotate("text", x=Inf, y=metadata(br.out)$knee, label="knee", color="dodgerblue", vjust=-1, hjust=1) +
    annotate("text", x=Inf, y=metadata(br.out)$inflection, label="inflection", color="forestgreen", vjust=-1, hjust=1)
  #Summary table -> may want to add a bunch of other summary metrics 
  #Extract counts with checks for NULL
  cell_counts <- as.data.frame(table(E170_gene_object@meta.data$is.cell))
  count_true_cells <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Cells"]) > 0,
                             cell_counts$Freq[cell_counts$Var1 == "Cells"], 0)
  count_empty_drops <- ifelse(length(cell_counts$Freq[cell_counts$Var1 == "Empty_drops"]) > 0,
                              cell_counts$Freq[cell_counts$Var1 == "Empty_drops"], 0)
  #Create the summary table
  summary_table <- data.frame(
    Description = c('fdr', 'lower Counts', 'number of true cells', 'number of empty drops'),
    Value = c(fdr_threshold, lower, count_true_cells, count_empty_drops))
  summary_grob <- tableGrob(summary_table, rows = NULL, cols = NULL)
  #Create the combined plot
  plot1 <- grid.arrange(
    rankplot,
    DimPlot(E170_gene_object, reduction = "umap", group.by = 'is.cell') + 
      labs(color = "is.cell", title = 'Seurat Object') + 
      theme(text = element_text(size = 10), plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nCount_RNA") + theme(plot.background = element_rect(fill = "white")),
    FeaturePlot(E170_gene_object, features = "nFeature_RNA") + theme(plot.background = element_rect(fill = "white")),
    summary_grob, ncol = 2, top = textGrob('Empty drops vs real cells'))
  #Subset the Seurat object to remove cells marked as empty drops
  E435_gene_object_rm_empty <- subset(E170_gene_object, subset = is.cell == 'Cells')

  
###Low quality cell filtering (mitochondrial content & feature numbers)
E170_qc <- E170_gene_object_rm_empty
E191_qc <- E191_gene_object_rm_empty
E226_qc <- E226_gene_object_rm_empty
E231_qc <- E231_gene_object_rm_empty
E333_qc <- E333_gene_object_rm_empty
E435_qc <- E435_gene_object_rm_empty

#Plot relationship between reads and unique genes per cell
E170_plot_scatter1 <- FeatureScatter(E170_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() +
  labs(title = "E170 Reads vs Unique Genes per Cell BEFORE Filtering")
  plot(E170_plot_scatter1)
E191_plot_scatter1 <- FeatureScatter(E191_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() +
  labs(title = "E191 Reads vs Unique Genes per Cell BEFORE Filtering")
  plot(E191_plot_scatter1)
E226_plot_scatter1 <- FeatureScatter(E226_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() +
  labs(title = "E226 Reads vs Unique Genes per Cell BEFORE Filtering")
  plot(E226_plot_scatter1)
E231_plot_scatter1 <- FeatureScatter(E231_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() +
  labs(title = "E231 Reads vs Unique Genes per Cell BEFORE Filtering")
  plot(E231_plot_scatter1)
E333_plot_scatter1 <- FeatureScatter(E333_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() +
  labs(title = "E333 Reads vs Unique Genes per Cell BEFORE Filtering")
  plot(E333_plot_scatter1)
E435_plot_scatter1 <- FeatureScatter(E435_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + NoLegend() +
  labs(title = "E435 Reads vs Unique Genes per Cell BEFORE Filtering")
  plot(E435_plot_scatter1)

#Add mitochondrial percentage
E170_qc[["joined"]] <- JoinLayers(E170_qc[["RNA"]])
E170_qc[["percent.mt"]] <- PercentageFeatureSet(E170_qc, pattern = "^MT-")
E170_qc$group <- "All"
plot <- VlnPlot(E170_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group") 
plot + plot_annotation(title = "E170 QC plots (gene level) BEFORE Filtering")

E191_qc[["joined"]] <- JoinLayers(E191_qc[["RNA"]])
E191_qc[["percent.mt"]] <- PercentageFeatureSet(E191_qc, pattern = "^MT-")
E191_qc$group <- "All"
plot <- VlnPlot(E191_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group") 
plot + plot_annotation(title = "E191 QC plots (gene level) BEFORE Filtering")

E226_qc[["joined"]] <- JoinLayers(E226_qc[["RNA"]])
E226_qc[["percent.mt"]] <- PercentageFeatureSet(E226_qc, pattern = "^MT-")
E226_qc$group <- "All"
plot <- VlnPlot(E226_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group") 
plot + plot_annotation(title = "E226 QC plots (gene level) BEFORE Filtering")

E231_qc[["joined"]] <- JoinLayers(E231_qc[["RNA"]])
E231_qc[["percent.mt"]] <- PercentageFeatureSet(E231_qc, pattern = "^MT-")
E231_qc$group <- "All"
plot <- VlnPlot(E231_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group") 
plot + plot_annotation(title = "E231 QC plots (gene level) BEFORE Filtering")

E333_qc[["joined"]] <- JoinLayers(E333_qc[["RNA"]])
E333_qc[["percent.mt"]] <- PercentageFeatureSet(E333_qc, pattern = "^MT-")
E333_qc$group <- "All"
plot <- VlnPlot(E333_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group") 
plot + plot_annotation(title = "E333 QC plots (gene level) BEFORE Filtering")

E435_qc[["joined"]] <- JoinLayers(E435_qc[["RNA"]])
E435_qc[["percent.mt"]] <- PercentageFeatureSet(E435_qc, pattern = "^MT-")
E435_qc$group <- "All"
plot <- VlnPlot(E435_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group") 
plot  + plot_annotation(title = "E435 QC plots (gene level) BEFORE Filtering")

#Filter cells based on feature and count thresholds - (change these based on your data)
  #max.features = 
  min.features = 200
  min.counts = 200
  #max.counts = 
  MT = 20
  
  #Filter the seurat object based on the QC params listed above 
  filtered_E170 <- subset(E170_qc, subset = nFeature_RNA > min.features & percent.mt < MT & nCount_RNA > min.counts)
  filtered_E191 <- subset(E191_qc, subset = nFeature_RNA > min.features & percent.mt < MT & nCount_RNA > min.counts)
  filtered_E226 <- subset(E226_qc, subset = nFeature_RNA > min.features & percent.mt < MT & nCount_RNA > min.counts)
  filtered_E231 <- subset(E231_qc, subset = nFeature_RNA > min.features & percent.mt < MT & nCount_RNA > min.counts)
  filtered_E333 <- subset(E333_qc, subset = nFeature_RNA > min.features & percent.mt < MT & nCount_RNA > min.counts)
  filtered_E435 <- subset(E435_qc, subset = nFeature_RNA > min.features & percent.mt < MT & nCount_RNA > min.counts)
  
  #Plot quality metrics after filtering
    plot <- VlnPlot(filtered_E170, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group")
    plot + plot_annotation(title = "E170 QC metrics gene level AFTER Filtering")
    plot <- VlnPlot(filtered_E191, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group")
    plot + plot_annotation(title = "E191 QC metrics gene level AFTER Filtering")
    plot <- VlnPlot(filtered_E226, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group")
    plot + plot_annotation(title = "E226 QC metrics gene level AFTER Filtering")
    plot <- VlnPlot(filtered_E231, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group")
    plot + plot_annotation(title = "E231 QC metrics gene level AFTER Filtering")
    plot <- VlnPlot(filtered_E333, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group")
    plot + plot_annotation(title = "E333 QC metrics gene level AFTER Filtering")
    plot <- VlnPlot(filtered_E435, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="group")
    plot + plot_annotation(title = "E435 QC metrics gene level AFTER Filtering")
    
  #Normalize data
  filtered_E170 <- NormalizeData(filtered_E170, normalization.method = "LogNormalize", scale.factor = 10000)
  filtered_E191 <- NormalizeData(filtered_E191, normalization.method = "LogNormalize", scale.factor = 10000)
  filtered_E226 <- NormalizeData(filtered_E226, normalization.method = "LogNormalize", scale.factor = 10000)
  filtered_E231 <- NormalizeData(filtered_E231, normalization.method = "LogNormalize", scale.factor = 10000)
  filtered_E333 <- NormalizeData(filtered_E333, normalization.method = "LogNormalize", scale.factor = 10000)
  filtered_E435 <- NormalizeData(filtered_E435, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #Identify highly variable features
  filtered_E170 <- FindVariableFeatures(filtered_E170, selection.method = "vst", nfeatures = 2000)
  filtered_E191 <- FindVariableFeatures(filtered_E191, selection.method = "vst", nfeatures = 2000)
  filtered_E226 <- FindVariableFeatures(filtered_E226, selection.method = "vst", nfeatures = 2000)
  filtered_E231 <- FindVariableFeatures(filtered_E231, selection.method = "vst", nfeatures = 2000)
  filtered_E333 <- FindVariableFeatures(filtered_E333, selection.method = "vst", nfeatures = 2000)
  filtered_E435 <- FindVariableFeatures(filtered_E435, selection.method = "vst", nfeatures = 2000)
  
  #Apply linear transformation
  E170_all_genes <- rownames(filtered_E170)
  filtered_E170 <- ScaleData(filtered_E170, features = E170_all_genes)
  E191_all_genes <- rownames(filtered_E191)
  filtered_E191 <- ScaleData(filtered_E191, features = E191_all_genes)
  E226_all_genes <- rownames(filtered_E226)
  filtered_E226 <- ScaleData(filtered_E226, features = E226_all_genes)
  E231_all_genes <- rownames(filtered_E231)
  filtered_E231 <- ScaleData(filtered_E231, features = E231_all_genes)
  E333_all_genes <- rownames(filtered_E333)
  filtered_E333 <- ScaleData(filtered_E333, features = E333_all_genes)
  E435_all_genes <- rownames(filtered_E435)
  filtered_E435 <- ScaleData(filtered_E435, features = E435_all_genes)
  
  #Perform PCA
  filtered_E170 <- RunPCA(filtered_E170, features = VariableFeatures(object =  filtered_E170))
  filtered_E191 <- RunPCA(filtered_E191, features = VariableFeatures(object =  filtered_E191))
  filtered_E226 <- RunPCA(filtered_E226, features = VariableFeatures(object =  filtered_E226))
  filtered_E231 <- RunPCA(filtered_E231, features = VariableFeatures(object =  filtered_E231))
  filtered_E333 <- RunPCA(filtered_E333, features = VariableFeatures(object =  filtered_E333))
  filtered_E435 <- RunPCA(filtered_E435, features = VariableFeatures(object =  filtered_E435))
  
  #Cluster cells
  filtered_E170 <- FindNeighbors(filtered_E170, dims=1:15)
  filtered_E170 <- FindClusters(filtered_E170, resolution=0.2)
  filtered_E191 <- FindNeighbors(filtered_E191, dims=1:15)
  filtered_E191 <- FindClusters(filtered_E191, resolution=0.2)
  filtered_E226 <- FindNeighbors(filtered_E226, dims=1:15)
  filtered_E226 <- FindClusters(filtered_E226, resolution=0.2)
  filtered_E231 <- FindNeighbors(filtered_E231, dims=1:15)
  filtered_E231 <- FindClusters(filtered_E231, resolution=0.2)
  filtered_E333 <- FindNeighbors(filtered_E333, dims=1:15)
  filtered_E333 <- FindClusters(filtered_E333, resolution=0.2)
  filtered_E435 <- FindNeighbors(filtered_E435, dims=1:15)
  filtered_E435 <- FindClusters(filtered_E435, resolution=0.2)

  
###Doublet Filtering
#E170
    #pK Identification (no ground-truth) 
    sweep.res.list_E170 <- paramSweep(filtered_E170, PCs=1:15, sct=FALSE)
    sweep.stats_E170 <- summarizeSweep(sweep.res.list_E170, GT=FALSE)
    bcmvn_E170 <- find.pK(sweep.stats_E170)
    pK_E170 <- bcmvn_E170 %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
    pK_E170 <- as.numeric(as.character(pK_E170[[1]]))
    #Homotypic Doublet Proportion Estimate
    annotations_E170 <- filtered_E170@meta.data$seurat_clusters
    homotypic.prop_E170 <- modelHomotypic(annotations_E170)
    nExp_poi_E170 <- round(0.016*nrow(filtered_E170@meta.data))
    nExp_poi.adj_E170 <- round(nExp_poi_E170*(1-homotypic.prop_E170))
    #Run doubletFinder 
    E170_gene_object_doublet <- doubletFinder(filtered_E170, PCs=1:15, pN=0.25, pK=pK_E170, nExp=nExp_poi.adj_E170, reuse.pANN=FALSE, sct=FALSE)
    #Summary doublets
    colnames(E170_gene_object_doublet@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(E170_gene_object_doublet@meta.data))
    statsDoublets_E170 <- E170_gene_object_doublet@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA = median(nCount_RNA), Median_nFeature_RNA = median(nFeature_RNA), Count = n())
    #Visualize doublets
    DimPlot(E170_gene_object_doublet, reduction='umap', group.by = "DF.classifications")

#E191
    #pK Identification (no ground-truth) 
    sweep.res.list_E191 <- paramSweep(filtered_E191, PCs=1:15, sct=FALSE)
    sweep.stats_E191 <- summarizeSweep(sweep.res.list_E191, GT=FALSE)
    bcmvn_E191 <- find.pK(sweep.stats_E191)
    pK_E191 <- bcmvn_E191 %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
    pK_E191 <- as.numeric(as.character(pK_E191[[1]]))
    #Homotypic Doublet Proportion Estimate
    annotations_E191 <- filtered_E191@meta.data$seurat_clusters
    homotypic.prop_E191 <- modelHomotypic(annotations_E191)
    nExp_poi_E191 <- round(0.016*nrow(filtered_E191@meta.data))
    nExp_poi.adj_E191 <- round(nExp_poi_E191*(1-homotypic.prop_E191))
    #Run doubletFinder 
    E191_gene_object_doublet <- doubletFinder(filtered_E191, PCs=1:15, pN=0.25, pK=pK_E191, nExp=nExp_poi.adj_E191, reuse.pANN=FALSE, sct=FALSE)
    #Summary doublets
    colnames(E191_gene_object_doublet@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(E191_gene_object_doublet@meta.data))
    statsDoublets_E191 <- E191_gene_object_doublet@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA = median(nCount_RNA), Median_nFeature_RNA=median(nFeature_RNA), Count=n())
    #Visualize doublets
    DimPlot(E191_gene_object_doublet, reduction='umap', group.by="DF.classifications")

#E226
    #pK Identification (no ground-truth) 
    sweep.res.list_E226 <- paramSweep(filtered_E226, PCs=1:15, sct=FALSE)
    sweep.stats_E226 <- summarizeSweep(sweep.res.list_E226, GT=FALSE)
    bcmvn_E226 <- find.pK(sweep.stats_E226)
    pK_E226 <- bcmvn_E226 %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
    pK_E226 <- as.numeric(as.character(pK_E226[[1]]))
    #Homotypic Doublet Proportion Estimate
    annotations_E226 <- filtered_E226@meta.data$seurat_clusters
    homotypic.prop_E226 <- modelHomotypic(annotations_E226)
    nExp_poi_E226 <- round(0.016*nrow(filtered_E226@meta.data))
    nExp_poi.adj_E226 <- round(nExp_poi_E226*(1-homotypic.prop_E226))
    #Run doubletFinder 
    E226_gene_object_doublet <- doubletFinder(filtered_E226, PCs=1:15, pN=0.25, pK=pK_E226, nExp=nExp_poi.adj_E226, reuse.pANN=FALSE, sct=FALSE)
    #Summary doublets
    colnames(E226_gene_object_doublet@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(E226_gene_object_doublet@meta.data))
    statsDoublets_E226 <- E226_gene_object_doublet@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA = median(nCount_RNA), Median_nFeature_RNA=median(nFeature_RNA), Count=n())
    #Visualize doublets
    DimPlot(E226_gene_object_doublet, reduction='umap', group.by="DF.classifications")
        
#E231
    #pK Identification (no ground-truth) 
    sweep.res.list_E231 <- paramSweep(filtered_E231, PCs=1:15, sct=FALSE)
    sweep.stats_E231 <- summarizeSweep(sweep.res.list_E231, GT=FALSE)
    bcmvn_E231 <- find.pK(sweep.stats_E231)
    pK_E231 <- bcmvn_E231 %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
    pK_E231 <- as.numeric(as.character(pK_E231[[1]]))
    #Homotypic Doublet Proportion Estimate
    annotations_E231 <- filtered_E231@meta.data$seurat_clusters
    homotypic.prop_E231 <- modelHomotypic(annotations_E231)
    nExp_poi_E231 <- round(0.016*nrow(filtered_E231@meta.data))
    nExp_poi.adj_E231 <- round(nExp_poi_E231*(1-homotypic.prop_E231))
    #Run doubletFinder 
    E231_gene_object_doublet <- doubletFinder(filtered_E231, PCs=1:15, pN=0.25, pK=pK_E231, nExp=nExp_poi.adj_E231, reuse.pANN=FALSE, sct=FALSE)
    #Summary doublets
    colnames(E231_gene_object_doublet@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(E231_gene_object_doublet@meta.data))
    statsDoublets_E231 <- E231_gene_object_doublet@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA = median(nCount_RNA), Median_nFeature_RNA=median(nFeature_RNA), Count=n())
    #Visualize doublets
    DimPlot(E231_gene_object_doublet, reduction='umap', group.by="DF.classifications")

#E333
    #pK Identification (no ground-truth) 
    sweep.res.list_E333 <- paramSweep(filtered_E333, PCs=1:15, sct=FALSE)
    sweep.stats_E333 <- summarizeSweep(sweep.res.list_E333, GT=FALSE)
    bcmvn_E333 <- find.pK(sweep.stats_E333)
    pK_E333 <- bcmvn_E333 %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
    pK_E333 <- as.numeric(as.character(pK_E333[[1]]))
    #Homotypic Doublet Proportion Estimate
    annotations_E333 <- filtered_E333@meta.data$seurat_clusters
    homotypic.prop_E333 <- modelHomotypic(annotations_E333)
    nExp_poi_E333 <- round(0.016*nrow(filtered_E333@meta.data))
    nExp_poi.adj_E333 <- round(nExp_poi_E333*(1-homotypic.prop_E333))
    #Run doubletFinder 
    E333_gene_object_doublet <- doubletFinder(filtered_E333, PCs=1:15, pN=0.25, pK=pK_E333, nExp=nExp_poi.adj_E333, reuse.pANN=FALSE, sct=FALSE)
    #Summary doublets
    colnames(E333_gene_object_doublet@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(E333_gene_object_doublet@meta.data))
    statsDoublets_E333 <- E333_gene_object_doublet@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA = median(nCount_RNA), Median_nFeature_RNA=median(nFeature_RNA), Count=n())
    #Visualize doublets
    DimPlot(E333_gene_object_doublet, reduction='umap', group.by="DF.classifications")

#E435
    #pK Identification (no ground-truth) 
    sweep.res.list_E435 <- paramSweep(filtered_E435, PCs=1:15, sct=FALSE)
    sweep.stats_E435 <- summarizeSweep(sweep.res.list_E435, GT=FALSE)
    bcmvn_E435 <- find.pK(sweep.stats_E435)
    pK_E435 <- bcmvn_E435 %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK) 
    pK_E435 <- as.numeric(as.character(pK_E435[[1]]))
    #Homotypic Doublet Proportion Estimate
    annotations_E435 <- filtered_E435@meta.data$seurat_clusters
    homotypic.prop_E435 <- modelHomotypic(annotations_E435)
    nExp_poi_E435 <- round(0.016*nrow(filtered_E435@meta.data))
    nExp_poi.adj_E435 <- round(nExp_poi_E435*(1-homotypic.prop_E435))
    #Run doubletFinder 
    E435_gene_object_doublet <- doubletFinder(filtered_E435, PCs=1:15, pN=0.25, pK=pK_E435, nExp=nExp_poi.adj_E435, reuse.pANN=FALSE, sct=FALSE)
    #Summary doublets
    colnames(E435_gene_object_doublet@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(E435_gene_object_doublet@meta.data))
    statsDoublets_E435 <- E435_gene_object_doublet@meta.data %>%
    group_by(DF.classifications) %>%
    summarize(Median_nCount_RNA=median(nCount_RNA), Median_nFeature_RNA=median(nFeature_RNA), Count=n())
    #Visualize doublets
    DimPlot(E435_gene_object_doublet, reduction='umap', group.by="DF.classifications")

#Save the seurat object with doublets listed 
E170_gene_object_doublets <- E170_gene_object_doublet
E191_gene_object_doublets <- E191_gene_object_doublet
E226_gene_object_doublets <- E226_gene_object_doublet
E231_gene_object_doublets <- E231_gene_object_doublet
E333_gene_object_doublets <- E333_gene_object_doublet
E435_gene_object_doublets <- E435_gene_object_doublet
    
#Generate new seurat objects excluding the doublets
E170_genes_filtered <- subset(E170_gene_object_doublet, subset = DF.classifications == 'Singlet')
E191_genes_filtered <- subset(E191_gene_object_doublet, subset = DF.classifications == 'Singlet')
E226_genes_filtered <- subset(E226_gene_object_doublet, subset = DF.classifications == 'Singlet')
E231_genes_filtered <- subset(E231_gene_object_doublet, subset = DF.classifications == 'Singlet')
E333_genes_filtered <- subset(E333_gene_object_doublet, subset = DF.classifications == 'Singlet')
E435_genes_filtered <- subset(E435_gene_object_doublet, subset = DF.classifications == 'Singlet')

#Generate figures 
#E170
  ggplot_list <- list(
    ElbowPlot(E170_genes_filtered) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10)),
    FeatureScatter(E170_genes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering"),
    DimPlot(E170_genes_filtered, reduction = "umap") + 
        labs(color = "Cluster \n(from PCA)", title = '') + 
        theme(text = element_text(size = 10)),
    FeaturePlot(E170_genes_filtered, reduction = "umap", features = 'nCount_RNA') + 
        labs(color = "UMI count", title = '') + 
        theme(text = element_text(size = 10)),
    FeaturePlot(E170_genes_filtered, reduction = "umap", features = 'nFeature_RNA') + 
        labs(color = str_wrap("Feature count (gene)", 15), title = '') + 
        theme(text = element_text(size = 10)))
    E170_combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)
    plot(E170_combined_plots)
    plot(DimPlot(E170_gene_object_doublets, reduction = 'umap', group.by = "DF.classifications")) +
      ggtitle("E170 DF.classifications")
    #Table Summary of doublets and singlets
      E170_tbl_sts1 <- tableGrob(statsDoublets_E170)
      grid.newpage()
      grid.draw(E170_tbl_sts1)
#E191
  ggplot_list <- list(
    ElbowPlot(E191_genes_filtered) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10)),
    FeatureScatter(E191_genes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering"),
    DimPlot(E191_genes_filtered, reduction = "umap") + 
        labs(color = "Cluster \n(from PCA)", title = '') + 
        theme(text = element_text(size = 10)),
    FeaturePlot(E191_genes_filtered, reduction = "umap", features = 'nCount_RNA') + 
        labs(color = "UMI count", title = '') + 
        theme(text = element_text(size = 10)),
    FeaturePlot(E191_genes_filtered, reduction = "umap", features = 'nFeature_RNA') + 
        labs(color = str_wrap("Feature count (gene)", 15), title = '') + 
        theme(text = element_text(size = 10)))
    E191_combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)
    plot(E191_combined_plots)
    plot(DimPlot(E191_gene_object_doublets, reduction = 'umap', group.by = "DF.classifications")) +
        ggtitle("E191 DF.classifications")
  #Table Summary of doublets and singlets
  E191_tbl_sts1 <- tableGrob(statsDoublets_E191)
  grid.newpage()
  grid.draw(E191_tbl_sts1)  
#E226
  ggplot_list <- list(
    ElbowPlot(E226_genes_filtered) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10)),
    FeatureScatter(E226_genes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering"),
    DimPlot(E226_genes_filtered, reduction = "umap") + 
      labs(color = "Cluster \n(from PCA)", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E226_genes_filtered, reduction = "umap", features = 'nCount_RNA') + 
      labs(color = "UMI count", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E226_genes_filtered, reduction = "umap", features = 'nFeature_RNA') + 
      labs(color = str_wrap("Feature count (gene)", 15), title = '') + 
      theme(text = element_text(size = 10)))
  E226_combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)
  plot(E226_combined_plots)
  plot(DimPlot(E226_gene_object_doublets, reduction = 'umap', group.by = "DF.classifications")) +
    ggtitle("E226 DF.classifications")
  #Table Summary of doublets and singlets
  E226_tbl_sts1 <- tableGrob(statsDoublets_E226)
  grid.newpage()
  grid.draw(E226_tbl_sts1) 
#E231
  ggplot_list <- list(
    ElbowPlot(E231_genes_filtered) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10)),
    FeatureScatter(E231_genes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering"),
    DimPlot(E231_genes_filtered, reduction = "umap") + 
      labs(color = "Cluster \n(from PCA)", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E231_genes_filtered, reduction = "umap", features = 'nCount_RNA') + 
      labs(color = "UMI count", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E231_genes_filtered, reduction = "umap", features = 'nFeature_RNA') + 
      labs(color = str_wrap("Feature count (gene)", 15), title = '') + 
      theme(text = element_text(size = 10)))
  E231_combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)
  plot(E231_combined_plots)
  plot(DimPlot(E231_gene_object_doublets, reduction = 'umap', group.by = "DF.classifications")) +
    ggtitle("E231 DF.classifications")
  #Table Summary of doublets and singlets
  E231_tbl_sts1 <- tableGrob(statsDoublets_E231)
  grid.newpage()
  grid.draw(E231_tbl_sts1)
#E333
  ggplot_list <- list(
    ElbowPlot(E333_genes_filtered) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10)),
    FeatureScatter(E333_genes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering"),
    DimPlot(E333_genes_filtered, reduction = "umap") + 
      labs(color = "Cluster \n(from PCA)", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E333_genes_filtered, reduction = "umap", features = 'nCount_RNA') + 
      labs(color = "UMI count", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E333_genes_filtered, reduction = "umap", features = 'nFeature_RNA') + 
      labs(color = str_wrap("Feature count (gene)", 15), title = '') + 
      theme(text = element_text(size = 10)))
  E333_combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)
  plot(E333_combined_plots)
  plot(DimPlot(E333_gene_object_doublets, reduction = 'umap', group.by = "DF.classifications")) +
    ggtitle("E333 DF.classifications")
  #Table Summary of doublets and singlets
  E333_tbl_sts1 <- tableGrob(statsDoublets_E333)
  grid.newpage()
  grid.draw(E333_tbl_sts1)
#E435
  ggplot_list <- list(
    ElbowPlot(E435_genes_filtered) + labs(title = 'SD explained by each PC') + theme(text = element_text(size = 10)),
    FeatureScatter(E435_genes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell AFTER filtering"),
    DimPlot(E435_genes_filtered, reduction = "umap") + 
      labs(color = "Cluster \n(from PCA)", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E435_genes_filtered, reduction = "umap", features = 'nCount_RNA') + 
      labs(color = "UMI count", title = '') + 
      theme(text = element_text(size = 10)),
    FeaturePlot(E435_genes_filtered, reduction = "umap", features = 'nFeature_RNA') + 
      labs(color = str_wrap("Feature count (gene)", 15), title = '') + 
      theme(text = element_text(size = 10)))
  E435_combined_plots <- plot_grid(plotlist = ggplot_list, ncol = 2)
  plot(E435_combined_plots)
  plot(DimPlot(E435_gene_object_doublets, reduction = 'umap', group.by = "DF.classifications")) +
    ggtitle("E435 DF.classifications")
  #Table Summary of doublets and singlets
  E435_tbl_sts1 <- tableGrob(statsDoublets_E435)
  grid.newpage()
  grid.draw(E435_tbl_sts1)

  
#Stats summary
#E170
E170_stats_sumary <- rbind("Sample ID" = E170_genes,
                           "Cells_before_filter" = dim(E170_qc)[2],
                           "Cells_after_filter" = dim(E170_gene_object_doublet)[2],
                           "Median Feature per Cell before filter" = median(E170_qc$nFeature_RNA),
                           "Median Reads per Gene before filter" = median(E170_qc$nCount_RNA),
                           "Median Feature per Cell" = median(E170_gene_object_doublet$nFeature_RNA),
                           "Median Reads per Gene" = median(E170_gene_object_doublet$nCount_RNA),
                           #"Max Features" = max.features,
                           "Min Features" = min.features,
                           "Min Counts" = min.counts,
                           #"Max Counts" = max.counts,
                           "MT Percentage" = MT,
                           "NPCs" = 15,
                           "Median Percent MT before Filter" = median(E170_qc@meta.data[["percent.mt"]]),
                           "Median Percent MT after Filter" = median(E170_gene_object_doublet@meta.data[["percent.mt"]]))
  E170_tbl_sts2 <- tableGrob(E170_stats_sumary, theme = ttheme_minimal())
  grid.newpage()
  grid.draw(E170_tbl_sts2)
#E191
E191_stats_sumary <- rbind("Sample ID" = E191_genes,
                           "Cells_before_filter" = dim(E191_qc)[2],
                           "Cells_after_filter" = dim(E191_gene_object_doublet)[2],
                           "Median Feature per Cell before filter" = median(E191_qc$nFeature_RNA),
                           "Median Reads per Gene before filter" = median(E191_qc$nCount_RNA),
                           "Median Feature per Cell" = median(E191_gene_object_doublet$nFeature_RNA),
                           "Median Reads per Gene" = median(E191_gene_object_doublet$nCount_RNA),
                           #"Max Features" = max.features,
                           "Min Features" = min.features,
                           "Min Counts" = min.counts,
                           #"Max Counts" = max.counts,
                           "MT Percentage" = MT,
                           "NPCs" = 15,
                           "Median Percent MT before Filter" = median(E191_qc@meta.data[["percent.mt"]]),
                           "Median Percent MT after Filter" = median(E191_gene_object_doublet@meta.data[["percent.mt"]]))
  E191_tbl_sts2 <- tableGrob(E191_stats_sumary, theme = ttheme_minimal())
  grid.newpage()
  grid.draw(E191_tbl_sts2)
#E226
E226_stats_sumary <- rbind("Sample ID" = E226_genes,
                           "Cells_before_filter" = dim(E226_qc)[2],
                           "Cells_after_filter" = dim(E226_gene_object_doublet)[2],
                           "Median Feature per Cell before filter" = median(E226_qc$nFeature_RNA),
                           "Median Reads per Gene before filter" = median(E226_qc$nCount_RNA),
                           "Median Feature per Cell" = median(E226_gene_object_doublet$nFeature_RNA),
                           "Median Reads per Gene" = median(E226_gene_object_doublet$nCount_RNA),
                           #"Max Features" = max.features,
                           "Min Features" = min.features,
                           "Min Counts" = min.counts,
                           #"Max Counts" = max.counts,
                           "MT Percentage" = MT,
                           "NPCs" = 15,
                           "Median Percent MT before Filter" = median(E226_qc@meta.data[["percent.mt"]]),
                           "Median Percent MT after Filter" = median(E226_gene_object_doublet@meta.data[["percent.mt"]]))
  E226_tbl_sts2 <- tableGrob(E226_stats_sumary, theme = ttheme_minimal())
  grid.newpage()
  grid.draw(E226_tbl_sts2)
#E231
E231_stats_sumary <- rbind("Sample ID" = E231_genes,
                           "Cells_before_filter" = dim(E231_qc)[2],
                           "Cells_after_filter" = dim(E231_gene_object_doublet)[2],
                           "Median Feature per Cell before filter" = median(E231_qc$nFeature_RNA),
                           "Median Reads per Gene before filter" = median(E231_qc$nCount_RNA),
                           "Median Feature per Cell" = median(E231_gene_object_doublet$nFeature_RNA),
                           "Median Reads per Gene" = median(E231_gene_object_doublet$nCount_RNA),
                           #"Max Features" = max.features,
                           "Min Features" = min.features,
                           "Min Counts" = min.counts,
                           #"Max Counts" = max.counts,
                           "MT Percentage" = MT,
                           "NPCs" = 15,
                           "Median Percent MT before Filter" = median(E231_qc@meta.data[["percent.mt"]]),
                           "Median Percent MT after Filter" = median(E231_gene_object_doublet@meta.data[["percent.mt"]]))
  E231_tbl_sts2 <- tableGrob(E231_stats_sumary, theme = ttheme_minimal())
  grid.newpage()
  grid.draw(E231_tbl_sts2)
#E333
E333_stats_sumary <- rbind("Sample ID" = E333_genes,
                           "Cells_before_filter" = dim(E333_qc)[2],
                           "Cells_after_filter" = dim(E333_gene_object_doublet)[2],
                           "Median Feature per Cell before filter" = median(E333_qc$nFeature_RNA),
                           "Median Reads per Gene before filter" = median(E333_qc$nCount_RNA),
                           "Median Feature per Cell" = median(E333_gene_object_doublet$nFeature_RNA),
                           "Median Reads per Gene" = median(E333_gene_object_doublet$nCount_RNA),
                           #"Max Features" = max.features,
                           "Min Features" = min.features,
                           "Min Counts" = min.counts,
                           #"Max Counts" = max.counts,
                           "MT Percentage" = MT,
                           "NPCs" = 15,
                           "Median Percent MT before Filter" = median(E333_qc@meta.data[["percent.mt"]]),
                           "Median Percent MT after Filter" = median(E333_gene_object_doublet@meta.data[["percent.mt"]]))
  E333_tbl_sts2 <- tableGrob(E333_stats_sumary, theme = ttheme_minimal())
  grid.newpage()
  grid.draw(E333_tbl_sts2)
#E435
E435_stats_sumary <- rbind("Sample ID" = E435_genes,
                           "Cells_before_filter" = dim(E435_qc)[2],
                           "Cells_after_filter" = dim(E435_gene_object_doublet)[2],
                           "Median Feature per Cell before filter" = median(E435_qc$nFeature_RNA),
                           "Median Reads per Gene before filter" = median(E435_qc$nCount_RNA),
                           "Median Feature per Cell" = median(E435_gene_object_doublet$nFeature_RNA),
                           "Median Reads per Gene" = median(E435_gene_object_doublet$nCount_RNA),
                           #"Max Features" = max.features,
                           "Min Features" = min.features,
                           "Min Counts" = min.counts,
                           #"Max Counts" = max.counts,
                           "MT Percentage" = MT,
                           "NPCs" = 15,
                           "Median Percent MT before Filter" = median(E435_qc@meta.data[["percent.mt"]]),
                           "Median Percent MT after Filter" = median(E435_gene_object_doublet@meta.data[["percent.mt"]]))
  E435_tbl_sts2 <- tableGrob(E435_stats_sumary, theme = ttheme_minimal())
  grid.newpage()
  grid.draw(E435_tbl_sts2)


#Save objects with doublets
saveRDS(E170_gene_object_doublets, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_with_doublets_umap_object.rds")
saveRDS(E191_gene_object_doublets, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_with_doublets_umap_object.rds")
saveRDS(E226_gene_object_doublets, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_with_doublets_umap_object.rds")
saveRDS(E231_gene_object_doublets, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_with_doublets_umap_object.rds")
saveRDS(E333_gene_object_doublets, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_with_doublets_umap_object.rds")
saveRDS(E435_gene_object_doublets, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_with_doublets_umap_object.rds")

#Save tables
write.table(E170_stats_sumary, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_stats.csv")
write.table(E191_stats_sumary, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_stats.csv")
write.table(E226_stats_sumary, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_stats.csv")
write.table(E231_stats_sumary, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_stats.csv")
write.table(E333_stats_sumary, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_stats.csv")
write.table(E435_stats_sumary, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_stats.csv")

#Save filtered objects
saveRDS(E170_genes_filtered, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_genes_filtered.rds")
saveRDS(E191_genes_filtered, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_genes_filtered.rds")
saveRDS(E226_genes_filtered, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_genes_filtered.rds")
saveRDS(E231_genes_filtered, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_genes_filtered.rds")
saveRDS(E333_genes_filtered, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_genes_filtered.rds")
saveRDS(E435_genes_filtered, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_genes_filtered.rds")
