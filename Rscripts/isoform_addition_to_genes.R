#Adding isoform counts to gene count seurat object

#Load required libraries
library(Matrix)
library(dplyr)
library(SeuratObject)
library(Seurat)
library(patchwork)
library(ggplot2)

#Read in the resource table (transcript_id, gene_id, gene_symbol)
reference <- read.csv("/data/gpfs/projects/punim1901/flames_v2/naming_reference.csv")

#Import the transcript data
#E170
  E170_counts <- readMM("/data/gpfs/projects/punim1901/flames_v2/flames_output/E170_matched_reads_dedup.count.mtx")
  E170_barcodes <- readLines("/data/gpfs/projects/punim1901/flames_v2/flames_output/E170_matched_reads_dedup.barcodes.txt")
  E170_features <- read.delim("/data/gpfs/projects/punim1901/flames_v2/flames_output/E170_matched_reads_dedup.features.txt", header = FALSE)
#E191
  E191_counts <- readMM("/data/gpfs/projects/punim1901/flames_v2/flames_output/E191_matched_reads_dedup.count.mtx")
  E191_barcodes <- readLines("/data/gpfs/projects/punim1901/flames_v2/flames_output/E191_matched_reads_dedup.barcodes.txt")
  E191_features <- read.delim("/data/gpfs/projects/punim1901/flames_v2/flames_output/E191_matched_reads_dedup.features.txt", header = FALSE)
#E226
  E226_counts <- readMM("/data/gpfs/projects/punim1901/flames_v2/flames_output/E226_matched_reads_dedup.count.mtx")
  E226_barcodes <- readLines("/data/gpfs/projects/punim1901/flames_v2/flames_output/E226_matched_reads_dedup.barcodes.txt")
  E226_features <- read.delim("/data/gpfs/projects/punim1901/flames_v2/flames_output/E226_matched_reads_dedup.features.txt", header = FALSE)
#E231
  E231_counts <- readMM("/data/gpfs/projects/punim1901/flames_v2/flames_output/E231_matched_reads_dedup.count.mtx")
  E231_barcodes <- readLines("/data/gpfs/projects/punim1901/flames_v2/flames_output/E231_matched_reads_dedup.barcodes.txt")
  E231_features <- read.delim("/data/gpfs/projects/punim1901/flames_v2/flames_output/E231_matched_reads_dedup.features.txt", header = FALSE)
#E333
  E333_counts <- readMM("/data/gpfs/projects/punim1901/flames_v2/flames_output/E333_matched_reads_dedup.count.mtx")
  E333_barcodes <- readLines("/data/gpfs/projects/punim1901/flames_v2/flames_output/E333_matched_reads_dedup.barcodes.txt")
  E333_features <- read.delim("/data/gpfs/projects/punim1901/flames_v2/flames_output/E333_matched_reads_dedup.features.txt", header = FALSE)
#E435
  E435_counts <- readMM("/data/gpfs/projects/punim1901/flames_v2/flames_output/E435_matched_reads_dedup.count.mtx")
  E435_barcodes <- readLines("/data/gpfs/projects/punim1901/flames_v2/flames_output/E435_matched_reads_dedup.barcodes.txt")
  E435_features <- read.delim("/data/gpfs/projects/punim1901/flames_v2/flames_output/E435_matched_reads_dedup.features.txt", header = FALSE)

#Transpose the matrix for Seurat compatibility
  E170_counts <- t(E170_counts)
  E191_counts <- t(E191_counts)
  E226_counts <- t(E226_counts)
  E231_counts <- t(E231_counts)
  E333_counts <- t(E333_counts)
  E435_counts <- t(E435_counts)

#Set row and column names
  rownames(E170_counts) <- E170_features$V1
    colnames(E170_counts) <- E170_barcodes
  rownames(E191_counts) <- E191_features$V1
    colnames(E191_counts) <- E191_barcodes
  rownames(E226_counts) <- E226_features$V1
    colnames(E226_counts) <- E226_barcodes
  rownames(E231_counts) <- E231_features$V1
    colnames(E231_counts) <- E231_barcodes
  rownames(E333_counts) <- E333_features$V1
    colnames(E333_counts) <- E333_barcodes
  rownames(E435_counts) <- E435_features$V1
    colnames(E435_counts) <- E435_barcodes

#Convert to a data frame
  E170_counts_df <- as.data.frame(as.matrix(E170_counts))
    E170_counts_df <- as.data.frame(E170_counts_df)
  E191_counts_df <- as.data.frame(as.matrix(E191_counts))
    E191_counts_df <- as.data.frame(E191_counts_df)
  E226_counts_df <- as.data.frame(as.matrix(E226_counts))
    E226_counts_df <- as.data.frame(E226_counts_df)
  E231_counts_df <- as.data.frame(as.matrix(E231_counts))
    E231_counts_df <- as.data.frame(E231_counts_df)
  E333_counts_df <- as.data.frame(as.matrix(E333_counts))
    E333_counts_df <- as.data.frame(E333_counts_df)
  E435_counts_df <- as.data.frame(as.matrix(E435_counts))
    E435_counts_df <- as.data.frame(E435_counts_df)

#Add transcript_id as the first column
  E170_counts_df$transcript_id <- rownames(E170_counts_df)
    E170_counts_df <- E170_counts_df[, c(ncol(E170_counts_df), 1:(ncol(E170_counts_df)-1))]
  E191_counts_df$transcript_id <- rownames(E191_counts_df)
    E191_counts_df <- E191_counts_df[, c(ncol(E191_counts_df), 1:(ncol(E191_counts_df)-1))]
  E226_counts_df$transcript_id <- rownames(E226_counts_df)
    E226_counts_df <- E226_counts_df[, c(ncol(E226_counts_df), 1:(ncol(E226_counts_df)-1))]
  E231_counts_df$transcript_id <- rownames(E231_counts_df)
    E231_counts_df <- E231_counts_df[, c(ncol(E231_counts_df), 1:(ncol(E231_counts_df)-1))]
  E333_counts_df$transcript_id <- rownames(E333_counts_df)
    E333_counts_df <- E333_counts_df[, c(ncol(E333_counts_df), 1:(ncol(E333_counts_df)-1))]
  E435_counts_df$transcript_id <- rownames(E435_counts_df)
    E435_counts_df <- E435_counts_df[, c(ncol(E435_counts_df), 1:(ncol(E435_counts_df)-1))]

#Merge with the resource table to add gene symbols
  E170_df_genesymbol <- E170_counts_df %>%
    left_join(reference, by = "transcript_id")
  E191_df_genesymbol <- E191_counts_df %>%
    left_join(reference, by = "transcript_id")
  E226_df_genesymbol <- E226_counts_df %>%
    left_join(reference, by = "transcript_id")
  E231_df_genesymbol <- E231_counts_df %>%
    left_join(reference, by = "transcript_id")
  E333_df_genesymbol <- E333_counts_df %>%
    left_join(reference, by = "transcript_id")
  E435_df_genesymbol <- E435_counts_df %>%
    left_join(reference, by = "transcript_id")

#Remove the gene_id column and reorder the columns
  E170_df_genesymbol$gene_id <- NULL
    E170_df_genesymbol <- E170_df_genesymbol[, c(ncol(E170_df_genesymbol), 1:(ncol(E170_df_genesymbol)-1))]
  E191_df_genesymbol$gene_id <- NULL
    E191_df_genesymbol <- E191_df_genesymbol[, c(ncol(E191_df_genesymbol), 1:(ncol(E191_df_genesymbol)-1))]
  E226_df_genesymbol$gene_id <- NULL
    E226_df_genesymbol <- E226_df_genesymbol[, c(ncol(E226_df_genesymbol), 1:(ncol(E226_df_genesymbol)-1))]
  E231_df_genesymbol$gene_id <- NULL
    E231_df_genesymbol <- E231_df_genesymbol[, c(ncol(E231_df_genesymbol), 1:(ncol(E231_df_genesymbol)-1))]
  E333_df_genesymbol$gene_id <- NULL
    E333_df_genesymbol <- E333_df_genesymbol[, c(ncol(E333_df_genesymbol), 1:(ncol(E333_df_genesymbol)-1))]
  E435_df_genesymbol$gene_id <- NULL
    E435_df_genesymbol <- E435_df_genesymbol[, c(ncol(E435_df_genesymbol), 1:(ncol(E435_df_genesymbol)-1))]

#Update row names to include gene symbol instead of transcript_id
  rownames(E170_df_genesymbol) <- paste0(E170_df_genesymbol$transcript_id, "_", E170_df_genesymbol$gene_symbol)
    E170_df_genesymbol$transcript_id <- NULL
    E170_df_genesymbol$gene_symbol <- NULL
  rownames(E191_df_genesymbol) <- paste0(E191_df_genesymbol$transcript_id, "_", E191_df_genesymbol$gene_symbol)
    E191_df_genesymbol$transcript_id <- NULL
    E191_df_genesymbol$gene_symbol <- NULL
  rownames(E226_df_genesymbol) <- paste0(E226_df_genesymbol$transcript_id, "_", E226_df_genesymbol$gene_symbol)
    E226_df_genesymbol$transcript_id <- NULL
    E226_df_genesymbol$gene_symbol <- NULL
  rownames(E231_df_genesymbol) <- paste0(E231_df_genesymbol$transcript_id, "_", E231_df_genesymbol$gene_symbol)
    E231_df_genesymbol$transcript_id <- NULL
    E231_df_genesymbol$gene_symbol <- NULL
  rownames(E333_df_genesymbol) <- paste0(E333_df_genesymbol$transcript_id, "_", E333_df_genesymbol$gene_symbol)
    E333_df_genesymbol$transcript_id <- NULL
    E333_df_genesymbol$gene_symbol <- NULL
  rownames(E435_df_genesymbol) <- paste0(E435_df_genesymbol$transcript_id, "_", E435_df_genesymbol$gene_symbol)
    E435_df_genesymbol$transcript_id <- NULL
    E435_df_genesymbol$gene_symbol <- NULL

#Write the output to a CSV file
  write.csv(E170_df_genesymbol, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_iso.csv")
  write.csv(E191_df_genesymbol, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_iso.csv")
  write.csv(E226_df_genesymbol, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_iso.csv")
  write.csv(E231_df_genesymbol, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_iso.csv")
  write.csv(E333_df_genesymbol, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_iso.csv")
  write.csv(E435_df_genesymbol, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_iso.csv")

E170_iso <- E170_df_genesymbol
E191_iso <- E191_df_genesymbol
E226_iso <- E226_df_genesymbol
E231_iso <- E231_df_genesymbol
E333_iso <- E333_df_genesymbol
E435_iso <- E435_df_genesymbol


#Import gene matrices
E170 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_genes_filtered.rds")
E191 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_genes_filtered.rds")
E226 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_genes_filtered.rds")
E231 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_genes_filtered.rds")
E333 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_genes_filtered.rds")
E435 <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_genes_filtered.rds")

E170_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_iso.csv", header=TRUE, row.names=1)
E191_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_iso.csv", header=TRUE, row.names=1)
E226_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_iso.csv", header=TRUE, row.names=1)
E231_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_iso.csv", header=TRUE, row.names=1)
E333_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_iso.csv", header=TRUE, row.names=1)
E435_iso <- read.csv("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_iso.csv", header=TRUE, row.names=1)


###Create Seurat objects with isoform expression data
  E170_iso_object <- CreateSeuratObject(counts=E170_iso, project="E170")
  E191_iso_object <- CreateSeuratObject(counts=E191_iso, project="E191")
  E226_iso_object <- CreateSeuratObject(counts=E226_iso, project="E226")
  E231_iso_object <- CreateSeuratObject(counts=E231_iso, project="E231")
  E333_iso_object <- CreateSeuratObject(counts=E333_iso, project="E333")
  E435_iso_object <- CreateSeuratObject(counts=E435_iso, project="E435")

VlnPlot(E170_iso_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + plot_annotation(title = "E170 QC plots (isoform level) BEFORE Filtering")
VlnPlot(E191_iso_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + plot_annotation(title = "E191 QC plots (isoform level) BEFORE Filtering")
VlnPlot(E226_iso_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + plot_annotation(title = "E226 QC plots (isoform level) BEFORE Filtering")
VlnPlot(E231_iso_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + plot_annotation(title = "E231 QC plots (isoform level) BEFORE Filtering")
VlnPlot(E333_iso_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + plot_annotation(title = "E333 QC plots (isoform level) BEFORE Filtering")
VlnPlot(E435_iso_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) + plot_annotation(title = "E435 QC plots (isoform level) BEFORE Filtering")

###Filter the new Seurat objects based on gene level information
#Filter the data so isoform and gene cells match
  E170_oarfish_iso_matched_gene <- subset(E170_iso_object, cells =E170@graphs[["RNA_nn"]]@Dimnames[[1]])
  E191_oarfish_iso_matched_gene <- subset(E191_iso_object, cells =E191@graphs[["RNA_nn"]]@Dimnames[[1]])
  E226_oarfish_iso_matched_gene <- subset(E226_iso_object, cells =E226@graphs[["RNA_nn"]]@Dimnames[[1]])
  E231_oarfish_iso_matched_gene <- subset(E231_iso_object, cells =E231@graphs[["RNA_nn"]]@Dimnames[[1]])
  E333_oarfish_iso_matched_gene <- subset(E333_iso_object, cells =E333@graphs[["RNA_nn"]]@Dimnames[[1]])
  E435_oarfish_iso_matched_gene <- subset(E435_iso_object, cells =E435@graphs[["RNA_nn"]]@Dimnames[[1]])

#E170
  #Rejoin data sets after integration
    E170_oarfish_iso_matched_gene <- JoinLayers(E170_oarfish_iso_matched_gene)
    E170_counts_table_iso <- E170_oarfish_iso_matched_gene[["RNA"]]$counts
    as.data.frame(E170_counts_table_iso) -> E170_df_iso
  #Remove rows where the sum is 0
    E170_df_iso <- E170_df_iso[rowSums(E170_df_iso) != 0, ]
  #Create iso assay
    E170[["iso"]] <- CreateAssay5Object(counts = E170_df_iso)

#E191
  #Rejoin data sets after integration
    E191_oarfish_iso_matched_gene <- JoinLayers(E191_oarfish_iso_matched_gene)
    E191_counts_table_iso <- E191_oarfish_iso_matched_gene[["RNA"]]$counts
    as.data.frame(E191_counts_table_iso) -> E191_df_iso
  #Remove rows where the sum is 0
    E191_df_iso <- E191_df_iso[rowSums(E191_df_iso) != 0, ]
  #Create iso assay
    E191[["iso"]] <- CreateAssay5Object(counts = E191_df_iso)

#E226
  #Rejoin data sets after integration
    E226_oarfish_iso_matched_gene <- JoinLayers(E226_oarfish_iso_matched_gene)
    E226_counts_table_iso <- E226_oarfish_iso_matched_gene[["RNA"]]$counts
    as.data.frame(E226_counts_table_iso) -> E226_df_iso
  #Remove rows where the sum is 0
    E226_df_iso <- E226_df_iso[rowSums(E226_df_iso) != 0, ]
  #Create iso assay
    E226[["iso"]] <- CreateAssay5Object(counts = E226_df_iso)

#E231
  #Rejoin data sets after integration
    E231_oarfish_iso_matched_gene <- JoinLayers(E231_oarfish_iso_matched_gene)
    E231_counts_table_iso <- E231_oarfish_iso_matched_gene[["RNA"]]$counts
    as.data.frame(E231_counts_table_iso) -> E231_df_iso
  #Remove rows where the sum is 0
    E231_df_iso <- E231_df_iso[rowSums(E231_df_iso) != 0, ]
  #Create iso assay
    E231[["iso"]] <- CreateAssay5Object(counts = E231_df_iso)

#E333
  #Rejoin data sets after integration
    E333_oarfish_iso_matched_gene <- JoinLayers(E333_oarfish_iso_matched_gene)
    E333_counts_table_iso <- E333_oarfish_iso_matched_gene[["RNA"]]$counts
    as.data.frame(E333_counts_table_iso) -> E333_df_iso
  #Remove rows where the sum is 0
    E333_df_iso <- E333_df_iso[rowSums(E333_df_iso) != 0, ]
  #Create iso assay
    E333[["iso"]] <- CreateAssay5Object(counts = E333_df_iso)

#E435
  #Rejoin data sets after integration
    E435_oarfish_iso_matched_gene <- JoinLayers(E435_oarfish_iso_matched_gene)
    E435_counts_table_iso <- E435_oarfish_iso_matched_gene[["RNA"]]$counts
    as.data.frame(E435_counts_table_iso) -> E435_df_iso
  #Remove rows where the sum is 0
    E435_df_iso <- E435_df_iso[rowSums(E435_df_iso) != 0, ]
  #Create iso assay
    E435[["iso"]] <- CreateAssay5Object(counts = E435_df_iso)


#Normalize the new assay data
#E170
    E170 <- NormalizeData(E170, assay = "iso")
    E170 <- FindVariableFeatures(E170, assay = "iso")
    E170 <- ScaleData(E170, assay = "iso")
    E170 <- RunPCA(E170, assay = "iso", reduction.name = "pca_iso")
  #Run UMAP
    E170 <- RunUMAP(E170, reduction = "pca_iso", dims = 1:15, assay = "iso", reduction.name = "umap_iso")
  #check to see that we have two assays
    E170

#E191
    E191 <- NormalizeData(E191, assay = "iso")
    E191 <- FindVariableFeatures(E191, assay = "iso")
    E191 <- ScaleData(E191, assay = "iso")
    E191 <- RunPCA(E191, assay = "iso", reduction.name = "pca_iso")
  #Run UMAP
    E191 <- RunUMAP(E191, reduction = "pca_iso", dims = 1:15, assay = "iso", reduction.name = "umap_iso")
  #Check to see that we have two assays
    E191

#E226
    E226 <- NormalizeData(E226, assay = "iso")
    E226 <- FindVariableFeatures(E226, assay = "iso")
    E226 <- ScaleData(E226, assay = "iso")
    E226 <- RunPCA(E226, assay = "iso", reduction.name = "pca_iso")
  #Run UMAP
    E226 <- RunUMAP(E226, reduction = "pca_iso", dims = 1:15, assay = "iso", reduction.name = "umap_iso")
  #Check to see that we have two assays
    E226

#E231
    E231 <- NormalizeData(E231, assay = "iso")
    E231 <- FindVariableFeatures(E231, assay = "iso")
    E231 <- ScaleData(E231, assay = "iso")
    E231 <- RunPCA(E231, assay = "iso", reduction.name = "pca_iso")
  #Run UMAP
    E231 <- RunUMAP(E231, reduction = "pca_iso", dims = 1:15, assay = "iso", reduction.name = "umap_iso")
  #Check to see that we have two assays
    E231

#E333
    E333 <- NormalizeData(E333, assay = "iso")
    E333 <- FindVariableFeatures(E333, assay = "iso")
    E333 <- ScaleData(E333, assay = "iso")
    E333 <- RunPCA(E333, assay = "iso", reduction.name = "pca_iso")
  #Run UMAP
    E333 <- RunUMAP(E333, reduction = "pca_iso", dims = 1:15, assay = "iso", reduction.name = "umap_iso")
  #Check to see that we have two assays
    E333

#E435
    E435 <- NormalizeData(E435, assay = "iso")
    E435 <- FindVariableFeatures(E435, assay = "iso")
    E435 <- ScaleData(E435, assay = "iso")
    E435 <- RunPCA(E435, assay = "iso", reduction.name = "pca_iso")
  #Run UMAP
    E435 <- RunUMAP(E435, reduction = "pca_iso", dims = 1:15, assay = "iso", reduction.name = "umap_iso")
  #Check to see that we have two assays
    E435

#Visualize the UMAP for gene and isoform
  DimPlot(E170, label = TRUE, reduction = "umap") + ggtitle("E170 UMAP gene level clustering") |
    DimPlot(E170, label = TRUE, reduction = "umap_iso") + ggtitle("E170 UMAP isoform level clustering")
  DimPlot(E191, label = TRUE, reduction = "umap") + ggtitle("E191 UMAP gene level clustering") |
    DimPlot(E191, label = TRUE, reduction = "umap_iso") + ggtitle("E191 UMAP isoform level clustering")
  DimPlot(E226, label = TRUE, reduction = "umap") + ggtitle("E226 UMAP gene level clustering") |
    DimPlot(E226, label = TRUE, reduction = "umap_iso") + ggtitle("E226 UMAP isoform level clustering")
  DimPlot(E231, label = TRUE, reduction = "umap") + ggtitle("E231 UMAP gene level clustering") |
    DimPlot(E231, label = TRUE, reduction = "umap_iso") + ggtitle("E231 UMAP isoform level clustering")
  DimPlot(E333, label = TRUE, reduction = "umap") + ggtitle("E333 UMAP gene level clustering") |
    DimPlot(E333, label = TRUE, reduction = "umap_iso") + ggtitle("E333 UMAP isoform level clustering")
  DimPlot(E435, label = TRUE, reduction = "umap") + ggtitle("E435 UMAP gene level clustering") |
    DimPlot(E435, label = TRUE, reduction = "umap_iso") + ggtitle("E435 UMAP isoform level clustering")
  
#Save RDS
  saveRDS(E170, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_gene_and_isoform_seurat.rds")
  saveRDS(E191, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_gene_and_isoform_seurat.rds")
  saveRDS(E226, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_gene_and_isoform_seurat.rds")
  saveRDS(E231, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_gene_and_isoform_seurat.rds")
  saveRDS(E333, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_gene_and_isoform_seurat.rds")
  saveRDS(E435, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_gene_and_isoform_seurat.rds")

#Save workspace
save.image("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/gene_and_isoform_overlay.RData")
