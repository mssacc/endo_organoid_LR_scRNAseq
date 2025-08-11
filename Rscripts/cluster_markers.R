#Find all marker genes for each cluster - use for pathway analysis

#Load required libraries
library(Seurat)

#Load RDS
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Write out markers per individual cluster
pre_unciliated_markers <- FindMarkers(all_samples_integrated, assay="RNA", ident.1 = "Pre-Unciliated")
unciliated_markers <- FindMarkers(all_samples_integrated, assay="RNA", ident.1 = "Unciliated")
ciliated_markers <- FindMarkers(all_samples_integrated, assay="RNA", ident.1 = "Ciliated")
secretory_markers <- FindMarkers(all_samples_integrated, assay="RNA", ident.1 = "Secretory")
pre_ciliated_markers <- FindMarkers(all_samples_integrated, assay="RNA", ident.1 = "Pre-Ciliated")
proliferative_markers <- FindMarkers(all_samples_integrated, assay="RNA", ident.1 = "Proliferative")

pre_unciliated_markers <- subset(pre_unciliated_markers, subset= p_val_adj<0.05)
unciliated_markers <- subset(unciliated_markers, subset= p_val_adj<0.05)
ciliated_markers <- subset(ciliated_markers, subset= p_val_adj<0.05)
secretory_markers <- subset(secretory_markers, subset= p_val_adj<0.05)
pre_ciliated_markers <- subset(pre_ciliated_markers, subset= p_val_adj<0.05)
proliferative_markers <- subset(proliferative_markers, subset= p_val_adj<0.05)

#Generate excels
write.csv(pre_unciliated_markers, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/pre_unciliated_markers.csv")
write.csv(unciliated_markers, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/unciliated_markers.csv")
write.csv(ciliated_markers, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/ciliated_markers.csv")
write.csv(secretory_markers, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/secretory_markers.csv")
write.csv(pre_ciliated_markers, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/pre_ciliated_markers.csv")
write.csv(proliferative_markers, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/proliferative_markers.csv")
