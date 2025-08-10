###sscomp Rscript

#Load required libraries
library(sccomp)
library(Seurat)
library(dplyr)
library(ggplot2)

#Load seurat object RDS (RNA + iso assays of all 6 samples)
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Run sccomp
sccomp_result = 
  all_samples_integrated |>
  sccomp_estimate(formula_composition = ~ fertility, 
                  sample = "endo.ID", 
                  cell_group = "cell_type", 
                  cores = 1,
                  verbose = FALSE) |> 
  sccomp_test()

#Save to CSV
write.csv(sccomp_result, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/sccomp_differential_composition_results.csv", row.names = FALSE)

#Box plot of cell subtypes
sccomp_result |> 
  sccomp_boxplot(factor = "fertility")
#Plot of estimates of differential composition and variability
sccomp_result |> 
  plot_1D_intervals()
#Plot relationship between abundance and variability
sccomp_result |> 
  plot_2D_intervals()

#Look for significant differences of inf vs fer
significant <- sccomp_result %>%
  filter(parameter == "fertilityInfertile", c_FDR < 0.05)


sccomp_result |> 
  sccomp_boxplot(factor = "fertility") +
  theme(
    axis.text = element_text(size = 12),        # Increase axis tick labels font size
    axis.title = element_text(size = 12),       # Increase axis titles font size
    legend.text = element_text(size = 12),      # Legend text size
    legend.title = element_text(size = 12),     # Legend title size
    plot.title = element_text(size = 12, face = "bold")  # Plot title font size
  )