#Rscript - Gene ID (ENSIGID) to Symbol

#Set library
.libPaths("/home/mssacc/R_libs_4.4")

#Load required libraries
library(data.table)
library(tidyr)
library(dplyr)
		
#Import gene count matrix
E170_gene_matrix <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/E170_gene_count.csv", header=T, row.names=1)
E191_gene_matrix <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/E191_gene_count.csv", header=T, row.names=1)
E226_gene_matrix <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/E226_gene_count.csv", header=T, row.names=1)
E231_gene_matrix <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/E231_gene_count.csv", header=T, row.names=1)
E333_gene_matrix <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/E333_gene_count.csv", header=T, row.names=1)
E435_gene_matrix <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_output/E435_gene_count.csv", header=T, row.names=1)

#Import empty droplet matrix
E170_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_droplet_analysis/output/E170_gene_count.csv", header=T, row.names = 1)
E191_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_droplet_analysis/output/E191_gene_count.csv", header=T, row.names = 1)
E226_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_droplet_analysis/output/E226_gene_count.csv", header=T, row.names = 1)
E231_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_droplet_analysis/output/E231_gene_count.csv", header=T, row.names = 1)
E333_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_droplet_analysis/output/E333_gene_count.csv", header=T, row.names = 1)
E435_emptydrops <- read.csv("/data/gpfs/projects/punim1901/flames_v2/flames_droplet_analysis/output/E435_gene_count.csv", header=T, row.names = 1)


#Import gene reference list
gene_info <- read.csv("/home/mssacc/genomes/gene_info.csv")
#Remove duplicates from reference list
single_genes <- unique(gene_info)
lookup_table <- setNames(single_genes$gene_symbol, single_genes$gene_id)
		

#Gene count renaming with symbols
#rownames(gene_symbols_matrix) <- lookup_table[rownames(gene_symbols_matrix)]
E170_gene_matrix$gene_symbol <- lookup_table[rownames(E170_gene_matrix)]
E170_gene_matrix <- E170_gene_matrix |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E170_gene_symbols_matrix)[colnames(E170_gene_matrix) != "gene_symbol"],
  sum)

E191_gene_matrix$gene_symbol <- lookup_table[rownames(E191_gene_matrix)]
E191_gene_matrix <- E191_gene_matrix |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E191_gene_matrix)[colnames(E191_gene_matrix) != "gene_symbol"],
  sum)

E226_gene_matrix$gene_symbol <- lookup_table[rownames(E226_gene_matrix)]
E226_gene_matrix <- E226_gene_symbols_matrix |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E226_gene_matrix)[colnames(E226_gene_matrix) != "gene_symbol"],
  sum)

E231_gene_symbols_matrix$gene_symbol <- lookup_table[rownames(E231_gene_symbols_matrix)]
E231_gene_symbols_matrix <- E231_gene_symbols_matrix |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E231_gene_symbols_matrix)[colnames(E231_gene_symbols_matrix) != "gene_symbol"],
  sum)

E333_gene_matrix$gene_symbol <- lookup_table[rownames(E333_gene_matrix)]
E333_gene_matrix <- E333_gene_matrix |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E333_gene_matrix)[colnames(E333_gene_matrix) != "gene_symbol"],
  sum)

E435_gene_matrix$gene_symbol <- lookup_table[rownames(E435_gene_matrix)]
E435_gene_matrix <- E435_gene_matrix |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E435_gene_matrix)[colnames(E435_gene_matrix) != "gene_symbol"],
  sum)	


#Empty droplet renaming with symbols
#rownames(gene_symbols_matrix) <- lookup_table[rownames(gene_symbols_matrix)]
E170_emptydrops$gene_symbol <- lookup_table[rownames(E170_emptydrops)]
E170_emptydrops <- E170_emptydrops |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E170_emptydrops)[colnames(E170_emptydrops) != "gene_symbol"],
               sum)

E191_emptydrops$gene_symbol <- lookup_table[rownames(E191_emptydrops)]
E191_emptydrops <- E191_emptydrops |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E191_emptydrops)[colnames(E191_emptydrops) != "gene_symbol"],
               sum)

E226_emptydrops$gene_symbol <- lookup_table[rownames(E226_emptydrops)]
E226_emptydrops <- E226_emptydrops |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E226_emptydrops)[colnames(E226_emptydrops) != "gene_symbol"],
               sum)

E231_emptydrops$gene_symbol <- lookup_table[rownames(E231_emptydrops)]
E231_emptydrops <- E231_emptydrops |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E231_emptydrops)[colnames(E231_emptydrops) != "gene_symbol"],
               sum)

E333_emptydrops$gene_symbol <- lookup_table[rownames(E333_emptydrops)]
E333_emptydrops <- E333_emptydrops |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E333_emptydrops)[colnames(E333_emptydrops) != "gene_symbol"],
               sum)

E435_emptydrops$gene_symbol <- lookup_table[rownames(E435_emptydrops)]
E435_emptydrops <- E435_emptydrops |> 
  group_by(gene_symbol) |>
  summarise_at(colnames(E435_emptydrops)[colnames(E435_emptydrops) != "gene_symbol"],
               sum)	


#Write output gene symbol files to a csv file
write.csv(E170_gene_symbols_matrix, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_gene_count_symbols.csv")
write.csv(E191_gene_symbols_matrix, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_gene_count_symbols.csv")
write.csv(E226_gene_symbols_matrix, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_gene_count_symbols.csv")
write.csv(E231_gene_symbols_matrix, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_gene_count_symbols.csv")
write.csv(E333_gene_symbols_matrix, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_gene_count_symbols.csv")
write.csv(E435_gene_symbols_matrix, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_gene_count_symbols.csv")

#Write output empty droplet symbol files to a csv file
write.csv(E170_emptydrops, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E170_empty_droplets.csv")
write.csv(E191_emptydrops, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E191_empty_droplets.csv")
write.csv(E226_emptydrops, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E226_empty_droplets.csv")
write.csv(E231_emptydrops, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E231_empty_droplets.csv")
write.csv(E333_emptydrops, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E333_empty_droplets.csv")
write.csv(E435_emptydrops, row.names=FALSE, "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/E435_empty_droplets.csv")
