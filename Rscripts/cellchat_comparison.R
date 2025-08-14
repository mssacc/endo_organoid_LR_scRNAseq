#Rscript - Comparative analysis of Fertile vs Infertile datasets using CellChat

#Load the required libraries
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(presto)
library(wordcloud)
ptm = Sys.time()


#Load CellChat object of each dataset and merge them together
cellchat.fer <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat/cellchat_fertile.rds")
cellchat.inf <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat/cellchat_infertile.rds")

#Merge
object.list <- list(Fertile = cellchat.fer, Infertile = cellchat.inf)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

execution.time = Sys.time() - ptm

#Define custom colours
colours.celltype <- c("Pre-Unciliated" = "#F8766D", 
                      "Unciliated"     = "#ABA300", 
                      "Ciliated"       = "#0CB702", 
                      "Secretory"      = "#00BFC4", 
                      "Pre-Ciliated"   = "#849AFF",
                      "Proliferative"  = "#FF61CC")
colours.fertility <- c("Fertile"   = "palegreen3",
                       "Infertile" = "#F8766D")


#Save merged cell chat object
save(object.list, file = "cellchat_object.fer_inf.RData")
save(cellchat, file = "cellchat_merged_fer_inf.RData")



###Part I: Identify altered interactions and cell populations
#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Compare the number of interactions and interaction strength among different cell populations
#(A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, color.use = colours.celltype, color.edge = c("#F8766D", "palegreen3"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = colours.celltype, color.edge = c("#F8766D", "palegreen3"))

#(B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat, color.use = colours.celltype)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = colours.celltype)
gg1 + gg2

#(C) Circle plot showing the number of interactions or interaction strength among different cell populations across multiple datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, color.use = colours.celltype, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction Weights - ", names(object.list)[i]))
}

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, color.use = colours.celltype, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of Interactions - ", names(object.list)[i]))
}

#Compare the major sources and targets in a 2D space
#(A) Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], color.use = colours.celltype, weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


#(B) Identify the signaling changes of specific cell populations
my_colors <- c("Shared" = "black",
               "Fertile specific" = "palegreen3",
               "Infertile specific" = "#F8766D")

gg0 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Pre-Unciliated") +
          scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Unciliated") +
          scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Ciliated") +
          scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Secretory") +
          scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Pre-Ciliated") +
          scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Proliferative") +
          scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)

patchwork::wrap_plots(plots = list(gg0))
patchwork::wrap_plots(plots = list(gg1))
patchwork::wrap_plots(plots = list(gg2))
patchwork::wrap_plots(plots = list(gg3))
patchwork::wrap_plots(plots = list(gg4))
patchwork::wrap_plots(plots = list(gg5))

combined_plot <- (gg0 + gg1 + gg2) / (gg3 + gg4 + gg5) + plot_layout(guides = "collect") &
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))
combined_plot


#Part II: Identify altered signaling with distinct network architecture and interaction strength
#Identify altered signaling with distinct interaction strength
#(A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
#CellChat can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network)
#When setting do.stat = TRUE, a paired Wilcoxon test is performed to determine whether there is a significant difference of the signaling information flow between two conditions. The top signaling pathways colored red are enriched in dataset 1, and these colored greens are enriched in dataset 2.
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, color.use = colours.fertility)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, color.use = colours.fertility)
gg1 + gg2


#Part III: Identify the up-regulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabilities
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Infertile", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Infertile", angle.x = 45, remove.isolate = T)
gg1 + gg2


#Identify dysfunctional signaling by using differential expression analysis
#Define a positive dataset (the dataset with positive fold change against the other dataset)
pos.dataset = "Infertile"

#Define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

#Perform differential expression analysis 
#CellChat v2 performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc. Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

#Map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)

#Extract the ligand-receptor pairs with upregulated ligands in infertile
net.up <- subsetCommunication(cellchat, net = net, datasets = "Infertile", ligand.logFC = 0.05, receptor.logFC = NULL)

#Extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Fertile", ligand.logFC = -0.05, receptor.logFC = NULL)

#Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


#Part V: Compare the signaling gene expression distribution between different datasets
#Plot gene expression distribution of signaling genes related to LR pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Fertile", "Infertile")) #set factor level


#Save the merged CellChat object
save(object.list, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_object_fer_v_inf.RData")
save(cellchat, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_fer_v_inf.RData")
