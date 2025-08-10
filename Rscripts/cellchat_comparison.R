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


# Users can now export the merged CellChat object and the list of the two separate objects for later use
# save(object.list, file = "cellchat_object.fer_inf.RData")
# save(cellchat, file = "cellchat_merged_fer_inf.RData")



###Part I: Identify altered interactions and cell populations
#Compare the total number of interactions and interaction strength
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
#Compare the number of interactions and interaction strength among different cell populations
#A) Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
#The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.  
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, color.use = colours.celltype, color.edge = c("#F8766D", "palegreen3"))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = colours.celltype, color.edge = c("#F8766D", "palegreen3"))


#(B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat, color.use = colours.celltype)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = colours.celltype)
#> Do heatmap based on a merged object
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
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

my_colors <- c("Shared" = "black",
               "Fertile specific" = "palegreen3",     # red - change to your preferred hex
               "Infertile specific" = "#F8766D")   # blue

#(B) Identify the signaling changes of specific cell populations
gg0 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Pre-Unciliated") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Unciliated") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Ciliated") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Secretory") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Pre-Ciliated") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Proliferative") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)

#> Visualizing differential outgoing and incoming signaling changes from fertile to infertile
patchwork::wrap_plots(plots = list(gg0))
patchwork::wrap_plots(plots = list(gg1))
patchwork::wrap_plots(plots = list(gg2))
patchwork::wrap_plots(plots = list(gg3))
patchwork::wrap_plots(plots = list(gg4))
patchwork::wrap_plots(plots = list(gg5))

combined_plot <- (gg0 + gg1 + gg2) / (gg3 + gg4 + gg5) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))
combined_plot

#Part II: Identify altered signaling with distinct network architecture and interaction strength
#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

#Compute and visualize the pathway distance in the learned joint manifold
#CellChat can identify the signaling networks with larger (or smaller) difference based on their Euclidean distance in the shared two-dimensions space. Larger distance implies larger difference of the communication networks between two datasets in terms of either functional or structure similarity. 
rankSimilarity(cellchat, type = "functional")

#Identify altered signaling with distinct interaction strength
#(A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
#CellChat can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network)
#When setting do.stat = TRUE, a paired Wilcoxon test is performed to determine whether there is a significant difference of the signaling information flow between two conditions. The top signaling pathways colored red are enriched in dataset 1, and these colored greens are enriched in dataset 2.
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, color.use = colours.fertility)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, color.use = colours.fertility)
gg1 + gg2

#(B) Compare outgoing (or incoming) signaling patterns associated with each cell population
#Combining all the identified signaling pathways from different datasets
for (i in seq_len(length(object.list) - 1)) {
  # Now 'i + 1' will always be within bounds
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)}
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#Part III: Identify the up-regulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabilities
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 2), angle.x = 45)
#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:6),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Infertile", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Infertile", angle.x = 45, remove.isolate = T)
gg1 + gg2


#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = c(1:6),  comparison = c(1, 2),signaling = pathways.show, max.dataset = 2, title.name = "Increased signaling in Infertile", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = c(1:6),  comparison = c(1, 2), signaling = pathways.show,max.dataset = 1, title.name = "Decreased signaling in Infertile", angle.x = 45, remove.isolate = T)
gg1 + gg2


#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, sources.use = c(1:6), targets.use = c(1:6),   signaling = pathways.show)

#Identify dysfunctional signaling by using differential expression anahttps://spartan-ood.hpc.unimelb.edu.au/rnode/spartan-bm131.hpc.unimelb.edu.au/36639/graphics/plot_zoom_png?width=499&height=831lysis
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

#Visualize the identified up-regulated and down-regulated signaling ligand-receptor pairs
#(A) Bubble plot to visualize up-regulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:6), targets.use = c(1:6), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:6), targets.use = c(1:6), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

#(B) Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 1, targets.use = c(2:3), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1, targets.use = c(2:3), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway


#Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
#CellChat can visually compare cell-cell communication networks using hierarchy plot, circle plot, chord diagram, or heatmap
pathways.show <- c("BMP") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      color.use = colours.celltype)}


#Part V: Compare the signaling gene expression distribution between different datasets
#Plot gene expression distribution of signaling genes related to LR pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Fertile", "Infertile")) # set factor level
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", color.use = colours.fertility, colors.ggplot = T, type = "violin")


#Save the merged CellChat object
save(object.list, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_object_fer_v_inf.RData")
save(cellchat, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_fer_v_inf.RData")
