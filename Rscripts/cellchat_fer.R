#Rscript - CellChat Fertile Samples

#Set library
.libPaths("/home/mssacc/R_libs_4.4")

#Load required libraries
library(presto)
library(ggplot2)
library(BiocNeighbors)
library(CellChat, lib="/home/mssacc/R_libs_4.4")
library(patchwork)
library(SeuratObject)
library(Seurat)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T)

###Part I: Data input & processing and initialization of CellChat object
  #Load data
  all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")
  all_samples_integrated$samples=all_samples_integrated$orig.ident
  all_samples_integrated$samples[all_samples_integrated$samples %in% "E170_genes"] <- "E170"
  all_samples_integrated$samples[all_samples_integrated$samples %in% "E191_genes"] <- "E191"
  all_samples_integrated$samples[all_samples_integrated$samples %in% "E226_genes"] <- "E226"
  all_samples_integrated$samples[all_samples_integrated$samples %in% "E231_genes"] <- "E231"
  all_samples_integrated$samples[all_samples_integrated$samples %in% "E333_genes"] <- "E333"
  all_samples_integrated$samples[all_samples_integrated$samples %in% "E435_genes"] <- "E435"

  #Subset fertile samples
  fertile_samples <- subset(all_samples_integrated, subset = samples %in% c("E191", "E333", "E435"))

  #Normalized data matrix
  data.input <- fertile_samples[["RNA"]]$data
  labels <- Idents(fertile_samples)
  #Create a dataframe of the cell labels
  meta <- data.frame(labels = labels, row.names = names(labels))

  #Create a CellChat object
  cellchat <- createCellChat(object = fertile_samples, group.by = "ident", assay = "RNA")

  #Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)

  #Use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
  CellChatDB.use <- subsetDB(CellChatDB)

  #Set the used database in the object
  cellchat@DB <- CellChatDB.use

  #Preprocessing the expression data for cell-cell communication analysis
  #Subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 8) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  options(future.globals.maxSize = 1000 * 1024^2)  # Set to 1 GB
  cellchat <- identifyOverExpressedInteractions(cellchat)


###Part II: Inference of cell-cell communication network
  #Compute the communication probability and infer cellular communication network
  ptm = Sys.time()
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  #took 1 hour to run with 8 cpus
  #> triMean is used for calculating the average gene expression per cell group. 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #Extract the inferred cellular communication network as a data frame
  #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
  df.net <- subsetCommunication(cellchat)

  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)

  #Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  execution.time = Sys.time() - ptm
  print(as.numeric(execution.time, units = "secs"))
  #CellChat can also visualize the aggregated cell-cell communication network.
  ptm = Sys.time()
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  #Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
  mat <- cellchat@net$weight
  par(mfrow = c(2,3), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }


###Part III: Visualization of cell-cell communication network
  #In all visualization plots, edge colors are consistent with the sources as sender, and edge weights are proportional to the interaction strength. Thicker edge line indicates a stronger signal. In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group. In the hierarchy plot, solid and open circles represent source and target, respectively. In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. The inner bar size is proportional to the signal strength received by the targets. Such inner bar is helpful for interpreting the complex chord diagram. Note that there exist some inner bars without any chord for some cell groups, please just igore it because this is an issue that has not been addressed by circlize package.
  #Here we take input of one signaling pathway as an example.
  pathways.show <- c("WNT") 
  #All the signaling pathways showing significant communications can be accessed by 
  cellchat@netP$pathways
  #Hierarchy plot
  #Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  vertex.receiver = seq(1,4) # a numeric vector. 
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  #Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

  #Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
  netAnalysis_contribution(cellchat, signaling = pathways.show)

  #We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
  pairLR.WNT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  LR.show <- pairLR.WNT[8,] #show one ligand-receptor pair, number indicates pair selected - listed in alphabetical/numerical order from plot generated above
  #Hierarchy plot
  vertex.receiver = seq(1,4) #a numeric vector
  netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
  #Circle plot
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

  #Access all the signaling pathways showing significant communications
  pathways.show.all <- cellchat@netP$pathways
  #Check the order of cell identity to set suitable vertex.receiver
  levels(cellchat@idents)
  vertex.receiver = seq(1,5)

  ##Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
  #Bubble plot
  #Show all the significant interactions (LR pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2,4), remove.isolate = FALSE)
  netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1), remove.isolate = FALSE)

  #Show all the significant interactions (LR pairs) associated with certain signaling pathways
  netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:5), signaling = c("WNT"), remove.isolate = FALSE)
  #Show all the significant interactions (LR pairs) based on user's input (defined by `pairLR.use`)
  pairLR.use <- extractEnrichedLR(cellchat, signaling = c("WNT","NOTCH"))
  netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:5), pairLR.use = pairLR.use, remove.isolate = TRUE)
  
  #Chord diagram
  #Show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(2), lab.cex = 0.5,legend.pos.y = 30)
  netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 5, legend.pos.x = 15)
  #Show all the significant interactions (L-R pairs) associated with certain signaling pathways
  netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5), signaling = c("WNT","NOTCH"),legend.pos.x = 8)
  #Show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

  #Plot the signaling gene expression distribution using violin/dot plot
  plotGeneExpression(cellchat, signaling = "WNT", enriched.only = TRUE, type = "violin")
  print(as.numeric(execution.time, units = "secs"))
  #By default, plotGeneExpression only shows the expression of signaling genes related to the inferred significant communications. USERS can show the expression of all signaling genes related to one signaling pathway by
  plotGeneExpression(cellchat, signaling = "WNT", enriched.only = FALSE)
  execution.time = Sys.time() - ptm


###Part IV: Systems analysis of cell-cell communication network
  #Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
  #(A) Compute and visualize the network centrality scores
  ptm = Sys.time()
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  #(B) Visualize dominant senders (sources) and receivers (targets) in a 2D space)
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  # Signaling role analysis on the cell-cell communication networks of interest
  gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("WNT", "NOTCH"))
  #> Signaling role analysis on the cell-cell communication network from user's input
  gg1 + gg2
  #(C) Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht1 + ht2

  #Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
  #(A) Identify and visualize outgoing communication pattern of secreting cells
  #Outgoing patterns reveal how the sender cells (i.e. cells as signal source) coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.
  library(NMF)
  library(ggalluvial)
  #Here we run selectK to infer the number of patterns.
  selectK(cellchat, pattern = "outgoing")
  #Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3
  nPatterns = 3
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  #river plot
  netAnalysis_river(cellchat, pattern = "outgoing")
  #> # dot plot
  netAnalysis_dot(cellchat, pattern = "outgoing")
  #(B) Identify and visualize incoming communication pattern of target cells
  #Incoming patterns show how the target cells (i.e. cells as signal receivers) coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signals.
  selectK(cellchat, pattern = "incoming")
  #Cophenetic values begin to drop when the number of incoming patterns is 4
  nPatterns = 4
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  #River plot
  netAnalysis_river(cellchat, pattern = "incoming")
  #> Please make sure you have load `library(ggalluvial)` when running this function
  #Dot plot
  netAnalysis_dot(cellchat, pattern = "incoming")


###Part V: Save the CellChat object
saveRDS(cellchat, file = "/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/cellchat_fertile.rds")
##Save the environment
