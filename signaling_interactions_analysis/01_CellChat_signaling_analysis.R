## IMPORT DEPENDENCIES
library(ggplot2)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
library(Seurat)
library(CellChat)
library(igraph)
library(RColorBrewer)
library(viridis)

set.seed(123)

## DEFINE PATHS
path2data = "data/"
path2fig = "02_CellChatDB_cell_signaling/results/"

## SOURCE ADDITIONAL SCRIPTS
source("src/CellChat_custom/netAnalysis_signalingRole_heatmap_custom.R")
source("src/colors.R")


##-----------------------------------------------
## IMPORT DATA
##-----------------------------------------------

## load datasets
# crest_data = readRDS(paste0(path2data, "export/crest_data.RDS"))
# motor_data = readRDS(paste0(path2data, "export/motor_data.RDS"))
motor_crest = readRDS(paste0(path2data, "export/motor_crest_integrated.RDS"))
motor_crest_E12 = readRDS(paste0(path2data, "export/motor_crest_E12_integrated.RDS"))

Idents(motor_crest) = motor_crest$assignments
Idents(motor_crest_E12) = motor_crest_E12$assignments

## subset data for time-matching devtime E12.5
# crest_data_E12 = subset(crest_data, subset = devtime == "E12.5")


##-----------------------------------------------
## CELL CHAT - CELL SIGNALING ANALYSIS
##-----------------------------------------------

## set ligand-receptor interaction database
CellChatDB = CellChatDB.mouse

## correct CellChatDB typos
CellChatDB$interaction["H2-BI_KIR3DL1","interaction_name"] = "H2-BL_KIR3DL1"
CellChatDB$interaction["H2-BI_KIR3DL1","ligand"] = "H2-Bl"
CellChatDB$interaction["H2-BI_KIR3DL1","interaction_name_2"] = "H2-bl - Kir3dl1"
rownames(CellChatDB$interaction)[rownames(CellChatDB$interaction) == "H2-BI_KIR3DL1"] = "H2-BL_KIR3DL1"

CellChatDB$interaction["H2-EA-PS_CD4","interaction_name"] = "H2-EA_CD4"
CellChatDB$interaction["H2-EA-PS_CD4","ligand"] = "H2-Ea"
CellChatDB$interaction["H2-EA-PS_CD4","interaction_name_2"] = "H2-ea - Cd4"
rownames(CellChatDB$interaction)[rownames(CellChatDB$interaction) == "H2-EA-PS_CD4"] = "H2-EA_CD4"

CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling")


##-----------------------------------------------
## CELL CHAT - STANDARD PROCESSING WORKFLOW
##-----------------------------------------------

runCellChat = function(cellchat) {
  cellchat@DB = CellChatDB.use
  cellchat = subsetData(cellchat)
  
  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)
  cellchat = projectData(cellchat, PPI.mouse)
  cellchat = computeCommunProb(cellchat)
  cellchat = filterCommunication(cellchat, min.cells = 10)
  cellchat = computeCommunProbPathway(cellchat)
  cellchat = aggregateNet(cellchat)
  
  cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  return(cellchat)
}

### CREATE CELLCHAT OBJECTS
## original object
# motor_data_CC = createCellChat(motor_data, assay = "RNA") # normalized motor_data
# crest_data_CC = createCellChat(crest_data, assay = "RNA") # normalized crest_data
motor_crest_CC = createCellChat(motor_crest, assay = "RNA") # merged object

## time-matched object
# crest_data_E12_CC = createCellChat(crest_data_E12, assay = "RNA") # normalized crest_data
motor_crest_E12_CC = createCellChat(motor_crest_E12, assay = "RNA") # merged time-matched object


### RUN CELLCHAT
## original object
# crest_data_CC = runCellChat(crest_data_CC)
# motor_data_CC = runCellChat(motor_data_CC)

## time-matched object
# crest_data_E12_CC = runCellChat(crest_data_E12_CC)
motor_crest_CC = runCellChat(motor_crest_CC)
motor_crest_E12_CC = runCellChat(motor_crest_E12_CC)


# ### MERGE CELLCHAT OBECTS: perform merging only after individual CellChat processing
# ## original object
# object.list = list(motor = motor_data_CC, crest = crest_data_CC)
# motor_crest_CCmerge = mergeCellChat(object.list = object.list, add.names = names(object.list))
# 
# ## time-matched object
# object.list = list(motor = motor_data_CC, crest_E12 = crest_data_E12_CC)
# motor_crest_E12_CCmerge = mergeCellChat(object.list = object.list, add.names = names(object.list))
# 
# motor_data_CC_lift = liftCellChat(motor_data_CC, group.new = levels(crest_data_CC@idents))
# crest_data_CC_lift = liftCellChat(crest_data_CC, group.new = levels(motor_data_CC@idents))




##-----------------------------------------------
## VISUALIZATION
##-----------------------------------------------

# cellchat = motor_crest_CC
cellchat = motor_crest_E12_CC

### import CellChat ComNet filters
motorcrest_PWfilter1 = as.vector(unlist(read.csv(paste0(path2data, "cellchat/motorcrest_PWfilter1.csv"), header = F))) # to keep
motorcrest_CTfilter1 = as.vector(unlist(read.csv(paste0(path2data, "cellchat/motorcrest_CTfilter1.csv"), header = F))) # to eliminate

motorcrest_CTfilter1 = setdiff(cellchat@meta$ident, motorcrest_CTfilter1)
setdiff(cellchat@netP$pathways, motorcrest_PWfilter1)

### visualize
plot_netVisual_circle = function(cellchat) {
  groupSize = as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
}

plot_netVisual_circle(cellchat)


plot_individual_netVisual_circle = function(cellchat) {
  groupSize = as.numeric(table(cellchat@idents))
  mat = cellchat@net$weight
  par(mfrow = c(5,4), xpd = TRUE)
  for (i in 1:nrow(mat)) {
    mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] = mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}

plot_individual_netVisual_circle(cellchat)
par(mfrow = c(1, 1))


## IO signaling with applied PW & CT Filter 1 (just filtering the plotting output)
IO_signaling_custom = function(cellchat, height, PW_filter, CT_filter, colors) {
  ht_opt$ROW_ANNO_PADDING = unit(2, "mm")
  ht_opt$COLUMN_ANNO_PADDING = unit(2, "mm")
  
  ht1 = netAnalysis_signalingRole_heatmap_custom(cellchat, pattern = "outgoing", height = height, cluster.rows = F, cluster.cols = T,
                                                 PW_filter = PW_filter, CT_filter = CT_filter, CT_colors = colors,
                                                 font.size.title = 12, font.size = 10)
  ht2 = netAnalysis_signalingRole_heatmap_custom(cellchat, pattern = "incoming", height = height, cluster.rows = F, cluster.cols = T,
                                                 PW_filter = PW_filter, CT_filter = CT_filter, CT_colors = colors,
                                                 font.size.title = 12, font.size = 10)
  draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
}

IO_signaling_custom(cellchat, height = 12, PW_filter = motorcrest_PWfilter1, CT_filter = motorcrest_CTfilter1, colors = CT_color_dict_man)


## top PW (LR-pair) contributions
netAnalysis_contribution(cellchat, signaling = "SEMA3")


## circle/chord plot
origin = c()
for (CT in motor_crest_CT) {
  orig = unique(cellchat@meta[cellchat@meta$assignments == CT,]$orig.ident)
  origin = c(origin, orig)
}
celltypes = intersect(motor_crest_CT, cellchat@meta$assignments)
names(origin) = celltypes
CT.show = intersect(celltypes, motorcrest_CTfilter1)

pathways.show = c("SEMA3", "NRG", "MIF", "ENHO", "GAS")
par(mfrow = c(1, length(pathways.show)), xpd=TRUE)#, mar = c(1, 1, 1, 1))
for (PW in pathways.show) {
  # netVisual_aggregate(cellchat, signaling = PW, layout = "circle", color.use = CT_color_dict_man, cell.order = motor_crest_CT,
  #                     edge.width.max = 5, remove.isolate = T, sources.use = CT.show, targets.use = CT.show, group = origin,
  #                     small.gap = 3, big.gap = 15, vertex.label.cex = 1)
  netVisual_chord_cell(cellchat, signaling = PW, group = origin, color.use = CT_color_dict_man, cell.order = motor_crest_CT,
                       small.gap = 3, big.gap = 15, remove.isolate = T, lab.cex = 1, sources.use = CT.show, targets.use = CT.show)
}
par(mfrow = c(1, 1))


## Feature Plots
features.use = c("Phox2b", "Sox8", "Sox10", "Sox11", "Itga4", "Nell2", "Serpine2", "Ets1", "Tfap2a", "Tfap2b", "Sparc")
FeaturePlot(motor_crest_harm, features = features.use, pt.size = 0.5, coord.fixed = T)

