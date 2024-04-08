## IMPORT DEPENDENCIES
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(Seurat)
library(SeuratDisk)

## DEFINE PATHS
path2data = "data/"
path2fig = "01_QC_and_data_preprocessing/results/"

## SOURCE ADDITIONAL SCRIPTS
source("src/celltype_marker_profiles.R")
source("src/plotting_functions.R")


##-----------------------------------------------
## IMPORT DATA
##-----------------------------------------------

## motor
motor_assignments = as.vector(unlist(read.csv(paste0(path2data, "motor_data_CT_assignments.csv"))))
motor_data = readRDS(paste0(path2data, "Motor.rds"))
motor_data = UpdateSeuratObject(motor_data)
motor_data@meta.data[["devtime"]] = "E12.5"
motor_data@meta.data[["orig.ident"]] = "Motor"
motor_data@meta.data[["assignments"]] = motor_assignments

## crest
Convert(paste0(path2data, "crest_adata_assigned.h5ad"), dest = "h5seurat", overwrite = TRUE)
crest_data = LoadH5Seurat(paste0(path2data, "crest_adata_assigned.h5seurat"), assays = "RNA")
# crest_data[["percent.mt"]] = PercentageFeatureSet(crest_data, pattern = "^MT-")
crest_data = UpdateSeuratObject(crest_data)
crest_data@meta.data[["orig.ident"]] = "Crest"

assignments_new = as.character(crest_data@meta.data$assignments)
assignments_new[assignments_new == "none"] = "Hub"
assignments_new = as.factor(assignments_new)
crest_data@meta.data$assignments = assignments_new

## set default assay to RNA
DefaultAssay(crest_data) = "RNA"
DefaultAssay(motor_data) = "RNA"


##-----------------------------------------------
## Quality Control
##-----------------------------------------------

VlnPlot(object = motor_data, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
VlnPlot(object = crest_data, features = c("total_counts", "n_genes_by_counts"))

Idents(crest_data) = crest_data@meta.data$assignments
Idents(motor_data) = motor_data@meta.data$assignments

(DimPlot(object = crest_data, label = T, group.by = "assignments") + labs(title = "Crest", color = "Assignment") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5))) +
(DimPlot(object = crest_data, label = T, group.by = "leiden") + labs(title = "Crest", color = "Leiden Clustering") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)))


##-----------------------------------------------
## SANITY CHECK: FEATURE PLOTS
##-----------------------------------------------

cellProfileFeatPlot(crest_data, crest_celltype_marker_profiles, "crest_sanity_FP/")

# patchwork::wrap_plots(FeaturePlot(object = crest_data, features = crest_celltype_marker_profiles[[1]], ncol = 4, coord.fixed = T)) + 
#   plot_annotation(title = names(crest_celltype_marker_profiles[1])) & theme(plot.title = element_text(face = 'bold', size = 20))


##-----------------------------------------------
## MANUAL CELLTYPE RE-LABELING
##-----------------------------------------------

# ## crest
# crest_meta_CT_subset = crest_data@meta.data[,11:22]
# crest_meta_CT_subset$index = rownames(crest_meta_CT_subset)
# 
# crest_CT = crest_meta_CT_subset %>%
#   gather(celltype, CT, -index) %>%
#   filter(CT == TRUE) %>%
#   select(-CT)
# 
# ind = which(rownames(crest_data@meta.data) %in% crest_CT$index)
# crest_data@meta.data[["celltype"]][ind] = crest_CT$celltype
# 
# 
# ## celltype substitution list according to leiden clustering
# CTsub = list("mesenchyme" = c(1, 17),
#              "cranial neural crest" = c(2),
#              "trunk neural crest" = c(6, 7, 9),
#              "sensory neurons" = c(10, 12, 22),
#              "schwann cells" = c(8, 13, 15, 23, 24, 26, 27),
#              "ChC" = c(5, 28),
#              "sympathetic neurons" = c(0, 21),
#              "enteric neurons" = c(20, 29),
#              "enteric glia" = c(18),
#              "satellite glia" = c(16),
#              "BCC" = c(),
#              "melanocytes" = c(25),
#              "endoneurial fibroblasts" = c(),
#              "Hub" = c(3, 4, 11, 14, 19)
#            )
# 
# leiden2CT = function(seu, CT_sub_list) {
#   leiden_clust = as.vector(seu@meta.data$leiden)
#   celltype = leiden_clust
#   for (i in 1:length(CT_sub_list)) {
#     ind = which(leiden_clust %in% CT_sub_list[[i]])
#     celltype[ind] = names(CT_sub_list[i])
#   }
#   seu@meta.data[["celltype"]] = celltype
#   return(seu)
# }
# 
# crest_data = leiden2CT(crest_data, CTsub)


##-----------------------------------------------
## NORMALIZE
##-----------------------------------------------

# ## SS2
# crest_data = NormalizeData(crest_data)
# crest_data = FindVariableFeatures(crest_data)
# crest_data = ScaleData(crest_data)
# 
# ## 10X
# # motor_data = SCTransform(motor_data)
# motor_data = NormalizeData(motor_data)
# motor_data = FindVariableFeatures(motor_data)
# motor_data = ScaleData(motor_data)


##-----------------------------------------------
## EXPORT DATA
##-----------------------------------------------

if(!dir.exists(paste0(path2data, "export/"))) {
  dir.create(paste0(path2data, "export/"))
}

saveRDS(crest_data, paste0(path2data, "export/crest_data.RDS"))
saveRDS(motor_data, paste0(path2data, "export/motor_data.RDS"))




