## IMPORT DEPENDENCIES
library(patchwork)
library(Seurat)
library(harmony)

set.seed(123)

## DEFINE PATHS
path2data = "data/export/"
path2fig = "01_QC_and_data_preprocessing/results/"

## SOURCE ADDITIONAL SCRIPTS
source("src/colors.R")


##-----------------------------------------------
## IMPORT DATA
##-----------------------------------------------

crest_data = readRDS(paste0(path2data, "crest_data.RDS"))
motor_data = readRDS(paste0(path2data, "motor_data.RDS"))

## subset data for time-matching devtime E12.5
crest_data_E12 = subset(crest_data, subset = devtime == "E12.5")

## merge datasets
motor_crest = merge(motor_data, crest_data, add.cell.ids = c("motor", "crest"), merge.data = T, project = "motor_crest")
motor_crest_E12 = merge(motor_data, crest_data_E12, add.cell.ids = c("motor", "crest"), merge.data = T, project = "motor_crest_E12")

Idents(motor_crest) = motor_crest$assignments
Idents(motor_crest_E12) = motor_crest_E12$assignments

## standard seurat processing workflow
seu_processing = function(seu) {
  seu = NormalizeData(seu)
  seu = FindVariableFeatures(seu)
  seu = ScaleData(seu)
  seu = RunPCA(seu)
  seu = RunUMAP(seu, dims = 1:50)
  seu = FindNeighbors(seu, dims = 1:50)
  seu = FindClusters(seu, resolution = 1.3)
  return(seu)
}

motor_crest = seu_processing(motor_crest)
motor_crest_E12 = seu_processing(motor_crest_E12)


## perform harmony batch effect correction
run_harmony = function(seu) {
  seu = RunHarmony(seu, "orig.ident")
  seu = RunUMAP(seu, dims = 1:50, reduction = "harmony")
  return(seu)
}

motor_crest_harm = run_harmony(motor_crest)
motor_crest_E12_harm = run_harmony(motor_crest_E12)

saveRDS(motor_crest, paste0(path2data, "motor_crest_merged.RDS"))
saveRDS(motor_crest_E12, paste0(path2data, "motor_crest_E12_merged.RDS"))


##-----------------------------------------------
## VISUALIZE
##-----------------------------------------------

# get colors from CT_color_dict_man found in 01_CellChat_signaling_analysis.R

## merged datasets
(DimPlot(motor_crest, group.by = "assignments", label = F, label.size = 3, cols = CT_color_dict_man) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5))) +
(DimPlot(motor_crest, group.by = "orig.ident", cols = "Set1") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)))
ggsave(paste0(path2fig, "motor_crest/03_motor_crest_merged_UMAP.pdf"), width = 30, height = 15, units = "cm", dpi = 300)

(DimPlot(motor_crest_E12, group.by = "assignments", label = F, label.size = 3, cols = CT_color_dict_man) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5))) +
(DimPlot(motor_crest_E12, group.by = "orig.ident", cols = "Set1") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)))
ggsave(paste0(path2fig, "motor_crest_E12/03_motor_crest_E12_merged_UMAP.pdf"), width = 30, height = 15, units = "cm", dpi = 300)

## harmony corrected merged datasets
(DimPlot(motor_crest_harm, group.by = "assignments", label = F, label.size = 3, cols = CT_color_dict_man) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5))) +
(DimPlot(motor_crest_harm, group.by = "orig.ident", cols = "Set1") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)))
ggsave(paste0(path2fig, "motor_crest/03_motor_crest_merged_BECharmony_UMAP.pdf"), width = 30, height = 15, units = "cm", dpi = 300)

(DimPlot(motor_crest_E12_harm, group.by = "assignments", label = F, label.size = 3, cols = CT_color_dict_man) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5))) +
(DimPlot(motor_crest_E12_harm, group.by = "orig.ident", cols = "Set1") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)))
ggsave(paste0(path2fig, "motor_crest_E12/03_motor_crest_E12_merged_BECharmony_UMAP.pdf"), width = 30, height = 15, units = "cm", dpi = 300)

