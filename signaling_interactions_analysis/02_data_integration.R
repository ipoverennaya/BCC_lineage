## IMPORT DEPENDENCIES
library(patchwork)
library(Seurat)
library(insight)

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
# motor_crest_int = readRDS(paste0(path2data, "motor_crest_integrated.RDS"))
# motor_crest_E12_int = readRDS(paste0(path2data, "motor_crest_E12_integrated.RDS"))

## subset data for time-matching devtime E12.5
crest_data_E12 = subset(crest_data, subset = devtime == "E12.5")

## merge datasets
motor_crest = merge(motor_data, crest_data, add.cell.ids = c("motor", "crest"), merge.data = T, project = "motor_crest")
motor_crest_E12 = merge(motor_data, crest_data_12, add.cell.ids = c("motor", "crest"), merge.data = T, project = "motor_crest_E12")
print_color("Data import successful.\n", "green")


##-----------------------------------------------
## DATA INTEGRATION
##-----------------------------------------------

seu = motor_crest
# seu = motor_crest_E12

## seurat v4 workflow
data_list = SplitObject(seu, split.by = "orig.ident")
print_color("Object list created.\n", "green")

features = SelectIntegrationFeatures(object.list = data_list)
print_color("Integration features selected.\n", "green")

int_data = FindIntegrationAnchors(object.list = data_list)
print_color("Integration anchors found.\n", "green")

data_integrate = IntegrateData(int_data)
print_color("Data integration successful.\n", "green")

DefaultAssay(data_integrate) = "integrated"
data_integrate_processed = ScaleData(data_integrate)
data_integrate_processed = RunPCA(data_integrate_processed)
data_integrate_processed = RunUMAP(data_integrate_processed, dims = 1:50)
data_integrate_processed = FindNeighbors(data_integrate_processed, dims = 1:50)
data_integrate_processed = FindClusters(data_integrate_processed, resolution = 1.3)
print_color("Processing successful.\n", "green")

saveRDS(data_integrate_processed, paste0(path2data, "motor_crest_integrated.RDS"))
# saveRDS(data_integrate_processed, paste0(path2data, "motor_crest_E12_integrated.RDS"))
print_color("Integrated object stored.\n", "green")

print_color("done\n", "green")


##-----------------------------------------------
## VISUALIZE
##-----------------------------------------------

(DimPlot(data_integrate_processed, group.by = "assignments", label = F, label.size = 3, cols = CT_color_dict_man) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5))) +
(DimPlot(data_integrate_processed, group.by = "orig.ident", cols = "Set1") + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)))

# ggsave(paste0(path2fig, "motor_crest/02_motor_crest_integrated_UMAP.pdf"), width = 30, height = 15, units = "cm", dpi = 300)
# ggsave(paste0(path2fig, "motor_crest_E12/02_motor_crest_E12_integrated_UMAP.pdf"), width = 30, height = 15, units = "cm", dpi = 300)

