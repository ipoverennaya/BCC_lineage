## IMPORT DEPENDENCIES
library(ggplot2)
library(patchwork)
library(cowplot)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(viridis)
library(ggtree)

## DEFINE PATHS
path2data = "data/"
path2fig = "02_CellChatDB_cell_signaling/results/"


##-----------------------------------------------
## IMPORT DATA
##-----------------------------------------------

# crest_data = readRDS(paste0(path2data, "export/crest_data.RDS"))
# motor_data = readRDS(paste0(path2data, "export/motor_data.RDS"))
motor_crest = readRDS(paste0(path2data, "export/motor_crest_merged.RDS"))
motor_crest_E12 = readRDS(paste0(path2data, "export/motor_crest_E12_merged.RDS"))

Idents(motor_crest) = motor_crest$assignments
Idents(motor_crest_E12) = motor_crest_E12$assignments


##-----------------------------------------------
## n_Hub per Devtime
##-----------------------------------------------

hub_data = subset(motor_crest, subset = assignments == "Hub")
devtimes = unique(hub_data@meta.data$devtime)
nHub = lapply(devtimes, function(x) dim(hub_data@meta.data[hub_data@meta.data$devtime == x,])[1])
hub_count = data.frame("devtime" = unlist(devtimes), "count" = unlist(nHub))

hub_count$devtime = factor(hub_count$devtime, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E16.5", "E18.5", "P0", "P2", "P6", "P10", "Adult"))
hub_count = hub_count[order(hub_count$devtime),]

ggplot(hub_count) +
  geom_col(aes(x = devtime, y = count)) +
  labs(y = "# Hub") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(path2fig, "02_motor_crest_merged_Hub_state_prop.pdf"), width = 14, height = 10, units = "cm", dpi = 300)


##-----------------------------------------------
## VISUALIZATION
##-----------------------------------------------

## specify PW genes of interest
PW.SEMA3 = c("Sema3a", "Sema3b", "Sema3c", "Sema3d", "Sema3e", "Sema3f", "Plxna1", "Plxna2", "Plxna3", "Plxna4", "Nrp1", "Nrp2")
PW.NRG = c("Nrg1", "Nrg2", "Nrg3", "Nrg4", "Erbb2", "Erbb3", "Erbb4", "Itgav", "Itgb3")

### Dot plots
## merged object
DotPlot(motor_crest, features = PW.SEMA3) + labs(title = "PW: SEMA3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest/02_motor_crest_merged_PW_GOIs_SEMA3.pdf"), width = 18, height = 15, units = "cm", dpi = 300)

DotPlot(motor_crest, features = PW.NRG) + labs(title = "PW: NRG") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest/02_motor_crest_merged_PW_GOIs_NRG.pdf"), width = 15, height = 15, units = "cm", dpi = 300)

## E12.5 time-matched object
DotPlot(motor_crest_E12, features = PW.SEMA3) + labs(title = "PW: SEMA3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest_E12/02_motor_crest_E12_merged_PW_GOIs_SEMA3.pdf"), width = 18, height = 15, units = "cm", dpi = 300)

DotPlot(motor_crest_E12, features = PW.NRG) + labs(title = "PW: NRG") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest_E12/02_motor_crest_E12_merged_PW_GOIs_NRG.pdf"), width = 15, height = 15, units = "cm", dpi = 300)



### clustered Dot plots
DotPlot_hclust = function(seu, feat, width) {
  df = DotPlot(object = seu, features = feat)
  data = df$data
  
  data$origin = c("motor")
  dict = list("motor" = c("LMCm", "immature", "PGCa", "PGCb", "LMCl", "MMC", "HMC"),
              "crest" = c("BCC", "ChC", "enFib", "Gut_glia", "Gut_neuron", "Hub", "Melanocytes", "Mesenchyme", "NCC", "SatGlia", "SC", "Sensory", "Symp"))
  data$origin[which(data$id %in% dict$crest)] = "crest"
  
  # dendrogram
  mat = data %>%
    select(-avg.exp, -avg.exp.scaled, -origin) %>%
    pivot_wider(names_from = features.plot, values_from = pct.exp) %>%
    data.frame()
  row.names(mat) = mat$id
  mat = mat[,-1]
  clust = hclust(dist(mat %>% as.matrix()))
  
  ddgram = as.dendrogram(clust)
  ggtree_plot = ggtree::ggtree(ddgram, hang = -1)
  
  # dotplot
  dotplot = data %>%
    filter(pct.exp > 0) %>%
    mutate(id = factor(id, levels = clust$labels[clust$order])) %>%
    ggplot(aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) +
      geom_point() +
      # scale_color_gradient(low = "lightgrey", high = "blue") +
      scale_color_gradientn(colours = viridis::viridis(20), limits = c(-1.5,2.5), oob = scales::squish) +
      theme_cowplot() + 
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ylab("") + xlab("") +
    scale_y_discrete(position = "right")
  
  # annotation labels
  labels = ggplot(data %>%
                    mutate("Origin" = origin, id = factor(id, levels = clust$labels[clust$order])),
                  aes(y = id, x = 1, fill = origin)) +
    geom_tile() +
    scale_fill_brewer(palette = 'Set1') +
    theme_nothing() +
    theme(legend.position = "right")
  
  ggtree_plot + labels + dotplot +
    plot_layout(widths = c(0.7, 0.1, width), guides = "collect")
}


## merged object
DotPlot_hclust(motor_crest, feat = PW.SEMA3, width = 4) + labs(title = "PW: SEMA3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest/02_motor_crest_merged_PW_GOIs_SEMA3_clust.pdf"), width = 20, height = 15, units = "cm", dpi = 300)

DotPlot_hclust(motor_crest, feat = PW.NRG, width = 3) + labs(title = "PW: NRG") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest/02_motor_crest_merged_PW_GOIs_NRG_clust.pdf"), width = 17, height = 15, units = "cm", dpi = 300)


## E12.5 time-matched object
DotPlot_hclust(motor_crest_E12, feat = PW.SEMA3, width = 4) + labs(title = "PW: SEMA3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest_E12/02_motor_crest_E12_merged_PW_GOIs_SEMA3_clust.pdf"), width = 20, height = 15, units = "cm", dpi = 300)

DotPlot_hclust(motor_crest_E12, feat = PW.NRG, width = 3) + labs(title = "PW: NRG") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(path2fig, "motor_crest_E12/02_motor_crest_E12_merged_PW_GOIs_NRG_clust.pdf"), width = 17, height = 15, units = "cm", dpi = 300)

