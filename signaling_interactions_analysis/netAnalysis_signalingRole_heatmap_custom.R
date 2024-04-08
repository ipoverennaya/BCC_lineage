netAnalysis_signalingRole_heatmap_custom = function (object, signaling = NULL, PW_filter = NULL, CT_filter = NULL,
                                                     pattern = c("outgoing", "incoming", "all"),
                                                     slot.name = "netP", color.use = NULL, CT_colors = NULL,
                                                     color.heatmap = "BuGn", title = NULL, width = 10, height = 8, 
                                                     font.size = 8, font.size.title = 10, cluster.rows = FALSE, 
                                                     cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  
  # subset matrix
  mat.ori <- mat
  
  PW_filter = intersect(PW_filter, rownames(mat))
  CT_filter = intersect(CT_filter, colnames(mat))
  
  row_indices = which(rownames(mat) %in% PW_filter)
  col_indices = which(colnames(mat) %in% CT_filter)

  if (!is.null(PW_filter) && !is.null(CT_filter)) {
    mat <- mat[row_indices, col_indices]
  } else if (!is.null(PW_filter) && is.null(CT_filter)) {
    mat <- mat[row_indices, ]
  } else if (!is.null(CT_filter) && is.null(PW_filter)) {
    mat <- mat[, col_indices]
  } else {
    mat <- mat
  }
  
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  
  if (!is.null(CT_colors)) {
    color.use <- CT_colors[names(CT_colors) %in% colnames(mat)]
    color.use <- color.use[match(colnames(mat), names(color.use))]
    names(color.use) <- colnames(mat)
  } else {
    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
      names(color.use) <- colnames(mat)
    }
  }

  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(100)
  
  df <- data.frame(celltype = colnames(mat))
  rownames(df) <- colnames(mat)
  
  # Adding orig.ident annotation
  orig.ident = c()
  for (CT in colnames(mat)) {
    orig.ident = c(orig.ident, unique(object@meta[object@meta$assignments == CT,]$orig.ident))
  }
  df$origin = orig.ident
  
  col_annotation <- HeatmapAnnotation(df = df, col = list(origin = c("Crest" = "#E41A1C", "Motor" = "#377EB8"),
                                                          celltype = color.use),
                                      gap = unit(0.25, "cm"), annotation_name_side = "left", which = "column", show_annotation_name = T,
                                      show_legend = c(F, T), annotation_legend_param = list(origin = list(title = "Origin")),
                                      annotation_name_gp = gpar(fontsize = 8, fontface = "bold"), simple_anno_size = grid::unit(0.25, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(if (!is.null(CT_filter)) colSums(mat.ori)[col_indices] else colSums(mat.ori), 
                          border = FALSE, gp = gpar(fill = color.use, col = color.use)), show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  ha1 = rowAnnotation(Strength = anno_barplot(if (!is.null(PW_filter)) pSum[row_indices] else pSum,
                                              border = FALSE), show_annotation_name = FALSE)
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                cluster_rows = cluster.rows, cluster_columns = cluster.cols,
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size, fontface = "bold"),
                column_names_rot = 60, column_names_gp = gpar(fontsize = font.size, fontface = "bold"),
                column_title = title, column_title_gp = gpar(fontsize = font.size.title, fontface = "bold"),
                width = unit(width, "cm"), height = unit(height, "cm"),
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                            border = NA, at = legend.break, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8),
                                            grid_width = unit(2, "mm")))
  return(ht1)
}

