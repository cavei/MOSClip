plotModuleReport <- function(pathwayObj) {
  summary <- formatModuleReport(pathwayObj)
  pheatmap(summary, display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F,
           cluster_cols = F,
           gaps_col = 1)
}
