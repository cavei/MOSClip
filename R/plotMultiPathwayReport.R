plotMultiPathwayReport <- function(multiPathwayList, top=25){
  if(!is.list(multiPathwayList))
    stop("multiPathwayList must be a list.")

  summary <- multiPathwayReport(multiOmicsFull)
  top <- min(top, NROW(summary))
  pheatmap(summary[seq_len(top),,drop=F], display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F, cluster_cols = F, gaps_col = c(1),
           fontsize_row=5)
}

plotModuleReport <- function(pathwayObj) {
  summary <- formatModuleReport(pathwayObj)
  pheatmap(summary, display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F,
           cluster_cols = F,
           gaps_col = 1)
}
