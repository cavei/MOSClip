plotMultiPathwayReport <- function(multiPathwayList, top=25, p.adjust.method="BH"){
  require(gridExtra)
  library(RColorBrewer)
  library(pheatmap)

  if(!is.list(multiPathwayList))
    stop("multiPathwayList must be a list.")

  summary <- multiPathwayReport(multiOmicsFull, p.adjust.method = p.adjust.method)
  top <- min(top, NROW(summary))
  pheatmap(summary[seq_len(top),,drop=F], display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F, cluster_cols = F, gaps_col = c(1,2),
           fontsize_row=5)
}

plotMultiPathwayModuleReport <- function(multiPathwayList, top=25, p.adjust.method="BH"){
  require(gridExtra)
  library(RColorBrewer)
  library(pheatmap)

  if(!is.list(multiPathwayList))
    stop("multiPathwayList must be a list.")

  summary <- multiPathwayReport(multiPathwayList, p.adjust.method = p.adjust.method)
  top <- min(top, NROW(summary))
  pheatmap(summary[seq_len(top),,drop=F], display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F, cluster_cols = F, gaps_col = c(1,2),
           fontsize_row=5)
}

