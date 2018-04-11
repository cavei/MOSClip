plotModuleInGraph <- function(pathway, moduleNumber,
                        makeLegend=NULL, fileName=NULL) {

  if (is.null(makeLegend))
    makeLegend <- c(paste("omic",seq_len(length(pathway@modulesView[[moduleNumber]]))))
  net <- igraph.from.graphNEL(pathway@graphNEL)
  moduleGenes <- pathway@modules[[moduleNumber]]
  net <- simplify(net, remove.multiple = T, remove.loops = T)
  color <- rep("grey40", length(V(net)))
  color[names(V(net)) %in% moduleGenes] <- "red"
  V(net)$label <- NA
  V(net)$size <- 15
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber)
  mark.groups=lapply(involved, function(x) {
    row.names(x$subset)
  })
  colLength <- length(mark.groups)
  if (colLength<3) {
    mark.col=brewer.pal(3, "Set3")[seq_len(colLength)]
  } else {
    mark.col=brewer.pal(colLength, "Set3")
  }
  mark.border=NA
  entrez <- gsub("ENTREZID:", "", names(V(net)))

  symbol <- select(org.Hs.eg.db, keys=entrez,
                   columns = c("SYMBOL"), keytype="ENTREZID")$SYMBOL

  if (!is.null(fileName)) {
    pdf(fileName)
  }
  plot(net, edge.arrow.size=.2,vertex.label=symbol, vertex.label.cex=.5, edge.curved=.1,
       vertex.color=color, mark.groups=mark.groups, mark.col=mark.col, mark.border=NA)
  legend(x=-1, y=-1, makeLegend, pch=21, horiz=TRUE,
         col="#777777", pt.bg=mark.col, pt.cex=2, cex=.8, bty="n", ncol=1)
  if (!is.null(fileName)) {
    dev.off()
  }
}
