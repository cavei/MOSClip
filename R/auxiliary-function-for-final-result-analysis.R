getPathFathers<- function(p, df.g, ord=3, plot=F) {
  if (is.null(roots)) {
    roots=names(V(df.g))
  }
  require(igraph)
  if (!(p %in% names(V(df.g)))){
    warning(paste0("Id ", p, " is not in the hierarchy."))
    return(p)
  }

  mm <- make_ego_graph(df.g, ord, nodes = p, mode="in")
  mmlist <- ego(df.g, ord, nodes = p, mode="in")

  if (plot)
    plot(mm[[1]])

  chain <- as_ids(mmlist[[1]])
  parents <- chain[-1]
  if (length(parents)==0)
    return(chain)

  dis <- distances(mm[[1]])
  idx <- which.max(dis[p, ])
  return(colnames(dis)[idx])

  # idx <- which(max(dis) == dis, arr.ind = T)[, 'row']
  # fathers <- colnames(dis)[idx]
  # fathers <- row.names(which(max(dis) == dis, arr.ind = T))
  # fathers <- setdiff(fathers, chain[1])
  # refinedFathers <- intersect(fathers, roots)
  # if (length(refinedFathers) > 0)
  #   fathers <- refinedFathers
  # fathers
}

findRoots <- function(pathHierarchy) {
  setdiff(pathHierarchy[,1], pathHierarchy[,2])
}

createIntersection <- function(list, universeDim) {
  require(SuperExactTest)
  require(RColorBrewer)
  tmp <- lapply(list, function(x) {
    names(which(x == 1))
  })
  a<-supertest(tmp, n=universeDim)
  plot(a,color.on="#409ec3",color.off="white",
       heatmapColor=colorRampPalette(brewer.pal(9,"OrRd"))(100))
  # plot(a, Layout="landscape", degree=2:4, sort.by="size", y.pos=c(0.025,0.95))
  a
}

mapReactomeIDfromGraphite <- function(reactome, pathwayNames=NULL) {
  if (is.null(pathwayNames)) {
    pathwayNames <- names(reactome)
  }
  l <- lapply(pathwayNames, function(p) {
    if (!(p %in% names(reactome))) {
      warning(paste0("No id found for ", p ))
      return(NULL)
    }
    data.frame(id=reactome[[p]]@id, pname=p, stringsAsFactors = F)
  })
  do.call(rbind, l)
}
