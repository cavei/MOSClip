availableOmicsMethods <- function() {
  return(c("summarizeModulesInCluster",
           "summarizeModulesToBinaryEvents",
           "summarizeModulesWithPca"))
}

summarizeModulesToBinaryEvents <- function(data, cliqueGenes, name="cov",
                                           binaryClassMin=10, cliques=NULL) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), cliqueGenes)
  if (length(genes)==0)
    return(NULL)

  dataClique <- t(data[genes, , drop=F])
  if (ncol(dataClique)==0)
    return(NULL)

  collapsed <- apply(dataClique>0, 1, any, na.rm=T)

  if (sum(collapsed) < binaryClassMin | sum(collapsed) > NROW(dataClique)-binaryClassMin)
    return(NULL)

  collapsed <- data.frame(collapsed, row.names = names(collapsed), stringsAsFactors = F)
  colnames(collapsed) <- name
  list(x=collapsed, dataModule=t(dataClique), namesCov=name, method="binary")
}

summarizeModulesInCluster <- function(data, cliqueGenes, name="clust", cliques=NULL) {
  datamat <- data$met ## modificare in modo che se non c'è un dizionario usi i geni stessi
  dict <- data$dict

  if (is.null(datamat))
    return(NULL)
  genes <- intersect(names(dict), cliqueGenes)

  if (length(genes)==0)
    return(NULL)

  used <- dict[genes]
  clusters <- unlist(dict[genes])
  datamatClique <- t(datamat[clusters, , drop=F])

  if (ncol(datamatClique)==0)
    return(NULL)

  hc <- hclust(dist(datamatClique, method = "euclidean"), method="ward.D2")
  # clusters <- kmeans(datamatClique, centers=2) # provioamo a implementare anche il Kmeans?

  if (ncol(datamatClique)<4){
    covs <- data.frame(factor(cutree(hc, k = 2)), stringsAsFactors = T)
    names(covs) <- paste0(name,"_2k")
  } else {
    covs <- data.frame(factor(cutree(hc, k = 3)), stringsAsFactors = T)
    names(covs) <- paste0(name,"_3k")
  }
  collapse=covs
  list(x=collapse, dataModule=t(datamatClique), namesCov=names(covs), cls=used, method="cluster")
}

summarizeModulesWithPca <- function(data, cliqueGenes, name="exprs", shrink=FALSE, method="regular", cliques=NULL, maxPCs=3, loadThr=0.6) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), cliqueGenes)
  if (length(genes)==0)
    return(NULL)

  dataClique <- t(data[genes, , drop=F])
  if (ncol(dataClique)==0)
    return(NULL)

  if (NCOL(dataClique)!=1) {
    pcs <- computePCs(dataClique, shrink=shrink, method=method, cliques=cliques, maxPCs=maxPCs)
  } else {
    colnames(dataClique) <- "PC1"
    pcs <- list(x=dataClique, sdev=sd(dataClique), loadings=1)
  }

  pcs$dataModule <- t(dataClique)
  pcs$method="pca"
  pcs$namesCov=colnames(pcs$x)
  pcs
}
