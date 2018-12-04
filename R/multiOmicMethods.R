#' Get available Omics Summarizing Methods
#'
#' Gives a vector of the available methods to summarize omics.
#'
#' @return character vector with the implemented methods.
#'
#' @export
availableOmicMethods <- function() {
  return(c("summarizeToBinaryEvents",
           "summarizeToNumberOfEvents",
           "summarizeWithPca",
           "summarizeInCluster",
           "summarizeInClusterWithoutDictionary",
           "summarizeToBinaryDirectionalEvents",
           "summarizeToNumberOfDirectionalEvents"))
}

#' Summarize To Binary Events
#'
#' Given a matrix it summarize to a 0 or 1
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param binaryClassMin the minimum number of event to include the covariate
#' @param cliques the features organized in cliques. Only use for topology.
#'
#' @return NULL
#'
#' @export
summarizeToBinaryEvents <- function(data, features, name="bin",
                                           binaryClassMin=10, cliques=NULL) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
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
  list(x=collapsed, dataModule=t(dataClique), namesCov=name, method="binary", omicName=name, eventThr=1)
}

#' Summarize To Number of Binary Events
#'
#' Given a matrix it summarize to a 0 or 1
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param min_prop minimal proportion in classes
#' @param cliques the features organized in cliques. Only use for topology.
#'
#' @return NULL
#'
#' @export
summarizeToNumberOfEvents <- function(data, features, name="event", min_prop=0.1, cliques=NULL) {
  if (is.null(data))
    return(NULL)
  
  genes <- intersect(row.names(data), features)
  if (length(genes)==0)
    return(NULL)
  
  dataClique <- t(data[genes, , drop=F])
  if (ncol(dataClique)==0)
    return(NULL)
  
  collapsed <- apply(dataClique>0, 1, sum, na.rm=T)
  
  # min <- ceiling(NCOL(data)*0.01)
  # if (sum(collapsed) <  min | (sum(collapsed)) > NROW(dataClique)-min)
  #   return(NULL)
  
  keep <- check_minimal_proportion(collapsed, min_prop=min_prop)
  if (!keep)
    return(NULL)
  
  collapsed <- data.frame(collapsed, row.names = names(collapsed), stringsAsFactors = F)
  colnames(collapsed) <- name
  list(x=collapsed, dataModule=t(dataClique), namesCov=name, method="count", omicName=name,
       eventThr = 1, min_prop=min_prop)
}

#' Summarize Using Cluster Analysis
#'
#' Given a matrix it summarize in classes
#' 
#' The user can define a maximum of classes. The function
#' guess the optimal number of clusters using NbClust methods.
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param dictionary translate features (genes) into sets (row.names of the data)
#' @param max_cluster_number the maximum number of cluster to evaluate
#' @param cliques the features organized in cliques. Only use for topology
#'
#' @return NULL
#' @importFrom stats cutree dist hclust
#' @importFrom NbClust NbClust
#' @export
summarizeInCluster <- function(data, features, name="clu", dictionary=NULL, max_cluster_number=3, cliques=NULL) {
  if (is.null(data) || (ncol(data)==0) || !(is.matrix(data)))
    return(NULL)

  if (is.null(dictionary)) {
    genes <- intersect(rownames(data), features)
    if (length(genes)==0)
      return(NULL)
    used <- genes
    names(used) <- genes
    datamatClique <- t(data[genes, ,drop=FALSE])
  } else {
    genes <- intersect(names(dictionary), features)
    if (length(genes)==0)
      return(NULL)
    used <- dictionary[genes]
    clusters <- unlist(dictionary[genes])
    clusters <- intersect(clusters, row.names(data))
    if (length(clusters)==0)
      return(NULL)
    datamatClique <- t(data[clusters, , drop=F])
  }

  if (ncol(datamatClique)==0)
    return(NULL)

  ## CREATE CLUSTERS
  if (TRUE){
    covs <- createOptiomalClusterClasses(datamatClique, name, max_cluster_number = max_cluster_number)
  } else {
    covs <- createClusterClassesOld(datamatClique, name)
  }
  
  if (any(table(covs[[1]])<2)){
    warning("Not meaningful class separation\n")
    return(NULL)
  }
    
  collapse=covs
  list(x=collapse, dataModule=t(datamatClique), namesCov=names(covs), cls=used, method="cluster", omicName=name)
}

createOptiomalClusterClasses <- function(datamatClique, name, max_cluster_number, index_method="silhouette") {
  nb <- sinkNbClust(data=datamatClique, min.nc=2, max.nc=max_cluster_number, method="ward.D2", index=index_method)
  covs <- data.frame(factor(nb$Best.partition), stringsAsFactors = T)
  optimalCLusterNumber <- length(table(nb$Best.partition))
  names(covs) <- paste0(name,optimalCLusterNumber, "k")
  covs
}

sinkNbClust <- function(data, min.nc=2, max.nc=6, method="ward.D2", index="silhouette"){
  if (index=="all")
    sink(file=tempfile()); pdf(file=NULL)

  nb <- NbClust(data=data, min.nc=min.nc, max.nc=max.nc, method=method, index=index)

  if (index=="all")
    sink(); dev.off()
  
  return(nb)
}

createClusterClassesOld <- function(datamatClique, name){
  warning("you are using an old way to evaluate cluster...\n
          this function will be deprecated soon...")
  
  md <- dist(datamatClique, method = "euclidean")
  if (any(is.na(md)))
    return(NULL)
  hc <- hclust(md, method="ward.D2")
  # clusters <- kmeans(datamatClique, centers=2) # provioamo a implementare anche il Kmeans?
  
  if (ncol(datamatClique)<4){
    covs <- data.frame(factor(cutree(hc, k = 2)), stringsAsFactors = T)
    names(covs) <- paste0(name,"2k")
  } else {
    covs <- data.frame(factor(cutree(hc, k = 3)), stringsAsFactors = T)
    names(covs) <- paste0(name,"3k")
  }
  return(covs)
}

#' Summarize Using Cluster Analysis with no dictionary to translate the matrix ids
#'
#' Given a quantitative matrix it clusters the samples 
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param cliques the features organized in cliques. Only use for topology.
#'
#' @return NULL
#' @importFrom stats cutree dist hclust
#' @export
summarizeInClusterWithoutDictionary <- function(data, features, name="clu", cliques=NULL) {
  warning("function summarizeInClusterWithoutDictionary has been deprecated...\n
          use summarizeInCluster")

  if (is.null(data) | (ncol(data)==0) | !(is.matrix(data)))
    return(NULL)
  genes <- intersect(rownames(data), features)
  
  if (length(genes)==0)
    return(NULL)
  
  datamatClique <- t(data[genes, ,drop=FALSE])
  
  used <- colnames(datamatClique)
  names(used) <- colnames(datamatClique)
  
  md <- dist(datamatClique, method = "euclidean")
  if (any(is.na(md)))
      return(NULL)
  hc <- hclust(md, method="ward.D2")
  # clusters <- kmeans(datamatClique, centers=2) # TO DO: add kmeans
  
  if (ncol(datamatClique)<4){
    covs <- data.frame(factor(cutree(hc, k = 2)), stringsAsFactors = T)
    names(covs) <- paste0(name,"2k")
  } else {
    covs <- data.frame(factor(cutree(hc, k = 3)), stringsAsFactors = T)
    names(covs) <- paste0(name,"3k")
  }
  
  if (any(table(covs[[1]])<2)){
    # warning("Not meaningful class separation\n")
    return(NULL)
  }
  
  list(x=covs, dataModule=t(datamatClique), namesCov=names(covs), cls=used, method="cluster", omicName=name)
}

#' Summarize Using PCA
#'
#' Given a matrix it summarize to a 0 or 1
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param shrink shirnk or not the covariance matrix.
#' @param method either "regular", "sparse" or "topological"
#' @param cliques the features organized in cliques. Only use for topology.
#' @param maxPCs maximum number of pcs to consider
#' @param loadThr loading threshold
#'
#' @return NULL
#' @importFrom stats sd
#' @importFrom houseOfClipUtility computePCs
#' @export
summarizeWithPca <- function(data, features, name="pca", shrink=FALSE, method="regular", cliques=NULL, maxPCs=3, loadThr=0.6) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
  if (length(genes)==0)
    return(NULL)

  dataClique <- t(data[genes, , drop=F])
  if (ncol(dataClique)==0)
    return(NULL)

  if (NCOL(dataClique)!=1) {
    pcs <- computePCs(dataClique, shrink=shrink, method=method, cliques=cliques, maxPCs=maxPCs)
    colnames(pcs$x) <- paste0(name,colnames(pcs$x))
    names(pcs$sdev) <- paste0(name,names(pcs$sdev))
    colnames(pcs$loading) <- paste0(name,colnames(pcs$loading))
  } else {
    colnames(dataClique) <- paste0(name,"PC1")
    pcs <- list(x=dataClique, sdev=sd(dataClique), loadings=1)
  }

  pcs$dataModule <- t(dataClique)
  pcs$method="pca"
  pcs$namesCov=colnames(pcs$x)
  pcs$omicName=name
  pcs
}

#' Summarize With Directed Sum
#'
#' Given a matrix it summarize the positive and negative in two vectors
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param eventThr the absolute value to threshold an event
#' @param min_prop minimal proportion in classes
#' @param cliques the features organized in cliques. Only use for topology.
#'
#' @return NULL
#'
#' @export
summarizeToNumberOfDirectionalEvents <- function(data, features, name="dCount", eventThr=2,
                                        min_prop=0.1, cliques=NULL) {
  if (is.null(data))
    return(NULL)
  
  genes <- intersect(row.names(data), features)
  if (length(genes)==0)
    return(NULL)
  
  dataClique <- t(data[genes, , drop=F])
  if (ncol(dataClique)==0)
    return(NULL)
  
  posDataClique <- extractPositivePortion(dataClique)
  negDataClique <- extractPositivePortion(dataClique, invert=TRUE)
  
  positive <- apply(posDataClique >= eventThr, 1, sum, na.rm=T)
  negative <- apply(negDataClique >= eventThr, 1, sum, na.rm=T)
  
  collapsed <- data.frame(positive=positive, negative=negative,
                          row.names = names(positive), stringsAsFactors = F)
  colnames(collapsed) <- paste0(name, c("POS","NEG"))
  
  # min <- ceiling(NCOL(data)*0.01)
  # keep <- colSums(collapsed>0) >=  min | colSums(collapsed>0) <= NROW(dataClique)-min
  keep <- sapply(collapsed, check_minimal_proportion, min_prop=min_prop)
  collapsed = collapsed[, keep, drop=F]
  
  if (NCOL(collapsed) == 0)
    return(NULL)
  
  list(x=collapsed, dataModule=t(dataClique), namesCov=names(collapsed),
       method="directedCount", omicName=name, eventThr=eventThr, min_prop=min_prop)
}

#' @importFrom stats quantile
check_minimal_proportion <- function(x, min_prop=0.1){
  min <- quantile(x, probs=c(min_prop))
  max <- quantile(x, probs=c(1-min_prop))
  if ((min==min(x)) && (max==min(x)))
    return(FALSE)
  
  if ((min==max(x)) && (max==max(x)))
    return(FALSE)
  
  TRUE
}

#' Summarize To Binary Directional Events
#'
#' Given a matrix it summarize the positive and negative in two vectors
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param binaryClassMin the minimum number of event to include the covariate
#' @param eventThr the absolute value to threshold an event
#' @param cliques the features organized in cliques. Only use for topology.
#'
#' @return NULL
#'
#' @export
summarizeToBinaryDirectionalEvents <- function(data, features, name="dirBin", binaryClassMin=10, 
                                        eventThr=2, cliques=NULL) {
  if (is.null(data))
    return(NULL)
  
  genes <- intersect(row.names(data), features)
  if (length(genes)==0)
    return(NULL)
  
  dataClique <- t(data[genes, , drop=F])
  if (ncol(dataClique)==0)
    return(NULL)
  
  posDataClique <- extractPositivePortion(dataClique)
  negDataClique <- extractPositivePortion(dataClique, invert=TRUE)
  
  positive <- apply(posDataClique >= eventThr, 1, any, na.rm=T)
  negative <- apply(negDataClique >= eventThr, 1, any, na.rm=T)
  
  collapsed <- data.frame(positive=positive, negative=negative,
                          row.names = names(positive), stringsAsFactors = F)
  colnames(collapsed) <- paste0(name, c("POS","NEG"))
  
  # if (sum(collapsed) < binaryClassMin | sum(collapsed) > NROW(dataClique)-binaryClassMin)
  #   return(NULL)
  keep <- sapply(collapsed, sum) >= binaryClassMin | sapply(collapsed, sum) <= NROW(dataClique)-binaryClassMin
  collapsed = collapsed[, keep, drop=F]
  
  if (NCOL(collapsed) == 0)
    return(NULL)
  
  list(x=collapsed, dataModule=t(dataClique), namesCov=names(collapsed),
       method="directedBinary", omicName=name, eventThr=eventThr)
}
