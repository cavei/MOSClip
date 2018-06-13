#' Get available Omics Summarizing Methods
#'
#' Gives a vector of the available methods to summarize omics.
#'
#' @return character vector with the implemented methods.
#'
#' @export
availableOmicMethods <- function() {
  return(c("summarizeToBinaryEvents",
           "summarizeWithPca",
           "summarizeInCluster",
           "summarizeInClusterWithoutDictionary"))
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
summarizeToBinaryEvents <- function(data, features, name="cov",
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
  list(x=collapsed, dataModule=t(dataClique), namesCov=name, method="binary")
}

#' Summarize Using CLuster Analysis
#'
#' Given a matrix it summarize to a 0 or 1
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param cliques the features organized in cliques. Only use for topology.
#'
#' @return NULL
#' @importFrom stats cutree dist hclust
#' @export
summarizeInCluster <- function(data, features, name="clust", cliques=NULL) {
  datamat <- data$met ## modificare in modo che se non c'è un dizionario usi i geni stessi
  dict <- data$dict

  if (is.null(datamat))
    return(NULL)
  genes <- intersect(names(dict), features)

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
summarizeInClusterWithoutDictionary <- function(datamat, features, name="clust", cliques=NULL) {
  
  if (is.null(datamat) | (ncol(datamat)==0) | !(is.matrix(datamat)))
    return(NULL)
  genes <- intersect(rownames(datamat), features)
  
  if (length(genes)==0)
    return(NULL)
  
  datamatClique <- t(datamat[genes, ,drop=FALSE])
  
  used <- colnames(datamatClique)
  names(used) <- colnames(datamatClique)
  
  hc <- hclust(dist(datamatClique, method = "euclidean"), method="ward.D2")
  # clusters <- kmeans(datamatClique, centers=2) # TO DO: add kmeans
  
  if (ncol(datamatClique)<4){
    covs <- data.frame(factor(cutree(hc, k = 2)), stringsAsFactors = T)
    names(covs) <- paste0(name,"_2k")
  } else {
    covs <- data.frame(factor(cutree(hc, k = 3)), stringsAsFactors = T)
    names(covs) <- paste0(name,"_3k")
  }
  
  list(x=covs, dataModule=t(datamatClique), namesCov=names(covs), cls=used, method="cluster")
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
summarizeWithPca <- function(data, features, name="exprs", shrink=FALSE, method="regular", cliques=NULL, maxPCs=3, loadThr=0.6) {
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
  } else {
    colnames(dataClique) <- "PC1"
    pcs <- list(x=dataClique, sdev=sd(dataClique), loadings=1)
  }

  pcs$dataModule <- t(dataClique)
  pcs$method="pca"
  pcs$namesCov=colnames(pcs$x)
  pcs
}
