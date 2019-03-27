#' Summarize pathways' info from a list of MultiOmicsPathway (MOP)
#'
#' Given the list of MOPs, it create the table.
#'
#' @param multiPathwayList MultiOmicsPathway pathway object list
#' @param priority_to a vector with the covariates (omic name) that should go first
#'
#' @return a data.frame with overal pvalue of the coxph, followed by covariates zs.
#'
#' @export
multiPathwayReport <- function(multiPathwayList, priority_to=NULL){
  if (!is(multiPathwayList,"list"))
    stop("A list of pathway results are expected.")

  pvalues <- sapply(multiPathwayList, function(p) {as.numeric(p@pvalue)})

  zs <- sort(unique(unlist(lapply(multiPathwayList, function(p) {
    names(p@zlist)
  }))))

  zMat <- do.call(rbind, lapply(multiPathwayList, function(p) {
    fixedCols=rep(NA, length(zs))
    names(fixedCols)<-zs
    fixedCols[names(p@zlist)] <- p@zlist
    fixedCols
  }))

  ord <- order(pvalues)
  df <- cbind(row.names=names(pvalues)[ord], pvalue=pvalues[ord], data.frame(zMat[ord, , drop=F]))
  order_by_covariates(df, 1, priority_to)
}

#' Provide a data.frame of pathways module test results from list of Multi Omics Module (MOM) objects
#'
#' Given the list of MOMs, it creates the table.
#'
#' @param multiPathwayModuleList a list of Multi Omics Modules resulting from a multi-omics module test.
#' @param priority_to a vector with the covariates (omic name) that should go first
#'
#' @return a data.frame, modules in rows, overall and covariate pvalues of the test in columns.
#' #'
#' @export
multiPathwayModuleReport <- function(multiPathwayModuleList, priority_to=NULL) {
  if (!is(multiPathwayModuleList,"list"))
    stop("A list of pathway results are expected.")

  multiMatrixRes <- lapply(multiPathwayModuleList, function(p) {
    summary <- formatModuleReport(p)
    data.frame(pathway=p@title, module=row.names(summary), summary,
               row.names=NULL, stringsAsFactors = F)
  })
  resDF <- mergeAll(multiMatrixRes)
  resDF <- resDF[order(resDF$pvalue), ]
  rownames(resDF) <- apply(resDF,1, function(r) paste(r["pathway"],r["module"],sep="."))
  
  resDF <- order_by_covariates(resDF, 3, priority_to)
  
  return(resDF)
}

formatModuleReport <- function(smObj){
  alphas <- smObj@alphas
  z  <- smObj@zlists
  idxs <- order(alphas)

  zcols <- sort(unique(unlist(lapply(z, function(x){
    names(x)
  }))))

  colDescription <- do.call(rbind, lapply(idxs, function(i){
    additionalCols <- rep(NA, length(zcols))
    names(additionalCols) <- zcols
    additionalCols[names(z[[i]])] <- z[[i]]
    additionalCols
  }))

  cbind(row.names=idxs, pvalue=alphas[idxs], data.frame(colDescription))
}

mergeAll <- function(list) {
  allColumnsNames <- sort(unique(unlist(lapply(list, function(o) {
    colnames(o)
  }))))

  matrix <- do.call(rbind,
                    lapply(list, function(o) {
                      fixedCols=matrix(NA, NROW(o),length(allColumnsNames))
                      colnames(fixedCols)<-allColumnsNames
                      for (col in colnames(o)) {
                        fixedCols[, col] <- o[,col]
                      }
                      fixedCols
                    }))
  removeCols <- match(c("pathway", "module","pvalue"), colnames(matrix))
  numericMat <- matrix[,-removeCols, drop=F]
  data.frame(pathway=matrix[, "pathway"], module=matrix[, "module"],
             pvalue=as.numeric(matrix[, "pvalue"]),
             apply(numericMat, 2, as.numeric), stringsAsFactors = F)
}

order_by_covariates <- function(dataF, skip_first_cols, priority_to=NULL) {
  # priority_to <- c("exp","cnv", "met", "mut")
  if (is.null(priority_to))
    return(dataF)
  
  covariates <- dataF[, seq_len(ncol(dataF)-skip_first_cols)+skip_first_cols, drop=F]
  fixed <- dataF[, seq_len(skip_first_cols), drop=F]
  omics_cov <- guessOmics(colnames(covariates))
  to_sort <- unique(omics_cov) %in% priority_to
  priority_to <- c(priority_to, unique(omics_cov)[!to_sort])
  cov_factor <- factor(omics_cov, levels=priority_to)
  cov_factor <- droplevels(cov_factor)
  covariates <- covariates[, order(cov_factor), drop=F]
  cbind(fixed, covariates)
}
