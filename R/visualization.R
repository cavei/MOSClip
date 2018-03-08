formatModuleReport <- function(smObj){
  alphas <- smObj@alphas
  z  <- smObj@zlists
  idxs <- order(alphas)

  zcols <- unique(unlist(lapply(z, function(x){
    names(x)
  })))

  colDescription <- do.call(rbind, lapply(idxs, function(i){
    additionalCols <- rep(NA, length(zcols))
    names(additionalCols) <- zcols
    additionalCols[names(z[[i]])] <- z[[i]]
    additionalCols
  }))

  cbind(row.names=idxs, pvalue=alphas[idxs], data.frame(colDescription))
}

multiPathwayReport <- function(multiPathway, p.adjust.method = p.adjust.methods){
  if (!is(multiPathway,"list"))
    stop("A list of pathway results are expected.")

  if (all(p.adjust.method==p.adjust.methods)) {
    p.adjust.method = "BH"
  }

  pvalues <- sapply(multiPathway, function(p) {as.numeric(p@pvalue)})
  p.adj <- p.adjust(pvalues, method=p.adjust.method)

  zs <- sort(unique(unlist(lapply(multiPathway, function(p) {
    names(p@zlist)
  }))))

  zMat <- do.call(rbind, lapply(multiPathway, function(p) {
    fixedCols=rep(NA, length(zs))
    names(fixedCols)<-zs
    fixedCols[names(p@zlist)] <- p@zlist
    fixedCols
  }))

  ord <- order(pvalues)
  cbind(row.names=names(pvalues)[ord], pvalue=pvalues[ord], p.adj=p.adj[ord], data.frame(zMat[ord, , drop=F]))
}

multiPathwayModuleReport <- function(multiPathwayModuleList, top=25) {
  if (!is(multiPathwayModuleList,"list"))
    stop("A list of pathway results are expected.")

  multiMatrixRes <- lapply(multiPathwayModuleList, function(p) {
    summary <- formatModuleReport(p)
    data.frame(pathway=p@title, module=row.names(summary), summary, row.names=NULL, stringsAsFactors = F)
  })

  resDF <- mergeAll(multiMatrixRes)

  if (top=="all") {
    top=NROW(resDF)
  }
  head(resDF[order(resDF$pvalue), ], min(top, NROW(resDF)))
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
  data.frame(pathway=matrix[, "pathway"], module=matrix[, "module"], pvalue=matrix[, "pvalue"],
             apply(numericMat, 2, as.numeric), stringsAsFactors = F)
}

