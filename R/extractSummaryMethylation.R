extractSummaryMethylation <-function(omic, n=3) {
  moduleMat <- t(omic$dataModule)
  classes <- omic$x[,1]
  KMsigMat <- KWtest(moduleMat, classes)
  g <- lapply(names(omic$cls), function(gene) {
    cbind(gene, omic$cls[[gene]])
  })
  g <- do.call(rbind, g)
  genes <- g[,1]
  names(genes) <- g[,2]

  involved <- head(KMsigMat, n)
  sigModule <- omic$dataModule[row.names(involved), , drop=F]
  topGenes <- genes[row.names(involved)]
  topGenes <- tapply(names(topGenes), topGenes, paste, collapse=";")

  list(sigModule=sigModule, discrete=omic$x, subset=data.frame(row.names=names(topGenes), metClust=topGenes), pvalues=involved,
       covsConsidered=omic$namesCov)
}

KWtest <- function(moduleMat, classes) {
  kwTest <- apply(moduleMat, 2, function(gene) {
    kruskal.test(x=gene, g=classes)[c("p.value", "statistic")]
  })
  res <- do.call(rbind, kwTest)
  res[order(as.numeric(res[,"p.value"])),]
}
