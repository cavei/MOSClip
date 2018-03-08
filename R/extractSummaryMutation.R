extractSummaryMutation <- function(omic, n=3) {
  moduleMat <- t(omic$dataModule)
  involved <- head(mostlyMutated(moduleMat), n)
  sigModule <- omic$dataModule[row.names(involved), , drop=F]
  discrete=data.frame(lapply(omic$x, as.numeric), row.names=row.names(omic$x))
  list(sigModule=sigModule, discrete=discrete, subset=involved, covsConsidered=omic$namesCov)
}

mostlyMutated <- function(moduleMat) {
  priority <- colSums(moduleMat, na.rm=T)
  priority <- data.frame(row.names = names(priority), patientsMutated=priority)
  priority[order(priority$patientsMutated, decreasing = T),, drop=F]
}

