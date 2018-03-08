plotClique <- function(r, clique=NULL, entrez=NULL, sortBy=NULL) {
  if(is.null(clique))
    return(getTopLoadGenes(r))

  gns <- r@cliques[[clique]]

  # portions <- extracFromAddiction(gns, addiction)
  # breaks = c(0,0.4,0.8,1)
  # colors = colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(3)
  # pheatmap(portions$metPromo, cluster_rows = F, breaks=breaks, color=colors,
  #          border_color = "black")

  full <- r@coxObjs[[clique]]
  gExp <- r@cliquesExpr[[clique]]
  if (!is.null(entrez))
    gExp <- gExp[, order(gExp[entrez,])]

  if (!is.null(full$metPromo))
    full$metPromo <- ifelse(full$metPromo, 1, 0)
  if (!is.null(full$hyper))
    full$hyper <- ifelse(full$hyper, 1, 0)
  if (!is.null(full$hypo))
    full$hypo <- ifelse(full$hypo, 1, 0)
  if (!is.null(full$mut))
    full$mut <- ifelse(full$mut, 1, 0)

  if(!is.null(sortBy))
    gExp <- gExp[, row.names(full)[order(full[,sortBy])], drop=F]

  pheatmap(log2(gExp+1), scale="row", cluster_cols=FALSE, annotation_col = full)

  # image(z = matrix((objPCA2[colnames(gExp),2]),ncol=1), yaxt = "n",
  #       col= colorRampPalette(rev(brewer.pal(n = 7, name ="YlOrRd")))(70))
  # image(z = matrix(objPCA2$metPromo[colnames(gExp)], ncol=1), yaxt="n", col=c("white", "grey"))
  # image(z = matrix(objPCA2$status(colnames(gExp)), ncol=1), yaxt="n", col=c("white", "black"))
}

plotCliquePCA <- function(res, clique, covariate="status", pcs=c("PC1","PC2")) {
  data <- res@coxObjs[[clique]]
  if (!is.null(data$mut))
    data$mut <- ifelse(data$mut,1,0)

  ggplot(data=data, aes_string(x=pcs[1],y=pcs[2], color=covariate)) +geom_point(size=3)
}

plotKM <- function(res, clique, covariate="PC1", createDiscrete=TRUE) {

  res <- reactAg[[p]]
  clique=4
  covariate="mut"
  createDiscrete=F
  obj <- res@coxObjs[[clique]]

  if (createDiscrete) {
    require(survminer)
    sc <- surv_cutpoint(obj, time="days", event="status", variables = covariate)
    obj <- surv_categorize(sc)
  }

  fm <- paste0("Surv(days, status==1) ~ ", covariate)
  kmFit <- survfit(as.formula(fm), data=obj)

  ggsurv <- ggsurvplot(kmFit, data=NULL, palette = "jco",
                       pval = TRUE, pval.coord = c(500, 0.4), risk.table = TRUE)

  # cs = sort(unique(obj[, covariate]))
  # kmDif <- survdiff(as.formula(fm), data=obj)
  # colors <- c('#990000','#3399FF')

  # plot(kmFit, col=colors, xlab="Time (years)", ylab="Probability", xscale=365.25, mark.time=T, cex.lab=0.7)
  # legend("topright", legend=cs, fill=c(colors))

  diff <- ggarrange(ggsurv$plot, ggsurv$table, heights = c(2, 0.7), ncol = 1, nrow = 2, align = "v")
  diff
}
