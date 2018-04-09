pathwayRadialSummary <- function() {
  wholeSummary <- multiPathwayReport(multiOmicsFull)
  wholeSummary.sig <-wholeSummary[wholeSummary$pvalue<0.06,]

  pathway2id <- mapReactomeIDfromGraphite(reactome) # codici
  pathwayDict <- pathway2id$pname
  names(pathwayDict) <- pathway2id$id

  wholePath2id <- mapReactomeIDfromGraphite(reactome, row.names(wholeSummary.sig))

  wfathers <- lapply(wholePath2id$id, getPathFathers, df.g, ord=length(wholePath2id$id)+1)
  names(wfathers) <- wholePath2id$id
  head(wfathers)

  wpathway2father <- lapply(wfathers, function(x) {
    unlist(unname(pathwayDict[x]))
  })

  if(identical(wholePath2id$id, names(wpathway2father))) {
    out <- lapply(seq_len(nrow(wholeSummary.sig)), function(idx) {
      suppressWarnings(cbind(wholeSummary.sig[idx, , drop=F], slim=wpathway2father[[idx]]))
    })
    moduleMat <- do.call(rbind, out)
  }

  row.names(moduleMat)<-paste(moduleMat$slim,1:nrow(moduleMat),sep="-")

  met<-apply(moduleMat[,c(2,3)],1,min,na.rm=T)
  expr<-apply(moduleMat[,c(5:7)],1,min,na.rm=T)
  names(expr)<-row.names(moduleMat)
  names(met)<-row.names(moduleMat)

  isSig.expr<-rep(0,nrow(moduleMat))
  isSig.expr[expr<0.05]<-1

  isSig.met<-rep(0,nrow(moduleMat))
  isSig.met[met<0.05]<-1

  isSig.mut<-rep(0,nrow(moduleMat))
  isSig.mut[moduleMat$mutationTRUE<0.05]<-1

  names(isSig.expr)<-row.names(moduleMat)
  names(isSig.met)<-row.names(moduleMat)
  names(isSig.mut)<-row.names(moduleMat)

  multiOmicIntersection <- MOSClip:::createIntersection(list(mut=isSig.mut, met=isSig.met, expr=isSig.expr), nrow(wholeSummary))
  multiOmicIntersectSummary <- summary(multiOmicIntersection)
  classes <- multiOmicIntersectSummary$Table$Elements  # tmp2
  names(classes)<-multiOmicIntersectSummary$Table$Intersections
  omicsIsSig <- data.frame(mut=isSig.mut,met=isSig.met,expr=isSig.expr)

  frequences<-c()
  outClass<-c()
  myNames<-c()

  splitClasses <- c()
  for(i in seq_along(classes)){
    class<-classes[i]
    subClasses <- unlist(strsplit(unlist(class),","))
    spClasses <-unlist(lapply(subClasses, function(subClass) strsplit(subClass,"-")[[1]][1]))
    spClasses <- gsub("^ ", "",spClasses)
    # splitClasses <- c(splitClasses, spClasses)
    frequences<-c(frequences,table(spClasses))
    outClass<-c(outClass,names(table(spClasses)))
    myNames<-c(myNames,rep(names(classes[i]),length(table(spClasses))))
  }

  df<-data.frame(categorie=outClass,freq=frequences, inter=myNames)

  df.multi <- df[grep("&", df$inter), , drop=F]
  library(ggplot2)
  ggplot(df.multi, aes(y = freq, x = reorder(categorie, inter),
                       group = inter, colour = inter)) + coord_polar() +
    geom_point(stat='identity') +
    geom_polygon(fill=NA)+
    geom_path() +
    theme(axis.text.x=element_text(size=3, angle=45))+
    labs(x = NULL)+
    theme_bw()


  df.single <- df[-grep("&",df$inter),]
  ggplot(df.single, aes(y = freq, x = reorder(categorie, inter),
                        group = inter, colour = inter)) + coord_polar() +
    geom_point(stat='identity') +
    geom_polygon(fill=NA)+
    geom_path() +
    theme(axis.text.x=element_text(size=3, angle=45))+
    labs(x = NULL)+
    theme_bw()



  ggplot(df, aes(y = freq, x = reorder(categorie, inter),
                 group = inter, colour = inter)) + coord_polar() +
    geom_point(stat='identity') +
    geom_polygon(fill=NA)+
    geom_path() +
    theme(axis.text.x=element_text(size=3, angle=45))+
    labs(x = NULL)+
    theme_bw()


}

createReactome <- function() {
  # load Hierarchy
  pathHierarchy <- read.table("ReactomePathwaysRelation-HSA.txt", stringsAsFactors = F)
  df.g <- graph.data.frame(d = pathHierarchy, directed = TRUE)

  pathway2id <- mapReactomeIDfromGraphite(reactome) # codici
  pathwayDict <- pathway2id$pname
  names(pathwayDict) <- pathway2id$id

  wholePath2id <- mapReactomeIDfromGraphite(reactome, wpath$row.names)

  wfathers <- lapply(wholePath2id$id, getPathFathers, df.g, ord=length(wholePath2id$id)+1)
  names(wfathers) <- wholePath2id$id
  head(wfathers)

  wpathway2father <- lapply(wfathers, function(x) {
    unlist(unname(pathwayDict[x]))
  })
}
