anlayze <- function() {
  # old wpath <- read.table(file="table-whole-pathways-pvalues.txt",sep="\t",header=T)
  module <- moduleSummary # old module
  wpath <- data.frame(sp, stringsAsFactors = F)

  module.sig<-module[module$pvalue<0.05,]

  met<-apply(module.sig[,c(4,5,6)],1 , function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    min(x, na.rm=T)
  })

  expr<-apply(module.sig[,c(8:10)],1, function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    min(x, na.rm=T)
  })

  names(expr)<-row.names(module.sig)
  names(met)<-row.names(module.sig)

  isSig.expr <- rep(0,nrow(module.sig))
  isSig.expr[expr<0.05]<-1

  isSig.met<-rep(0,nrow(module.sig))
  isSig.met[met<0.05]<-1

  isSig.mut<-rep(0,nrow(module.sig))
  isSig.mut[module.sig$mutationTRUE<0.05]<-1

  names(isSig.expr)<-row.names(module.sig)
  names(isSig.met)<-row.names(module.sig)
  names(isSig.mut)<-row.names(module.sig)

  multiOmicIntersection <- MOSClip:::createIntersection(list(mut=isSig.mut, met=isSig.met, expr=isSig.expr), nrow(module))

  omicsIsSig <- data.frame(mut=isSig.mut,met=isSig.met,expr=isSig.expr)

  ## Create reactome DB Not needed
  # require(reactome.db)
  # xx <- as.list(reactomePATHNAME2ID)
  ##

  pathway2id <- mapReactomeIDfromGraphite(reactome) # codici
  pathwayDict <- pathway2id$pname
  names(pathwayDict) <- pathway2id$id

  modules2id <- mapReactomeIDfromGraphite(reactome, module.sig$pathway)
  myIdDict <- modules2id$pname
  names(myIdDict) <- modules2id$id

  pathHierarchy <- read.table("ReactomePathwaysRelation-HSA.txt", stringsAsFactors = F)
  df.g <- graph.data.frame(d = pathHierarchy, directed = TRUE)
  # roots <- findRoots(pathHierarchy)

  # getPathFathers("R-HSA-110373", df.g, ord=2500, plot=T)
  r <- getPathFathers("R-HSA-450294", df.g, ord=3000, plot=T)

  fathers <- lapply(modules2id$id, getPathFathers, df.g, ord=10000)
  names(fathers) <- modules2id$id
  head(fathers)

  module2father <- lapply(fathers, function(x) {
    unlist(unname(pathwayDict[x]))
  })

  if(identical(modules2id$id, names(module2father))) {
    out <- lapply(seq_len(nrow(module.sig)), function(idx) {
      suppressWarnings(cbind(module.sig[idx, , drop=F], slim=module2father[[idx]]))
    })
    moduleMat <- do.call(rbind, out)
  }

  row.names(moduleMat)<-paste(moduleMat$slim,1:nrow(moduleMat),sep="-")
  met<-apply(moduleMat[,c(4,5,6)],1,min,na.rm=T)
  expr<-apply(moduleMat[,c(8:10)],1,min,na.rm=T)
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

  multiOmicIntersection <- MOSClip:::createIntersection(list(mut=isSig.mut, met=isSig.met, expr=isSig.expr), nrow(module))
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
    geom_point(, stat='identity') +
    geom_polygon(fill=NA)+
    geom_path() +
    theme(axis.text.x=element_text(size=3, angle=45))+
    labs(x = NULL)+
    theme_bw()


  df.single <- df[-grep("&",df$inter),]
  ggplot(df.single, aes(y = freq, x = reorder(categorie, inter),
                  group = inter, colour = inter)) + coord_polar() +
    geom_point(, stat='identity') +
    geom_polygon(fill=NA)+
    geom_path() +
    theme(axis.text.x=element_text(size=3, angle=45))+
    labs(x = NULL)+
    theme_bw()



  ggplot(df, aes(y = freq, x = reorder(categorie, inter),
                 group = inter, colour = inter)) + coord_polar() +
    geom_point(, stat='identity') +
    geom_polygon(fill=NA)+
    geom_path() +
    theme(axis.text.x=element_text(size=3, angle=45))+
    labs(x = NULL)+
    theme_bw()

  ############### mappo le categorie sui risultati dei pathway
  wholePath2id <- mapReactomeIDfromGraphite(reactome, wpath$row.names)
  wfathers <- lapply(wholePath2id$id, getPathFathers, df.g, ord=length(wholePath2id$id)+1)
  names(wfathers) <- wholePath2id$id
  head(wfathers)

  wpathway2father <- lapply(wfathers, function(x) {
    unlist(unname(pathwayDict[x]))
  })

  tmp<-table(unlist(wpathway2father))
  # names(tmp)<-names(table(unlist(wpathway2father)))
  tmp<-data.frame(tmp, stringsAsFactors = F)
  tmp<-tmp[tmp$Freq>1,]

  ggplot(tmp) +
    geom_bar(aes(x = Var1, y = Freq),
             stat="identity", position = "dodge", width = 0.7) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size=16),
      axis.text.x = element_text(size=14, angle=45, hjust=1, vjust=1),
    )


}

