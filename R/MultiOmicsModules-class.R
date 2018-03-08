setClassUnion("characterOrNULL", c("character", "NULL"))

MultiOmicsModules <- setClass("MultiOmicsModules", package = "MOSClip",
                        slots = c(alphas  = "numeric",
                                  zlists  = "list",
                                  coxObjs = "list",
                                  modulesView  = "list",
                                  modules     = "list",
                                  modulesData = "list",
                                  formulas = "list",
                                  graphNEL = "graphNEL",
                                  title = "characterOrNULL"),
                        contains = "list"
)

setMethod("show",
          signature = "MultiOmicsModules",
          definition = function(object) {
            sthis <- seq_len(min(length(object@alphas), 3))
            sthis <- order(object@alphas)[sthis]

            sigCliquesIdx = which(object@alphas <= 0.05)

            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }

            for (i in sthis) {
              cat(paste0("Module ",i, ": pvalue ", object@alphas[i], "\n"))
              covs <- names(which(object@zlists[[i]] <=0.05))
              if (length(covs)!=0)
                cat("The following covariates are implicated:\n",paste(covs, collapse=", "),"\n")
              cat("Module is composed by the followings:\n")
              cat(paste(object@modules[[i]], collapse=", "))
              cat("\n-+-\n")
            }

            if (length(sthis) < length(sigCliquesIdx)) {
              cat(paste0("There are other ", length(sigCliquesIdx)-length(sthis), " cliques with pvalue <= 0.05"))
            }

            invisible(NULL)
          })


