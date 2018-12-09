context("Test Supertest accessory functions")
library(houseOfClipUtility)
library(graphite)
pathHierarchy <- houseOfClipUtility::downloadPathwayRelationFromReactome()
pathHierarchyGraph <- igraph::graph.data.frame(d = pathHierarchy, directed = TRUE)

pathwayName = c("Circadian Clock", "Signaling Pathways")
reactome <- pathways("hsapiens", "reactome")[pathwayName]

omicsClasses2pathways <- list(exp=pathwayName, 'cnv;exp'="Not a pathway")
test_that("generateWarn", {
  expect_warning(omicsClasses2fathers <- lapply(omicsClasses2pathways, MOSClip::annotePathwayToFather, 
                                 graphiteDB=reactome, hierarchy=pathHierarchyGraph))
  expect_warning(correspondence <- lapply(names(omicsClasses2pathways), match_pathway_to_fathers, omicsClasses2pathways=omicsClasses2pathways,
                           omicsClasses2fathers=omicsClasses2fathers))
  
})



