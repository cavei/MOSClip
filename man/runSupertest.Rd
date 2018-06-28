% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supertest.R
\name{runSupertest}
\alias{runSupertest}
\title{Performs statistical tests of the pathway intersections among omics.}
\usage{
runSupertest(multiPathwayReportData, pvalueThr = 0.05, plot = c("circular",
  "landscape", "noplot"), sort.by = c("set", "size", "degree", "p-value"),
  excludeColumns = NULL, ...)
}
\arguments{
\item{multiPathwayReportData}{data.frame, the output of 
the \code{\link{multiPathwayReport}} or \code{\link{multiPathwayModuleReport}} functions.}

\item{pvalueThr}{integer indicating the p-value cut-off.}

\item{plot}{character indicating the layout for plotting. 
It is one of \code{circular}, \code{landscape} or \code{no}. 
By default, plot="circular", if plot="no" no plot will be provided.}

\item{sort.by}{character indicating how to sort the 
intersections in the plot. It is one of "set" (by omics), "size" 
(by intersection size), "degree" (by number of intersected omics), 
and "p-value".}

\item{excludeColumns}{a vector of characters listing the columns of 
\code{multiPathwayReportData} object to be excluded by the analysis. 
In the case \code{multiPathwayReportData} derives from \code{\link{multiPathwayModuleReport}} 
you should set \code{excludeColumns = c("pathway","module")}.}
}
\value{
a data.frame containing all the numeric information of the plot included
 the pathways shared by different omics.
}
\description{
This function calculates intersection sizes among the pathway sets 
significative in each omic and performs statistical intersection test.
}
\details{
This function calculates intersection sizes between multiple set of pathways or modules 
and performs statistical test of the intersections using the total amout of 
analyzed pathways or modules as background. The super exact test of this function 
was described in Wang et al 2015.
}
\references{
Minghui Wang, Yongzhong Zhao, and Bin Zhang (2015). 
Efficient Test and Visualization of Multi-Set Intersections. 
Scientific Reports 5: 16923.
}