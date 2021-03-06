% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plottingMO.R
\name{plotModuleReport}
\alias{plotModuleReport}
\title{Summarize and plot pathways' info from a MultiOmicsModule (MOM) object}
\usage{
plotModuleReport(pathwayObj, MOcolors = NULL, priority_to = NULL, ...)
}
\arguments{
\item{pathwayObj}{MultiOmicsModule of pathway object}

\item{MOcolors}{character vector with the omic colors.
The colors should be among the colors in \code{showMOSpalette()}}

\item{priority_to}{a vector with the covariates (omic name) that should go first}

\item{...}{additional argument to be passed to pheatmap}
}
\description{
Given a MOM, it plots the table.
}
