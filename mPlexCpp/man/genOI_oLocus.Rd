% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_OutputAnalysis.R
\name{genOI_oLocus}
\alias{genOI_oLocus}
\title{Generate genotypes of interest for mPlex-oLocus}
\usage{
genOI_oLocus(outputFile, alleles, collapse)
}
\arguments{
\item{outputFile}{Name of file to output. Must end in .csv}

\item{alleles}{A list of lists that contain the genotypes of interest at each locus. Default is all genotypes}

\item{collapse}{A list of lists containing TRUE/FALSE for each locus. If TRUE, the genotypes of interest at that locus are collapsed and the output is the sum of all of them.}
}
\description{
This function generates a *.csv of genotypes of interest and the grouping scheme, ie,
if one locus isn't interesting and should be considered the same. This file will be used
for data analysis.
}
