% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3_DriveReferenceSetups.R
\name{MakeReference_DaisyDrive}
\alias{MakeReference_DaisyDrive}
\title{Daisy Drive Offspring Reference}
\usage{
MakeReference_DaisyDrive(H, R, S, d, eta=NULL, phi=NULL,
omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
}
\arguments{
\item{H}{Vector of homing rates for each drive piece}

\item{R}{Vector of deleterious allele generation rates}

\item{S}{Vector of neutral allele generation rates}

\item{d}{Vector of background mutation rates at each locus}

\item{eta}{Named vector of mating fitness}

\item{phi}{Named vector of sex ratios at emerence}

\item{omega}{Named vector of adult mortality increase}

\item{xiF}{Named vector of female pupatory success}

\item{xiM}{Named vector of male pupatory success}

\item{s}{Named vector of genotype-dependent fertility reduction}
}
\value{
list of homing, cutting, and mendelian genotypes and rates
}
\description{
Create a list specifying the offspring probability distribution.
}
\details{
This function creates a reference list for \code{\link{DaisyOffspring}}.
The number of drive elements is specified by the length of H. The final homing
rate doesn't matter, as the last piece of the drive has nothing to drive, but
it must be there for compatibility. R,S, and d must be the same length as H,
but will generally be the same number replicated that many times. It is assumed
that S is R/3, but this can be varied. This function is similar to \code{\link{MakeReference_Multiplex_mLoci}}
and \code{\link{MakeReference_Multiplex_oLocus}}
}
\examples{
H <- c(0.9,0.4,0) # This drive has 3 pieces
R <- c(0.001, 0.002, 0.003)
S <- c(0.0003, 0.0006, 0.001)
d <- c(0.0001, 0.0001, 0.0001)

MakeReference_DaisyDrive(H,R,S,d)

}