% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_PatchSimulation.R
\name{oneDay_Mate_Patch}
\alias{oneDay_Mate_Patch}
\title{Daily Mating}
\usage{
oneDay_Mate_Patch()
}
\description{
Freshly matured pupa that become females and female releases exist as unmated
females. This function mates unmated_females with any member of the current male
population, with the possibility that males can mate multiple times per mating.
After mating, unmated_females are put into the general adult_female population.
}
