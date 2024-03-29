% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_PatchClass.R
\docType{class}
\name{Patch}
\alias{Patch}
\title{Patch Class Definition}
\format{An \code{\link{R6Class}} generator object}
\usage{
Patch
}
\description{
A Patch represents a single, well-mixed population of \code{\link{Mosquito}}
}
\section{\strong{Constructor}}{

\itemize{
\item patchID: Integer ID of this patch
\item simTime: Maximum time of simulation
\item eggs_t0: List of Mosquito objects as the initial egg population, \eqn{L_{eq}}
\item larva_t0: List of Mosquito objects as the initial larval population, \eqn{L_{eq}}
\item pupa_t0: List of Mosquito objects as the initial pupae population, \eqn{L_{eq}}
\item adult_male_t0: List of Mosquito objects as the initial adult male population, \eqn{Ad_{eq}}
\item unmated_female_t0: List of Mosquito objects as the initial adult female population, \eqn{Ad_{eq}}
\item maleReleases: Male release schedule for this patch, see \code{\link{Release_basicRepeatedReleases}}
\item femaleReleases: Female release schedule for this patch, see \code{\link{Release_basicRepeatedReleases}}
\item larvaeReleases: Larvae release schedule for this patch, see \code{\link{Release_basicRepeatedReleases}}
}
}

\section{\strong{Methods}}{

\itemize{
\item get_patchID: Return ID of current patch
\item get_egg: Return current egg population
\item get_larva: Return current larva population
\item get_pupa: Return current pupa population
\item get_adult_male: Return current adult male population
\item get_adult_female: Return current adult female population
\item get_unmated_female: Return current adult, unmated female population
\item get_maleMigration: Return migratory males (nPatch length list of lists)
\item get_femaleMigration: Return migratory females (nPatch length list of lists)
\item set_NetworkPointer: Set a referenced to the enclosing \code{\link{Network}} object
\item get_NetworkPointer: Return a reference to the enclosing \code{\link{Network}} object
\item reset: see \code{\link{reset_Patch}}
\item oneDay_initOutput: see \code{\link{oneDay_initOutput_Patch}}
\item oneDay_writeOutput: see \code{\link{oneDay_writeOutput_Patch}}
\item oneDay_maleReleases: see \code{\link{oneDay_maleReleases_Patch}}
\item oneDay_femaleReleases: see \code{\link{oneDay_femaleReleases_Patch}}
\item oneDay_larvaeReleases: see \code{\link{oneDay_larvaeReleases_Patch}}
\item oneDay_migrationOut: see \code{\link{oneDay_migrationOut_Patch}}
\item oneDay_migrationIn: see \code{\link{oneDay_migrationIn_Patch}}
\item oneDay_PopDynamics: see \code{\link{oneDay_PopDynamics_Patch}}
\item oneDay_Age: see \code{\link{oneDay_Age_Patch}}
\item oneDay_EggDeath: see\code{\link{oneDay_EggDeath_Patch}}
\item oneDay_LarvalDeath: see \code{\link{oneDay_LarvalDeath_Patch}}
\item oneDay_PupaDeath: see \code{\link{oneDay_PupaDeath_Patch}}
\item oneDay_AdultDeath: see \code{\link{oneDay_AdultDeath_Patch}}
\item oneDay_EggMaturation: see \code{\link{oneDay_EggMature_Patch}}
\item oneDay_LarvaMaturation: see \code{\link{oneDay_LarvaMature_Patch}}
\item oneDay_PupaMature: see \code{\link{oneDay_PupaMature_Patch}}
\item oneDay_Mate: see \code{\link{oneDay_Mate_Patch}}
\item oneDay_Reproduction: See \code{\link{oneDay_Reproduction_Patch}}
}
}

\section{\strong{Fields}}{

\itemize{
\item patchID: Integer ID of current patch
\item genericCounter: Generic counter used for indexing
\item eggs_t0: List of Mosquitoes comprising the initial egg population
\item larva_t0: List of Mosquitoes comprising the initial larval population
\item pupa_t0: List of Mosquitoes comprising the initial pupa population
\item adult_male_t0: List of Mosquitoes comprising the initial adult male population
\item adult_female_t0: List of Mosquitoes comprising the initial adult female population
\item unmated_female_t0: List of Mosquitoes comprising the initial unmated female population
\item eggs: List of Mosquitoes comprising the current egg population
\item larva: List of Mosquitoes comprising the current larval population
\item pupa: List of Mosquitoes comprising the current pupa population
\item adult_male: List of Mosquitoes comprising the current adult male population
\item adult_female: List of Mosquitoes comprising the current adult female population
\item unmated_female: List of Mosquitoes comprising the current unmated female population
\item maleMigration: List of outbound male Mosquitoes, length nPatch
\item femaleMigration: List of outbound female Mosquitoes, length nPatch
\item numMigration: Holder for the number of migratory Mosquitoes, integer
\item migrationDist: Holder for the distribution of Mosquitoes over patches, vector length nPatch-1
\item whoMigrate: Holder for which mosquitoes migrate, vector of indices
\item numMigrateWhere: Holder for where and how many mosquitoes migrate, vector of integers
\item maleReleases: Male release schedule and list of Mosquitoes to release
\item femaleReleases: female release schedule and list of Mosquitoes to release
\item larvaeReleases: Larval release schedule and list of Mosquitoes to release
\item NetworkPointer: a reference to enclosing \code{\link{Network}}
\item DenDep: Density-dependent parameter for larvae
\item death: Vector of T/F death at each stage. Slowly grows with population size
\item meanAge: Holder for maturation functions
\item sdAge: Holder for maturation functions
\item ages: Holder for maturation functions
\item matured: Holder for maturation functions
\item numUnweds: Number of unmated females. Holder for mating function
\item numMates: Number of adult males. Holder for mating function
\item mates: Character vector of male genotypes. Holder for mating function
\item offspring: List(Alleles, Probabilities) specifying new offspring genotype and distribution.
\item eggNumber: Integer vector of how many eggs to lay. Holder for reproduction function
\item newEggs: list of new Mosquitoes laid. Holder for reproduction function
}
}

\keyword{R6}
\keyword{class}
