



#' Run one run of mPlex
#' 
#' R interface to C++ simulation code.
#' 
#' @param seed An integer seed for random number generator
#' @param networkParameters A list of simulation parameters
#' @param reproductionReference A list of reproduction and genotype specific parameters
#' @param migrationMale A matrix specifying male migration rates
#' @param migrationFemale A matrix specifying female migration rates
#' @param migrationBatch A list specifing batch migration probabilities, rates, and sex ratios
#' @param reproductionType String specifying type of reproduction model
#' @param verbose Boolean, be chatty?
#' 
#' @export
mPlex_oneRun <- function(seed, networkParameters, reproductionReference,
                         migrationMale, migrationFemale, migrationBatch,
                         reproductionType, verbose){
  
  
  
  
  mPlexCpp:::run_mPlex_Cpp(seed,
                           networkParameters,
                           reproductionReference,
                           migrationMale,
                           migrationFemale,
                           migrationBatch,
                           reproductionType,
                           verbose)
  
  

  
  
}