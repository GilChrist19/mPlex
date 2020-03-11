###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
# RUN WRAPPER
###############################################################################

#' Run mPlex
#' 
#' R interface to C++ simulation code.
#' 
#' @param seed An integer seed for random number generator. Default is 1
#' @param numReps An integer specifying the number of repetitions. Default is 1
#' @param numThreads An integer specifying the number of threads to parallelize over. Default is 1
#' @param networkParameters A list of simulation parameters
#' @param reproductionReference A list of reproduction and genotype specific parameters
#' @param initAlleles A list of length 1 or numPatches specifying the initial genotype distributions
#' @param patchReleases A list of releases
#' @param migrationMale A matrix specifying male migration rates
#' @param migrationFemale A matrix specifying female migration rates
#' @param migrationBatch A list specifing batch migration probabilities, rates, and sex ratios
#' @param reproductionType String specifying type of reproduction model
#' @param outputDirectory String folder name to write output, default is current directory
#' @param verbose Boolean, be chatty? Default is FALSE
#' 
#' @export
runMPlex <- function(seed = 1, numReps = 1, numThreads = 1,
                     networkParameters, reproductionReference,
                     initAlleles = NULL,
                     patchReleases, migrationMale, migrationFemale,
                     migrationBatch, reproductionType,
                     outputDirectory = "./", verbose = FALSE){
  
  ##########
  # Directory Check
  ##########
  outputDirectory = path.expand(outputDirectory)
  if(!dir.exists(outputDirectory)) stop("Output directory doesn't exist.")
  
  ##########
  # Migration Check
  ##########
  # must be square, must be the same size, size must be equal to the number of patches
  # dim gets rowxcol dimensions, as.matrix protects for single patch run
  migTest <- c(dim(as.matrix(migrationMale)), dim(as.matrix(migrationFemale)))
  if(!all(migTest == netPar$nPatch)) {
    stop("migrationMale or migrationFemale are not properly specified. \nPlease ensure that both matrices are the same size, and the dimensions are equal to the number of patches in the simulation.")
  }
  
  ##########
  # Batch Migration Check
  ##########
  # this is taken care of in basicBatchMigration()
  
  
  ##########
  # Reproduction Type Check
  ##########
  if(!(reproductionType %in% c("DaisyDrive","mPlex_oLocus","mPlex_mLoci","Family"))){
    stop("reproductionType must match one of these choices:\n DaisyDrive\n mPlex_oLocus\n mPlex_mLoci\n Family")
  }
  
  
  ##########
  # reproductionReference Check
  ##########
  # this is taken care of in the MakeReference_* functions
  
  
  ##########
  # Network Parameters Check
  ##########
  # this is taken care of in NetworkParameters()
  
  
  ##########
  # initAlleles Check
  ##########
  # only do this if not family inheritance
  if(!(reproductionType %in% c("Family"))){
    # check and set length
    if(length(initAlleles) == 1){
      initAlleles <- replicate(n = networkParameters$nPatch,
                               expr = unlist(x = initAlleles, recursive = FALSE),
                               simplify = FALSE)
    } else if(length(initAlleles) != networkParameters$nPatch){
      stop("length of initAlleles list must equal the number of patches")
    }
    
    
    # check that the genotypes are the proper length for inheritance
    # options for reproduction reference are: mendelian, homing, mendelianAlleles, homingAlleles
    genCheck <- all(lengths(initAlleles) == length(x = reproductionReference$mendelianAlleles$female))
    if(!genCheck) stop("The number of initial alleles does not match the number of alleles present in the reproduction reference.")
    
    
    # check that alleles are in the possible alleles
    # 1 - unlist the alleles, just need the allele names, not the structure
    # 2 - pull out the "alleles" unit from every list, it's the first element
    # 3 - see if the alleles from each list are one of the 4 allowed ones
    # 4 - unlist because the logical operator only works on vectors
    # 5 - all given alleles must be in the approved set
    alleleCheck <- all(unlist(lapply(X = lapply(X = unlist(initAlleles, recursive = FALSE),
                                                FUN = "[[", 1),
                                     FUN = '%in%', c("W", "H", "R", "S"))
                              )
                       )
    if(!alleleCheck) stop("Not all initial alleles are in the allowed set.")
    
    
    # check that probs are actual probabilities
    # 1 - unlist the alleles, just need the allele names, not the structure
    # 2 - pull out the "probs" unit from every list, it's the second element
    # 3 - see if the probs from each list are greater than or equal to 0
    # 4 - unlist because the logical operator only works on vectors
    # 5 - all given alleles must be in the approved set
    probCheck <- all(unlist(lapply(X = lapply(X = unlist(initAlleles, recursive = FALSE),
                                              FUN = "[[", 2),
                                   FUN = '>=', 0)
                            )
                     )
    if(!probCheck) stop("Not all initial allele probabilities are >= 0.")
    
  } # end initAlleles check
  
  
  ##########
  # Releases Check
  ##########
  # make sure the releases vector is the same length as we have patches
  if(length(patchReleases) != networkParameters$nPatch){
    stop("length of patchReleases list must equal the number of patches")
  }
  
  # check if they are/aren't all null, 
  #  if no releases, skip this.
  # also, if family inheritance, skip this
  allNull <- all(unlist(x = lapply(X = patchReleases, FUN = function(x){lapply(X = x, FUN = is.null)}),
                     use.names = FALSE) == TRUE)
  if(!allNull && !(reproductionType %in% c("Family"))){
    # check that the genotypes are the proper length for inheritance
    
    # Get all genotypes from the releases
    #  Need this for both checks, so getting it here
    # 1 - unlist from each patch
    # 2 - unlist from male/female/eggs/mated females - but not their mates!
    # 3 - get genotypes and mate genotypes vector from the genotypes, ages, and time
    # 4 - unlist all of the genotypes into a basic character vector
    allGenos <- unlist(x = lapply(X = unlist(x = unlist(x = patchReleases,
                                                        recursive = FALSE),
                                             recursive = FALSE),
                                  FUN = function(x){c(x[["genVec"]],x[["mateVec"]])}),
                       recursive = TRUE, use.names = FALSE)
    
    # get length of every genotype in the releases
    # 1 - get length of each of the strings in the character vector
    # options for reproduction reference are: mendelian, homing, mendelianAlleles, homingAlleles
    #  times 2 because everything is diploid, this will fail if anyting becomes monoploid
    genLength <- nchar(x = allGenos)
    genCheck <- all(genLength == (length(x = reproductionReference$mendelianAlleles$female)*2) )
    if(!genCheck) stop("The genotypes of released individuals do not match the number of alleles present in the reproduction reference.")
    
    
    # check that alleles are in the possible alleles
    alleleCheck <- all(unlist(x = strsplit(x = allGenos, split = '', fixed = TRUE, useBytes = TRUE)
                              ) %in% c("W", "H", "R", "S") )
    if(!alleleCheck) stop("Not all genotypes of released individuals are in the allowed set.")
    
  } # end releases check
  
  
  ##########
  # Run
  ##########
  # pass all down for simulation
  run_mPlex(seed_ = seed, 
            numReps_ = numReps, 
            numThreads_ = numThreads,
            networkParameters_ = networkParameters,
            reproductionReference_ = reproductionReference, 
            initAlleles_ = initAlleles,
            patchReleases_ = patchReleases, 
            migrationMale_ = migrationMale,
            migrationFemale_ = migrationFemale, 
            migrationBatch_ = migrationBatch,
            reproductionType_ = reproductionType, 
            outputDirectory_ = outputDirectory,
            verbose_ = verbose)
}

