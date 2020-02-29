###############################################################################
#                            ____  __          ______          
#                 ____ ___  / __ \/ /__  _  __/ ____/___  ____ 
#                / __ `__ \/ /_/ / / _ \| |/_/ /   / __ \/ __ \
#               / / / / / / ____/ /  __/>  </ /___/ /_/ / /_/ /
#              /_/ /_/ /_/_/   /_/\___/_/|_|\____/ .___/ .___/ 
#                                               /_/   /_/      
###############################################################################
###############################################################################
# DAISY
###############################################################################

#' Daisy Drive Offspring Reference
#'
#' Create a list specifying the offspring probability distribution.
#'
#' @usage MakeReference_DaisyDrive(cRate, hRate, rRate, d, eta=NULL, phi=NULL,
#' omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param cRate Vector of cutting rates for each drive piece
#' @param hRate Vector of homing allele generation rates
#' @param rRate Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{DaisyOffspring}}.
#' The number of drive elements is specified by the length of cRate. The first element 
#' of each {c,h,r}Rate must be 0, because the first element behaves in a mendelian 
#' fashion, and the rates of the final piece don't matter, because there is not 
#' another piece to drive. cRate, hRate, and rRate must be the same length. This 
#' function is similar to \code{\link{MakeReference_Multiplex_mLoci}}
#' and \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return List of homing, cutting, and mendelian genotypes and rates.
#'
#' @examples
#' cRate <- c(0,0.4,0.7) # This drive targets 3 loci
#' hRate <- c(0, 0.2, 0.3)
#' rRate <- c(0, 0.006, 0.01)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_DaisyDrive(H,R,S,d)
#'
#' @export
MakeReference_DaisyDrive <- function(cRate=c(0, 1.0, 1.0), hRate=c(0,1.0,1.0),
                                          rRate=c(0,0,0), d=c(0.0001, 0.0001, 0.0001),
                                          eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  #cRate is cutting rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #hRate is the proper hoing rate for H alleles. Must be same length as cRate,
  # can be the same or different values.
  
  #rRate is the NHEJ rate for neutral alleles. Must be the same length as cRate,
  # can be the same or different values.
  
  #d is the background mutation rate. Must be the same length as cRate, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(cRate),length(hRate), length(rRate), length(d)) != length(cRate))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(cRate) == length(hRate) == length(rRate) == length(d)"))
  }
  if(any(cRate>1, hRate>1, rRate>1, d>1) || any(cRate<0, hRate<0, rRate<0, d<0)){
    return(cat("All rates must satisfy 0 <= rate <= 1\n"))
  }
  if(any(c(cRate[1],hRate[1],rRate[1]) != 0)){
    return(cat("First element is Mendelian, so all rates must be 0.\n"))
  }
  

  if(!is.null(eta) && any(eta<0)){
    return(cat("All elements of eta must be >= 0"))
  }
  if(!is.null(phi) && any(phi<0)){
    return(cat("All elements of phi must be >= 0"))
  }
  if(!is.null(omega) && any(omega<0)){
    return(cat("All elements of omega must be >= 0"))
  }
  if(!is.null(xiF) && any(xiF<0)){
    return(cat("All elements of xiF must be >= 0"))
  }
  if(!is.null(xiM) && any(xiM<0)){
    return(cat("All elements of xiM must be >= 0"))
  }
  if(!is.null(s) && any(s<0)){
    return(cat("All elements of s must be >= 0"))
  }
  

  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Neutral resistant
  #S = Deleterious resistant
  gtype <- c("W", "H", "R", "S")

  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(cRate), dimnames = list(gtype, NULL))
  cuttingProbs <- matrix(data = 0, nrow = 3, ncol = length(cRate), dimnames = list(gtype[-2], NULL))
  
  homingProbs[1, ] <- (1-d)*(1-cRate) #chance to stay W is (1-background mutation) * (1-cutting rate)
  homingProbs[2, ] <- (1-d)*cRate*hRate #chance to become H is (1-background)*cutting*homing
  homingProbs[3, ] <- d + (1-d)*cRate*(1-hRate)*rRate #NHEJ caused good resistance, 
  homingProbs[4, ] <- (1-d)*cRate*(1-hRate)*(1-rRate) #bad resistant allele, from NHEJ and background mutation rate
  
  cuttingProbs[1, ] <- (1-d)*(1-cRate) #chance to stay W is (1-background mutation) * (1-cutting rate)
  cuttingProbs[2, ] <- d + (1-d)*cRate*rRate #NHEJ caused good resistance, 
  cuttingProbs[3, ] <- (1-d)*cRate*(1-rRate) #bad resistant allele, from NHEJ and background mutation rate
  
  
  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(cRate))
  cutProbsList <- vector(mode = "list", length = length(cRate))
  homProbsList <- vector(mode = "list", length = length(cRate))
  
  mendAlleleList <- vector(mode = "list", length = length(cRate))
  cutAlleleList <- vector(mode = "list", length = length(cRate))
  homAlleleList <- vector(mode = "list", length = length(cRate))
  
  #fill the lists
  for(i in 1:length(cRate)){
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "R"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
      mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
    }
    
    cutProbsList[[i]]$W <- cuttingProbs[ ,i]
    cutProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    cutProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    cutProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = cutProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      cutProbsList[[i]][[j]] <- cutProbsList[[i]][[j]][ logicalHold[[j]] ]
      cutAlleleList[[i]][[j]] <- names(cutProbsList[[i]][[j]])
    }
    
    homProbsList[[i]]$W <- homingProbs[ ,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    homProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
      homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
    }
  }
  
  
  #set genotype-specific parameters
  eta_ = phi_ = omega_ = xiF_ = xiM_ = s_ = numeric(length = 0)
  
  if(!is.null(eta)){eta_ <- unlist(eta)} # unlist protects against list input. doesn't harm vectors
  if(!is.null(phi)) {phi_ <- unlist(phi)}
  if(!is.null(omega)){omega_ <- unlist(omega)}
  if(!is.null(xiF)){xiF_ <- unlist(xiF)}
  if(!is.null(xiM)){xiM_ <- unlist(xiM)}
  if(!is.null(s)){s_ <- unlist(s)}
    
    
  return(list(
    mendelian = mendProbsList,
    cutting = cutProbsList,
    homing = homProbsList,
    mendelianAlleles = mendAlleleList,
    homingAlleles = homAlleleList,
    cuttingAlleles = cutAlleleList,
    eta = eta_,
    phi = phi_,
    omega = omega_,
    xiF = xiF_,
    xiM = xiM_,
    s = s_))
  
}

###############################################################################
# MLOCI
###############################################################################

#' mPlex Multiple Loci Offspring Reference
#'
#' Create a list specifying the offspring probability distribution.
#'
#' @usage MakeReference_Multiplex_mLoci(cRate, hRate, rRate, d, eta=NULL, phi=NULL,
#' omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param cRate Vector of cutting rates for each drive piece
#' @param hRate Vector of homing allele generation rates
#' @param rRate Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{MultiplexOffspring_mLoci}}.
#' Each drive element targets one locus. Each locus is independent of the rest.
#' The number of targets is specified by the length of H. R, S, and d must
#' be the same length as H, but will generally be the same number replicated that
#' many times. This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return List of homing and mendelian genotypes and rates.
#'
#' @examples
#' cRate <- c(0.9,0.4,0.7) # This drive targets 3 loci
#' hRate <- c(0.1, 0.2, 0.3)
#' rRate <- c(0.003, 0.006, 0.01)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_mLoci(cRate,hRate,rRate,d)
#'
#' @export
MakeReference_Multiplex_mLoci <- function(cRate=c(1.0, 1.0, 1.0), hRate=c(1.0,1.0,1.0),
                                          rRate=c(0,0,0), d=c(0.0001, 0.0001, 0.0001),
                                          eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  #cRate is cutting rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #hRate is the proper hoing rate for H alleles. Must be same length as cRate,
  # can be the same or different values.
  
  #rRate is the NHEJ rate for neutral alleles. Must be the same length as cRate,
  # can be the same or different values.
  
  #d is the background mutation rate. Must be the same length as cRate, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(cRate),length(hRate), length(rRate), length(d)) != length(cRate))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(cRate) == length(hRate) == length(rRate) == length(d)"))
  }
  if(any(cRate>1, hRate>1, rRate>1, d>1) || any(cRate<0, hRate<0, rRate<0, d<0)){
    return(cat("All rates must satisfy 0 <= rate <= 1\n"))
  }
  
  
  if(!is.null(eta) && any(eta<0)){
    return(cat("All elements of eta must be >= 0"))
  }
  if(!is.null(phi) && any(phi<0)){
    return(cat("All elements of phi must be >= 0"))
  }
  if(!is.null(omega) && any(omega<0)){
    return(cat("All elements of omega must be >= 0"))
  }
  if(!is.null(xiF) && any(xiF<0)){
    return(cat("All elements of xiF must be >= 0"))
  }
  if(!is.null(xiM) && any(xiM<0)){
    return(cat("All elements of xiM must be >= 0"))
  }
  if(!is.null(s) && any(s<0)){
    return(cat("All elements of s must be >= 0"))
  }
  
  
  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Neutral resistant
  #S = Deleterious resistant
  gtype <- c("W", "H", "R", "S")
  
  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(cRate), dimnames = list(gtype, NULL))
  
  homingProbs[1, ] <- (1-d)*(1-cRate) #chance to stay W is (1-background mutation) * (1-cutting rate)
  homingProbs[2, ] <- (1-d)*cRate*hRate #chance to become H is (1-background)*cutting*homing
  homingProbs[3, ] <- d + (1-d)*cRate*(1-hRate)*rRate #NHEJ caused good resistance, 
  homingProbs[4, ] <- (1-d)*cRate*(1-hRate)*(1-rRate) #bad resistant allele, from NHEJ and background mutation rate
  
  
  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(cRate))
  homProbsList <- vector(mode = "list", length = length(cRate))
  
  mendAlleleList <- vector(mode = "list", length = length(cRate))
  homAlleleList <- vector(mode = "list", length = length(cRate))
  
  #fill the lists
  for(i in 1:length(cRate)){
    #Mendelian Probabilities
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "R"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
      mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
    }
    
    #Homing Probabilities
    homProbsList[[i]]$W <- homingProbs[ ,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    homProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
      homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
    }
  }
  
  
  
  #set genotype-specific parameters
  eta_ = phi_ = omega_ = xiF_ = xiM_ = s_ = numeric(length = 0)
  
  if(!is.null(eta)){eta_ <- unlist(eta)} # unlist protects against list input. doesn't harm vectors
  if(!is.null(phi)) {phi_ <- unlist(phi)}
  if(!is.null(omega)){omega_ <- unlist(omega)}
  if(!is.null(xiF)){xiF_ <- unlist(xiF)}
  if(!is.null(xiM)){xiM_ <- unlist(xiM)}
  if(!is.null(s)){s_ <- unlist(s)}
  
  
  return(list(
    mendelian = mendProbsList,
    homing = homProbsList,
    mendelianAlleles = mendAlleleList,
    homingAlleles = homAlleleList,
    eta = eta_,
    phi = phi_,
    omega = omega_,
    xiF = xiF_,
    xiM = xiM_,
    s = s_))
  
}

###############################################################################
# OLOCUS
###############################################################################

#' mPlex One Locus Offspring Reference
#'
#' Create a list specifying the offspring probability distribution.
#'
#' @usage MakeReference_Multiplex_oLocus(cRate, hRate, rRate, d, eta=NULL, phi=NULL,
#' omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param cRate Vector of cutting rates for each drive piece
#' @param hRate Vector of homing allele generation rates
#' @param rRate Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{MultiplexOffspring_oLocus}}.
#' It assumes multiple targeting gRNAs for 1 locus, all targets segregate together.
#' The length of each allele is specified by the length of H. R, S, and d must
#' be the same length as H, but will generally be the same number replicated that
#' many times. 
#' This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_mLoci}}
#'
#' @return List of homing and mendelian genotypes and rates.
#'
#' @examples
#' cRate <- c(0.9,0.4,0.7) # This drive targets 3 loci
#' hRate <- c(0.1, 0.2, 0.3)
#' rRate <- c(0.003, 0.006, 0.01)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_oLocus(H,R,S,d)
#'
#' @export
MakeReference_Multiplex_oLocus <- function(cRate=c(1.0, 1.0, 1.0), hRate=c(1.0,1.0,1.0),
                                           rRate=c(0,0,0), d=c(0.0001, 0.0001, 0.0001),
                                           eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  #cRate is cutting rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #hRate is the proper hoing rate for H alleles. Must be same length as cRate,
  # can be the same or different values.
  
  #rRate is the NHEJ rate for neutral alleles. Must be the same length as cRate,
  # can be the same or different values.
  
  #d is the background mutation rate. Must be the same length as cRate, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(cRate),length(hRate), length(rRate), length(d)) != length(cRate))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(cRate) == length(hRate) == length(rRate) == length(d)"))
  }
  if(any(cRate>1, hRate>1, rRate>1, d>1) || any(cRate<0, hRate<0, rRate<0, d<0)){
    return(cat("All rates must satisfy 0 <= rate <= 1\n"))
  }
  
  
  if(!is.null(eta) && any(eta<0)){
    return(cat("All elements of eta must be >= 0"))
  }
  if(!is.null(phi) && any(phi<0)){
    return(cat("All elements of phi must be >= 0"))
  }
  if(!is.null(omega) && any(omega<0)){
    return(cat("All elements of omega must be >= 0"))
  }
  if(!is.null(xiF) && any(xiF<0)){
    return(cat("All elements of xiF must be >= 0"))
  }
  if(!is.null(xiM) && any(xiM<0)){
    return(cat("All elements of xiM must be >= 0"))
  }
  if(!is.null(s) && any(s<0)){
    return(cat("All elements of s must be >= 0"))
  }
  
  
  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Neutral resistant
  #S = Deleterious resistant
  gtype <- c("W", "H", "R", "S")
  
  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(cRate), dimnames = list(gtype, NULL))
  
  homingProbs[1, ] <- (1-d)*(1-cRate) #chance to stay W is (1-background mutation) * (1-cutting rate)
  homingProbs[2, ] <- (1-d)*cRate*hRate #chance to become H is (1-background)*cutting*homing
  homingProbs[3, ] <- d + (1-d)*cRate*(1-hRate)*rRate #NHEJ caused good resistance, 
  homingProbs[4, ] <- (1-d)*cRate*(1-hRate)*(1-rRate) #bad resistant allele, from NHEJ and background mutation rate
  
  
  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(cRate))
  homProbsList <- vector(mode = "list", length = length(cRate))
  
  mendAlleleList <- vector(mode = "list", length = length(cRate))
  homAlleleList <- vector(mode = "list", length = length(cRate))
  
  #fill the lists
  for(i in 1:length(cRate)){
    #Mendelian Probabilities
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "R"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
      mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
    }
    
    #Homing Probabilities
    homProbsList[[i]]$W <- homingProbs[ ,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "R"))
    homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    homProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
      homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
    }
  }
  
  
  #set genotype-specific parameters
  eta_ = phi_ = omega_ = xiF_ = xiM_ = s_ = numeric(length = 0)
  
  if(!is.null(eta)){eta_ <- unlist(eta)} # unlist protects against list input. doesn't harm vectors
  if(!is.null(phi)) {phi_ <- unlist(phi)}
  if(!is.null(omega)){omega_ <- unlist(omega)}
  if(!is.null(xiF)){xiF_ <- unlist(xiF)}
  if(!is.null(xiM)){xiM_ <- unlist(xiM)}
  if(!is.null(s)){s_ <- unlist(s)}
  
  
  return(list(
    mendelian = mendProbsList,
    homing = homProbsList,
    mendelianAlleles = mendAlleleList,
    homingAlleles = homAlleleList,
    eta = eta_,
    phi = phi_,
    omega = omega_,
    xiF = xiF_,
    xiM = xiM_,
    s = s_))
  
}

