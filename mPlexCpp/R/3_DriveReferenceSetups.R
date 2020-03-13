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
#' @usage MakeReference_DaisyDrive(cRateM, hRateM, rRateM, dM,
#' cRateF=cRateM, hRateF=hRateM, rRateF=rRateM, dF=dM,
#' eta=NULL, phi=NULL, omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param cRateM Vector of cutting rates for each drive piece in males
#' @param hRateM Vector of homing allele generation rates in males
#' @param rRateM Vector of neutral allele generation rates in males
#' @param cRateF Vector of cutting rates for each drive piece in females
#' @param hRateF Vector of homing allele generation rates in females
#' @param rRateF Vector of neutral allele generation rates in females
#' @param dM Vector of background mutation rates at each locus in males
#' @param dF Vector of background mutation rates at each locus in females
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{DaisyOffspring}}.
#' The number of drive elements is specified by the length of cRateM. The first element 
#' of each {c,h,r}Rate must be 0, because the first element behaves in a mendelian 
#' fashion, and the rates of the final piece don't matter, because there is not 
#' another piece to drive. cRateM, hRateM, and rRateM must be the same length. Each 
#' rate can be different between males and females. The default is for females to be 
#' the same as males, and the user only needs to specify male rates. This 
#' function is similar to \code{\link{MakeReference_Multiplex_mLoci}}
#' and \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return List of homing, cutting, and mendelian genotypes and rates.
#'
#' @examples
#' cRateM <- c(0,0.4,0.7) # This drive targets 3 loci
#' hRateM <- c(0, 0.2, 0.3)
#' rRateM <- c(0, 0.006, 0.01)
#' dM <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_DaisyDrive(cRateM=cRateM, hRateM=hRateM, rRateM=rRateM,dM=dM)
#'
#' @export
MakeReference_DaisyDrive <- function(cRateM=c(0, 1.0, 1.0), hRateM=c(0,1.0,1.0), rRateM=c(0,0,0),
                                     cRateF=cRateM, hRateF=hRateM, rRateF=rRateM,
                                     dM=c(0.0001, 0.0001, 0.0001), dF=dM,
                                     eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  #cRateM is cutting rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #hRateM is the proper hoing rate for H alleles. Must be same length as cRateM,
  # can be the same or different values.
  
  #rRateM is the NHEJ rate for neutral alleles. Must be the same length as cRateM,
  # can be the same or different values.
  
  #dM is the background mutation rate. Must be the same length as cRateM, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(cRateM), length(hRateM), length(rRateM), length(dM), length(cRateF),length(hRateF), length(rRateF), length(dF)) != length(cRateM))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(cRateM) == length(hRateM) == length(rRateM) == length(dM)"))
  }
  if(any(cRateM>1, hRateM>1, rRateM>1, dM>1, cRateF>1, hRateF>1, rRateF>1, dF>1) || any(cRateM<0, hRateM<0, rRateM<0, dM<0, cRateF<0, hRateF<0, rRateF<0, dF<0)){
    return(cat("All rates must satisfy 0 <= rate <= 1\n"))
  }
  if(any(c(cRateM[1],hRateM[1],rRateM[1],cRateF[1],hRateF[1],rRateF[1]) != 0)){
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

  #matrix to hold homing probs
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(cRateM), dimnames = list(gtype, NULL))
  cuttingProbs <- matrix(data = 0, nrow = 3, ncol = length(cRateM), dimnames = list(gtype[-2], NULL))
  
  # setup return lists
  mendProbsRet <- list("female"=NULL,"male"=NULL)
  cutProbsRet <- list("female"=NULL,"male"=NULL)
  homProbsRet <- list("female"=NULL,"male"=NULL)
  mendAlleleRet <- list("female"=NULL,"male"=NULL)
  homAlleleRet <- list("female"=NULL,"male"=NULL)
  cutAlleleRet <- list("female"=NULL,"male"=NULL)
  
  
  # loop over both sexes
  sexC <- list(cRateF,cRateM)
  sexH <- list(hRateF,hRateM)
  sexR <- list(rRateF,rRateM)
  sexD <- list(dF,dM)
  for(sex in 1:2){
    # fill probs matrices
    homingProbs[1, ] <- (1-sexD[[sex]])*(1-sexC[[sex]]) #chance to stay W is (1-background mutation) * (1-cutting rate)
    homingProbs[2, ] <- (1-sexD[[sex]])*sexC[[sex]]*sexH[[sex]] #chance to become H is (1-background)*cutting*homing
    homingProbs[3, ] <- sexD[[sex]] + (1-sexD[[sex]])*sexC[[sex]]*(1-sexH[[sex]])*sexR[[sex]] #NHEJ caused good resistance, 
    homingProbs[4, ] <- (1-sexD[[sex]])*sexC[[sex]]*(1-sexH[[sex]])*(1-sexR[[sex]]) #bad resistant allele, from NHEJ and background mutation rate
    
    cuttingProbs[1, ] <- (1-sexD[[sex]])*(1-sexC[[sex]]) #chance to stay W is (1-background mutation) * (1-cutting rate)
    cuttingProbs[2, ] <- sexD[[sex]] + (1-sexD[[sex]])*sexC[[sex]]*sexR[[sex]] #NHEJ caused good resistance, 
    cuttingProbs[3, ] <- (1-sexD[[sex]])*sexC[[sex]]*(1-sexR[[sex]]) #bad resistant allele, from NHEJ and background mutation rate
    
    
    #set up lists to hold probabilities
    mendProbsList <- vector(mode = "list", length = length(cRateM))
    cutProbsList <- vector(mode = "list", length = length(cRateM))
    homProbsList <- vector(mode = "list", length = length(cRateM))
    
    mendAlleleList <- vector(mode = "list", length = length(cRateM))
    cutAlleleList <- vector(mode = "list", length = length(cRateM))
    homAlleleList <- vector(mode = "list", length = length(cRateM))
    
    #fill the lists
    for(i in 1:length(cRateM)){
      mendProbsList[[i]]$W <- setNames(object = c(1-dM[i], dM[i]), nm = c("W", "R"))
      mendProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
        mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
      }
      
      cutProbsList[[i]]$W <- cuttingProbs[ ,i]
      cutProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      cutProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      cutProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = cutProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        cutProbsList[[i]][[j]] <- cutProbsList[[i]][[j]][ logicalHold[[j]] ]
        cutAlleleList[[i]][[j]] <- names(cutProbsList[[i]][[j]])
      }
      
      homProbsList[[i]]$W <- homingProbs[ ,i]
      homProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      homProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
        homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
      }
    } # end loop over loci
    
    
    # fill return list
    mendProbsRet[[sex]] <- mendProbsList
    cutProbsRet[[sex]] <- cutProbsList
    homProbsRet[[sex]] <- homProbsList
    mendAlleleRet[[sex]] <- mendAlleleList
    homAlleleRet[[sex]] <- homAlleleList
    cutAlleleRet[[sex]] <- cutAlleleList
    
  } # end loop over male/female
  

  #set genotype-specific parameters
  eta_ = phi_ = omega_ = xiF_ = xiM_ = s_ = numeric(length = 0)
  
  if(!is.null(eta)){eta_ <- unlist(eta)} # unlist protects against list input. doesn't harm vectors
  if(!is.null(phi)) {phi_ <- unlist(phi)}
  if(!is.null(omega)){omega_ <- unlist(omega)}
  if(!is.null(xiF)){xiF_ <- unlist(xiF)}
  if(!is.null(xiM)){xiM_ <- unlist(xiM)}
  if(!is.null(s)){s_ <- unlist(s)}
    
    
  return(list(
    mendelian = mendProbsRet,
    cutting = cutProbsRet,
    homing = homProbsRet,
    mendelianAlleles = mendAlleleRet,
    homingAlleles = homAlleleRet,
    cuttingAlleles = cutAlleleRet,
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
#' @usage MakeReference_Multiplex_mLoci(cRateM, hRateM, rRateM, dM,
#' cRateF=cRateM, hRateF=hRateM, rRateF=rRateM, dF=dM,
#' eta=NULL, phi=NULL, omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param cRateM Vector of cutting rates for each drive piece in males
#' @param hRateM Vector of homing allele generation rates in males
#' @param rRateM Vector of neutral allele generation rates in males
#' @param cRateF Vector of cutting rates for each drive piece in females
#' @param hRateF Vector of homing allele generation rates in females
#' @param rRateF Vector of neutral allele generation rates in females
#' @param dM Vector of background mutation rates at each locus in males
#' @param dF Vector of background mutation rates at each locus in females
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{MultiplexOffspring_mLoci}}.
#' Each drive element targets one locus. Each locus is independent of the rest.
#' The number of targets is specified by the length of cRateM. hRateM, rRateM, and dM must
#' be the same length as cRateM. Each rate can be different between males and females. 
#' The default is for females to be the same as males, and the user only needs to 
#' specify male rates. This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return List of homing and mendelian genotypes and rates.
#'
#' @examples
#' cRateM <- c(0.9,0.4,0.7) # This drive targets 3 loci
#' hRateM <- c(0.1, 0.2, 0.3)
#' rRateM <- c(0.003, 0.006, 0.01)
#' dM <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_mLoci(cRateM=cRateM, hRateM=hRateM, rRateM=rRateM,dM=dM)
#'
#' @export
MakeReference_Multiplex_mLoci <- function(cRateM=c(1.0, 1.0, 1.0), hRateM=c(1.0,1.0,1.0), rRateM=c(0,0,0),
                                          cRateF=cRateM, hRateF=hRateM, rRateF=rRateM,
                                          dM=c(0.0001, 0.0001, 0.0001), dF=dM,
                                          eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  #cRateM is cutting rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #hRateM is the proper hoing rate for H alleles. Must be same length as cRateM,
  # can be the same or different values.
  
  #rRateM is the NHEJ rate for neutral alleles. Must be the same length as cRateM,
  # can be the same or different values.
  
  #dM is the background mutation rate. Must be the same length as cRateM, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(cRateM), length(hRateM), length(rRateM), length(dM), length(cRateF),length(hRateF), length(rRateF), length(dF)) != length(cRateM))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(cRateM) == length(hRateM) == length(rRateM) == length(dM)"))
  }
  if(any(cRateM>1, hRateM>1, rRateM>1, dM>1, cRateF>1, hRateF>1, rRateF>1, dF>1) || any(cRateM<0, hRateM<0, rRateM<0, dM<0, cRateF<0, hRateF<0, rRateF<0, dF<0)){
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
  
  #matrix to hold homing probs
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(cRateM), dimnames = list(gtype, NULL))
  
  # setup return lists
  mendProbsRet <- list("female"=NULL,"male"=NULL)
  homProbsRet <- list("female"=NULL,"male"=NULL)
  mendAlleleRet <- list("female"=NULL,"male"=NULL)
  homAlleleRet <- list("female"=NULL,"male"=NULL)

  
  # loop over both sexes
  sexC <- list(cRateF,cRateM)
  sexH <- list(hRateF,hRateM)
  sexR <- list(rRateF,rRateM)
  sexD <- list(dF,dM)
  for(sex in 1:2){
    # fill probs matrices
    homingProbs[1, ] <- (1-sexD[[sex]])*(1-sexC[[sex]]) #chance to stay W is (1-background mutation) * (1-cutting rate)
    homingProbs[2, ] <- (1-sexD[[sex]])*sexC[[sex]]*sexH[[sex]] #chance to become H is (1-background)*cutting*homing
    homingProbs[3, ] <- sexD[[sex]] + (1-sexD[[sex]])*sexC[[sex]]*(1-sexH[[sex]])*sexR[[sex]] #NHEJ caused good resistance, 
    homingProbs[4, ] <- (1-sexD[[sex]])*sexC[[sex]]*(1-sexH[[sex]])*(1-sexR[[sex]]) #bad resistant allele, from NHEJ and background mutation rate
    
    
    #set up lists to hold probabilities
    mendProbsList <- vector(mode = "list", length = length(cRateM))
    homProbsList <- vector(mode = "list", length = length(cRateM))
    
    mendAlleleList <- vector(mode = "list", length = length(cRateM))
    homAlleleList <- vector(mode = "list", length = length(cRateM))
    
    #fill the lists
    for(i in 1:length(cRateM)){
      mendProbsList[[i]]$W <- setNames(object = c(1-dM[i], dM[i]), nm = c("W", "R"))
      mendProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
        mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
      }
      
      #Homing probabilities
      homProbsList[[i]]$W <- homingProbs[ ,i]
      homProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      homProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
        homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
      }
    } # end loop over loci
    
    
    # fill return list
    mendProbsRet[[sex]] <- mendProbsList
    homProbsRet[[sex]] <- homProbsList
    mendAlleleRet[[sex]] <- mendAlleleList
    homAlleleRet[[sex]] <- homAlleleList

  } # end loop over male/female
  
  
  #set genotype-specific parameters
  eta_ = phi_ = omega_ = xiF_ = xiM_ = s_ = numeric(length = 0)
  
  if(!is.null(eta)){eta_ <- unlist(eta)} # unlist protects against list input. doesn't harm vectors
  if(!is.null(phi)) {phi_ <- unlist(phi)}
  if(!is.null(omega)){omega_ <- unlist(omega)}
  if(!is.null(xiF)){xiF_ <- unlist(xiF)}
  if(!is.null(xiM)){xiM_ <- unlist(xiM)}
  if(!is.null(s)){s_ <- unlist(s)}
  
  
  return(list(
    mendelian = mendProbsRet,
    homing = homProbsRet,
    mendelianAlleles = mendAlleleRet,
    homingAlleles = homAlleleRet,
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
#' @usage MakeReference_Multiplex_oLocus(cRateM, hRateM, rRateM, dM,
#' cRateF=cRateM, hRateF=hRateM, rRateF=rRateM, dF=dM,
#' eta=NULL, phi=NULL, omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param cRateM Vector of cutting rates for each drive piece in males
#' @param hRateM Vector of homing allele generation rates in males
#' @param rRateM Vector of neutral allele generation rates in males
#' @param cRateF Vector of cutting rates for each drive piece in females
#' @param hRateF Vector of homing allele generation rates in females
#' @param rRateF Vector of neutral allele generation rates in females
#' @param dM Vector of background mutation rates at each locus in males
#' @param dF Vector of background mutation rates at each locus in females
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{MultiplexOffspring_oLocus}}.
#' It assumes multiple targeting gRNAs for 1 locus, all targets segregate together.
#' The length of each allele is specified by the length of cRateM. hRateM, rRateM, and dM must
#' be the same length as cRateM. Each rate can be different between males and females. 
#' The default is for females to be the same as males, and the user only needs to 
#' specify male rates.
#' This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_mLoci}}
#'
#' @return List of homing and mendelian genotypes and rates.
#'
#' @examples
#' cRateM <- c(0.9,0.4,0.7) # This drive targets 3 loci
#' hRateM <- c(0.1, 0.2, 0.3)
#' rRateM <- c(0.003, 0.006, 0.01)
#' dM <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_oLocus(cRateM=cRateM, hRateM=hRateM, rRateM=rRateM,dM=dM)
#'
#' @export
MakeReference_Multiplex_oLocus <- function(cRateM=c(1.0, 1.0, 1.0), hRateM=c(1.0,1.0,1.0), rRateM=c(0,0,0),
                                           cRateF=cRateM, hRateF=hRateM, rRateF=rRateM,
                                           dM=c(0.0001, 0.0001, 0.0001), dF=dM,
                                           eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  #cRateM is cutting rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #hRateM is the proper hoing rate for H alleles. Must be same length as cRateM,
  # can be the same or different values.
  
  #rRateM is the NHEJ rate for neutral alleles. Must be the same length as cRateM,
  # can be the same or different values.
  
  #dM is the background mutation rate. Must be the same length as cRateM, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(cRateM), length(hRateM), length(rRateM), length(dM), length(cRateF),length(hRateF), length(rRateF), length(dF)) != length(cRateM))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(cRateM) == length(hRateM) == length(rRateM) == length(dM)"))
  }
  if(any(cRateM>1, hRateM>1, rRateM>1, dM>1, cRateF>1, hRateF>1, rRateF>1, dF>1) || any(cRateM<0, hRateM<0, rRateM<0, dM<0, cRateF<0, hRateF<0, rRateF<0, dF<0)){
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
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(cRateM), dimnames = list(gtype, NULL))
  
  # setup return lists
  mendProbsRet <- list("female"=NULL,"male"=NULL)
  homProbsRet <- list("female"=NULL,"male"=NULL)
  mendAlleleRet <- list("female"=NULL,"male"=NULL)
  homAlleleRet <- list("female"=NULL,"male"=NULL)
  
  
  # loop over both sexes
  sexC <- list(cRateF,cRateM)
  sexH <- list(hRateF,hRateM)
  sexR <- list(rRateF,rRateM)
  sexD <- list(dF,dM)
  for(sex in 1:2){
    # fill probs matrices
    homingProbs[1, ] <- (1-sexD[[sex]])*(1-sexC[[sex]]) #chance to stay W is (1-background mutation) * (1-cutting rate)
    homingProbs[2, ] <- (1-sexD[[sex]])*sexC[[sex]]*sexH[[sex]] #chance to become H is (1-background)*cutting*homing
    homingProbs[3, ] <- sexD[[sex]] + (1-sexD[[sex]])*sexC[[sex]]*(1-sexH[[sex]])*sexR[[sex]] #NHEJ caused good resistance, 
    homingProbs[4, ] <- (1-sexD[[sex]])*sexC[[sex]]*(1-sexH[[sex]])*(1-sexR[[sex]]) #bad resistant allele, from NHEJ and background mutation rate
    
    
    #set up lists to hold probabilities
    mendProbsList <- vector(mode = "list", length = length(cRateM))
    homProbsList <- vector(mode = "list", length = length(cRateM))
    
    mendAlleleList <- vector(mode = "list", length = length(cRateM))
    homAlleleList <- vector(mode = "list", length = length(cRateM))
    
    #fill the lists
    for(i in 1:length(cRateM)){
      mendProbsList[[i]]$W <- setNames(object = c(1-dM[i], dM[i]), nm = c("W", "R"))
      mendProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
        mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
      }
      
      #Homing probabilities
      homProbsList[[i]]$W <- homingProbs[ ,i]
      homProbsList[[i]]$H <- setNames(object = c(1-dM[i], dM[i]), nm = c("H", "R"))
      homProbsList[[i]]$R <- setNames(object = 1, nm = "R")
      homProbsList[[i]]$S <- setNames(object = 1, nm = "S")
      
      #remove 0 probs things, and set allele names
      logicalHold <- lapply(X = homProbsList[[i]], FUN = '!=', 0)
      for(j in 1:4){
        homProbsList[[i]][[j]] <- homProbsList[[i]][[j]][ logicalHold[[j]] ]
        homAlleleList[[i]][[j]] <- names(homProbsList[[i]][[j]])
      }
    } # end loop over loci
    
    
    # fill return list
    mendProbsRet[[sex]] <- mendProbsList
    homProbsRet[[sex]] <- homProbsList
    mendAlleleRet[[sex]] <- mendAlleleList
    homAlleleRet[[sex]] <- homAlleleList
    
  } # end loop over male/female
  
  
  #set genotype-specific parameters
  eta_ = phi_ = omega_ = xiF_ = xiM_ = s_ = numeric(length = 0)
  
  if(!is.null(eta)){eta_ <- unlist(eta)} # unlist protects against list input. doesn't harm vectors
  if(!is.null(phi)) {phi_ <- unlist(phi)}
  if(!is.null(omega)){omega_ <- unlist(omega)}
  if(!is.null(xiF)){xiF_ <- unlist(xiF)}
  if(!is.null(xiM)){xiM_ <- unlist(xiM)}
  if(!is.null(s)){s_ <- unlist(s)}
  
  
  return(list(
    mendelian = mendProbsRet,
    homing = homProbsRet,
    mendelianAlleles = mendAlleleRet,
    homingAlleles = homAlleleRet,
    eta = eta_,
    phi = phi_,
    omega = omega_,
    xiF = xiF_,
    xiM = xiM_,
    s = s_))
  
}

