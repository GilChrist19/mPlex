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
#' @usage MakeReference_DaisyDrive(H, R, S, d, eta=NULL, phi=NULL,
#' omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
#' @param d Vector of background mutation rates at each locus
#' @param eta Named vector of mating fitness 
#' @param phi Named vector of sex ratios at emerence
#' @param omega Named vector of adult mortality increase
#' @param xiF Named vector of female pupatory success
#' @param xiM Named vector of male pupatory success
#' @param s Named vector of genotype-dependent fertility reduction
#'
#' @details This function creates a reference list for \code{\link{DaisyOffspring}}.
#' The number of drive elements is specified by the length of H. The final homing
#' rate doesn't matter, as the last piece of the drive has nothing to drive, but
#' it must be there for compatibility. R,S, and d must be the same length as H,
#' but will generally be the same number replicated that many times. It is assumed
#' that S is R/3, but this can be varied. This function is similar to \code{\link{MakeReference_Multiplex_mLoci}}
#' and \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return list of homing, cutting, and mendelian genotypes and rates
#'
#' @examples
#' H <- c(0.9,0.4,0) # This drive has 3 pieces
#' R <- c(0.001, 0.002, 0.003)
#' S <- c(0.0003, 0.0006, 0.001)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_DaisyDrive(H,R,S,d)
#'
#' @export
MakeReference_DaisyDrive <- function(H=c(0.9, 0.4, 0.7),R=c(0.0, 0.0, 0.0), S=R/3, d=c(0.0001, 0.0001, 0.0001),
                                     eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  #H is homing rate. The length of this vector determines the number of pieces
  # in the daisy drive. Each drive can have the same or different rates.
  
  #R is the NHEJ rate for deleterious alleles. Must be same length as H,
  # can be the same or different values.
  
  #S is the NHEJ rate for neutral alleles. Must be the same length as H,
  # can be the same or different values.
  
  #d is the background mutation rate. Must be the same length as H, can
  # have the same or different values.
  
  
  
  #Safety checks
  if(any( c(length(H),length(R), length(S), length(d)) != length(H))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(H) == length(R) == length(S) == length(d)"))
  }
  if(any(H[-length(H)]>1, R>1, S>1, d>1)){
    #last H doesn't matter because it may/may not exist and isnt' used
    return(cat("All rates must be less than or equal to 1\n"))
  }
  if(any((d+c(0, H[-length(H)])) > 1)){
    #need the driving piece's haming rate plus the background mutation rate
    #  of the piece being driven into.
    return(cat("Homing rates plus background mutation rates must be <= 1\n",
               "i.e. H+d <= 1\n"))
  }
  if(any((R+S) > 1)){
    return(cat("Negative and neutral repair rates must sum to <= 1\n",
               "i.s. R+S <= 1"))
  }
  

  if(!is.null(eta) && any(eta<0)){
    return(cat("All elements of eta must be greater than 0"))
  }
  if(!is.null(phi) && any(phi<0)){
    return(cat("All elements of phi must be greater than 0"))
  }
  if(!is.null(omega) && any(omega<0)){
    return(cat("All elements of omega must be greater than 0"))
  }
  if(!is.null(xiF) && any(xiF<0)){
    return(cat("All elements of xiF must be greater than 0"))
  }
  if(!is.null(xiM) && any(xiM<0)){
    return(cat("All elements of xiM must be greater than 0"))
  }
  if(!is.null(s) && any(s<0)){
    return(cat("All elements of s must be greater than 0"))
  }
  

  
  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Deleterious resistant
  #S = Neutral resistant
  gtype <- c("W", "H", "R", "S")
  Hshift <- c(0, H[-length(H)]) #because each pieces relies on the efficiency of the previous piece
  
  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))
  cuttingProbs <- matrix(data = 0, nrow = 3, ncol = length(H), dimnames = list(gtype[-2], NULL))
  
  homingProbs[1, ] <- 1-d-Hshift #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- Hshift*(1-R-S) #chance to become homing is H*1-H*R-H*S
  homingProbs[3, ] <- Hshift*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+Hshift*S #good resistant allele, from NHEJ and background mutation rate
  
  cuttingProbs[1, ] <- 1-d-Hshift*(R+S) #chance to stay W, 1-homingrate*mutations-background
  cuttingProbs[2, ] <- Hshift*R #chance to become R, cutting*NHEJ deleterious rate
  cuttingProbs[3, ] <- Hshift*S+d #become S, cutting*NHEJ neutral + background
  
  
  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(H))
  cutProbsList <- vector(mode = "list", length = length(H))
  homProbsList <- vector(mode = "list", length = length(H))
  
  mendAlleleList <- vector(mode = "list", length = length(H))
  cutAlleleList <- vector(mode = "list", length = length(H))
  homAlleleList <- vector(mode = "list", length = length(H))
  
  #fill the lists
  for(i in 1:length(H)){
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "S"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    mendProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    mendProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = mendProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      mendProbsList[[i]][[j]] <- mendProbsList[[i]][[j]][ logicalHold[[j]] ]
      mendAlleleList[[i]][[j]] <- names(mendProbsList[[i]][[j]])
    }
    
    cutProbsList[[i]]$W <- cuttingProbs[ ,i]
    cutProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
    cutProbsList[[i]]$R <- setNames(object = 1, nm = "R")
    cutProbsList[[i]]$S <- setNames(object = 1, nm = "S")
    
    #remove 0 probs things, and set allele names
    logicalHold <- lapply(X = cutProbsList[[i]], FUN = '!=', 0)
    for(j in 1:4){
      cutProbsList[[i]][[j]] <- cutProbsList[[i]][[j]][ logicalHold[[j]] ]
      cutAlleleList[[i]][[j]] <- names(cutProbsList[[i]][[j]])
    }
    
    homProbsList[[i]]$W <- homingProbs[ ,i]
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
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
#' @usage MakeReference_Multiplex_mLoci(H, R, S, d, eta=NULL, phi=NULL,
#' omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
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
#' The number of targets is specified by the length of H. R,S, and d must
#' be the same length as H, but will generally be the same number replicated that
#' many times. It is assumed that S is R/3, but this can be varied.
#' This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_oLocus}}
#'
#' @return list of homing and mendelian genotypes and rates
#'
#' @examples
#' H <- c(0.9,0.4,0.7) # This drive targets 3 loci
#' R <- c(0.001, 0.002, 0.003)
#' S <- c(0.0003, 0.0006, 0.001)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_mLoci(H,R,S,d)
#'
#' @export
MakeReference_Multiplex_mLoci <- function(H=c(0.9, 0.4, 0.7),R=c(0.0, 0.0, 0.0), S=R/3, d=c(0.0001, 0.0001, 0.0001),
                                          eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  #H is homing rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #R is the NHEJ rate for deleterious alleles. Must be same length as H,
  # can be the same or different values.
  
  #S is the NHEJ rate for neutral alleles. Must be the same length as H,
  # can be the same or different values.
  
  #d is the background mutation rate. Must be the same length as H, can
  # have the same or different values.
  
  
  #Safety checks
  if(any( c(length(H),length(R), length(S), length(d)) != length(H))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(H) == length(R) == length(S) == length(d)"))
  }
  if(any(H>1, R>1, S>1, d>1)){
    return(cat("All rates must be less than or equal to 1\n"))
  }
  if(any((d+H) > 1)){
    return(cat("Homing rates plus background mutation rates must be <= 1\n",
               "i.e. H+d <= 1\n"))
  }
  if(any((R+S) > 1)){
    return(cat("Negative and neutral repair rates must sum to <= 1\n",
               "i.s. R+S <= 1"))
  }
  
  
  if(!is.null(eta) && any(eta<0)){
    return(cat("All elements of eta must be greater than 0"))
  }
  if(!is.null(phi) && any(phi<0)){
    return(cat("All elements of phi must be greater than 0"))
  }
  if(!is.null(omega) && any(omega<0)){
    return(cat("All elements of omega must be greater than 0"))
  }
  if(!is.null(xiF) && any(xiF<0)){
    return(cat("All elements of xiF must be greater than 0"))
  }
  if(!is.null(xiM) && any(xiM<0)){
    return(cat("All elements of xiM must be greater than 0"))
  }
  if(!is.null(s) && any(s<0)){
    return(cat("All elements of s must be greater than 0"))
  }
  
  
  
  
  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Deleterious resistant
  #S = Neutral resistant
  gtype <- c("W", "H", "R", "S")
  
  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))
  
  homingProbs[1, ] <- 1-d-H #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- H*(1-R-S) #chance to become homing is H*1-H*R-H*S
  homingProbs[3, ] <- H*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+H*S #good resistant allele, from NHEJ and background mutation rate
  
  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(H))
  homProbsList <- vector(mode = "list", length = length(H))
  
  mendAlleleList <- vector(mode = "list", length = length(H))
  homAlleleList <- vector(mode = "list", length = length(H))
  
  #fill the lists
  for(i in 1:length(H)){
    #Mendelian Probabilities
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "S"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
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
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
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
#' @usage MakeReference_Multiplex_oLocus(H, R, S, d, eta=NULL, phi=NULL,
#' omega=NULL, xiF=NULL, xiM=NULL, s=NULL)
#'
#' @param H Vector of homing rates for each drive piece
#' @param R Vector of deleterious allele generation rates
#' @param S Vector of neutral allele generation rates
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
#' many times. It is assumed that S is R/3, but this can be varied.
#' This function is similar to \code{\link{MakeReference_DaisyDrive}} and
#' \code{\link{MakeReference_Multiplex_mLoci}}
#'
#' @return list of homing and mendelian genotypes and rates
#'
#' @examples
#' H <- c(0.9,0.4,0.7) # This drive has 3 targets at 1 locus
#' R <- c(0.001, 0.002, 0.003)
#' S <- c(0.0003, 0.0006, 0.001)
#' d <- c(0.0001, 0.0001, 0.0001)
#'
#' MakeReference_Multiplex_oLocus(H,R,S,d)
#'
#' @export
MakeReference_Multiplex_oLocus <- function(H=c(0.9),R=c(0.0), S=R/3, d=c(0.0001),
                                           eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  #H is homing rate. The length of this vector determines the number of loci
  # in the multiplex drive. Each drive can have the same or different rates.
  
  #R is the NHEJ rate for deleterious alleles. Must be same length as H,
  # can be the same or different values.
  
  #S is the NHEJ rate for neutral alleles. Must be the same length as H,
  # can be the same or different values.
  
  #d is the background mutation rate. Must be the same length as H, can
  # have the same or different values.
  
  
  
  #Safety checks
  if(any( c(length(H),length(R), length(S), length(d)) != length(H))){
    return(cat("All inputs must be the same length!\n",
               "i.e. length(H) == length(R) == length(S) == length(d)"))
  }
  if(any(H>1, R>1, S>1, d>1)){
    return(cat("All rates must be less than or equal to 1\n"))
  }
  if(any((d+H) > 1)){
    return(cat("Homing rates plus background mutation rates must be <= 1\n",
               "i.e. H+d <= 1\n"))
  }
  if(any((R+S) > 1)){
    return(cat("Negative and neutral repair rates must sum to <= 1\n",
               "i.s. R+S <= 1"))
  }
  
  
  if(!is.null(eta) && any(eta<0)){
    return(cat("All elements of eta must be greater than 0"))
  }
  if(!is.null(phi) && any(phi<0)){
    return(cat("All elements of phi must be greater than 0"))
  }
  if(!is.null(omega) && any(omega<0)){
    return(cat("All elements of omega must be greater than 0"))
  }
  if(!is.null(xiF) && any(xiF<0)){
    return(cat("All elements of xiF must be greater than 0"))
  }
  if(!is.null(xiM) && any(xiM<0)){
    return(cat("All elements of xiM must be greater than 0"))
  }
  if(!is.null(s) && any(s<0)){
    return(cat("All elements of s must be greater than 0"))
  }
  
  
  
  
  #setup allele letters
  #W = Wild-type
  #H = Homing
  #R = Deleterious resistant
  #S = Neutral resistant
  gtype <- c("W", "H", "R", "S")
  
  #matrix to hold homing probs, then fill it
  homingProbs <- matrix(data = 0, nrow = 4, ncol = length(H), dimnames = list(gtype, NULL))
  
  homingProbs[1, ] <- 1-d-H #chance to stay W is 1-homing-background mutation
  homingProbs[2, ] <- H*(1-R-S) #chance to become homing is H*1-H*R-H*S
  homingProbs[3, ] <- H*R   #NHEJ caused resistance, detrimentalt allele
  homingProbs[4, ] <- d+H*S #good resistant allele, from NHEJ and background mutation rate
  
  #set up lists to hold probabilities
  mendProbsList <- vector(mode = "list", length = length(H))
  homProbsList <- vector(mode = "list", length = length(H))
  
  mendAlleleList <- vector(mode = "list", length = length(H))
  homAlleleList <- vector(mode = "list", length = length(H))
  
  #fill the lists
  for(i in 1:length(H)){
    #Mendelian Probabilities
    mendProbsList[[i]]$W <- setNames(object = c(1-d[i], d[i]), nm = c("W", "S"))
    mendProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
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
    homProbsList[[i]]$H <- setNames(object = c(1-d[i], d[i]), nm = c("H", "S"))
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
