###############################################################################
#                                    ____  _
#                          _ __ ___ |  _ \| | _____  __
#                         | '_ ` _ \| |_) | |/ _ \ \/ /
#                         | | | | | |  __/| |  __/>  <
#                         |_| |_| |_|_|   |_|\___/_/\_\
#
###############################################################################
#    __  __                       _ _           ____ _
#   |  \/  | ___  ___  __ _ _   _(_) |_ ___    / ___| | __ _ ___ ___
#   | |\/| |/ _ \/ __|/ _` | | | | | __/ _ \  | |   | |/ _` / __/ __|
#   | |  | | (_) \__ \ (_| | |_| | | || (_) | | |___| | (_| \__ \__ \
#   |_|  |_|\___/|___/\__, |\__,_|_|\__\___/   \____|_|\__,_|___/___/
#                        |_|
###############################################################################

#' Mosquito Class Definition
#'
#' Mosquitos are the individual units being simulated in mPlex
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#'
#' @section **Constructor**:
#'  * genotype: The genotype of the Mosquito to create
#'  * age: The age of the Mosquito at creation
#'
#' @section **Methods**:
#'  * set_age: Set the age of the Mosquito
#'  * set_mate: Set the mate of the Mosquito
#'  * set_genotype: Set the genotype of the Mosquito
#'  * get_age: Return the age of the Mosquito
#'  * get_mate: Return the mate of the Mosquito
#'  * get_genotype: Return the genotype of the Mosquito
#'  * age_one_day: Increase the age of the Mosquito by 1
#'  * print_male: Return age and genotype as comma separated values
#'  * print_female: Return age, genotype, and mate as comma separated values
#'
#' @section **Fields**:
#'  * age: Age of Mosquito
#'  * mate: Genotype of the mate of Mosquito
#'  * genotype: Genotype of Mosquito
#'
#' @export
Mosquito <- function(genotype = NULL, age = NULL){
  #last attribute
  mate = NULL
  genotype <<- genotype
  age <<- age

  #getters/setters
  set_mate <- function(matGen=NULL){mate <<- matGen}



  #functions
  age_one_day <- function(){age <<- age + 1L}
  print_female = function(){file.path(age, genotype, mate, fsep = ",")}
  print_male = function(){file.path(age, genotype, fsep = ",")}

  #this is an R environment
  environment()
}


# Mosquito <- R6::R6Class(classname = "mosquito",
#                     portable = TRUE,
#                     cloneable = FALSE,
#                     lock_class = FALSE,
#                     lock_objects = FALSE,
#                     class = FALSE,
#
#                     # public memebers
#                     public = list(
#
#                       # constructor
#                       initialize = function(genotype=NULL, age=NULL){
#                         private$age = age
#                         private$mate = NULL
#                         private$genotype = genotype
#
#                       }, # end constructor
#
#                       # setters
#                       set_age = function(age=NULL){private$age = age},
#                       set_mate = function(mate=NULL){private$mate = mate},
#                       set_genotype = function(genotype=NULL){private$genotype = genotype},
#
#                       age_one_day = function() {private$age = private$age + 1},
#                       print_female = function(){file.path(private$age, private$genotype, private$mate, fsep = ",")},
#                       print_male = function(){file.path(private$age, private$genotype, fsep = ",")},
#                       #file.path can't handle nulls
#
#                       #getters
#                       get_age = function(){private$age},
#                       get_mate = function(){private$mate},
#                       get_genotype = function(){private$genotype}
#
#                     ), # end public
#
#                     private = list(
#
#                       # fields
#                       age = NULL,
#                       mate = NULL,
#                       genotype = NULL
#
#                     ) # end private
#
# ) # end class definition


#new
# exact replica of Mosquito class
# drop in, ready to go
# NewMosquito <- function(genotype=NULL, age=NULL){
#   #last attribute
#   mate = NULL
#
#   #getters/setters
#   set_age <- function(newAge=NULL){age <<- newAge}
#   set_mate <- function(matGen=NULL){mate <<- matGen}
#   set_genotype <- function(newGen=NULL){genotype <<- newGen}
#
#   get_age <- function(){age}
#   get_mate <- function(){mate}
#   get_genotype <- function(){genotype}
#
#   #functions
#   age_one_day <- function(){age <<- age + 1}
#   print_female = function(){file.path(age, genotype, mate, fsep = ",")}
#   print_male = function(){file.path(age, genotype, fsep = ",")}
#
#   #this is an R environment
#   environment()
# }
#
# # Smallest, fastest
# NewMosquitoSMALLEST <- function(genotype=NULL, age=NULL){
#   #last attribute
#   mate = NULL
#
#   # no getters/setters required
#
#   #functions
#   age_one_day = function(){age <<- age + 1}
#   print_female = function(){file.path(age, genotype, mate, fsep = ",")}
#   print_male = function(){file.path(age, genotype, fsep = ",")}
#
#   #this is an R environment
#   environment()
# }
