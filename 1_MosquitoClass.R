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

Mosquito <- R6::R6Class(classname = "mosquito",
                    portable = TRUE,
                    cloneable = FALSE,
                    lock_class = FALSE,
                    lock_objects = FALSE,

                    # public memebers
                    public = list(

                      # constructor
                      initialize = function(genotype=NULL, age=NULL){
                        private$age = age
                        private$mate = NULL
                        private$genotype = genotype

                      }, # end constructor

                      # setters
                      set_age = function(age=NULL){private$age = age},
                      set_mate = function(mate=NULL){private$mate = mate},
                      set_genotype = function(genotype=NULL){private$genotype = genotype},

                      age_one_day = function() {private$age = private$age + 1},
                      print_female = function(){file.path(private$age, private$genotype, private$mate, fsep = ",")},
                      print_male = function(){file.path(private$age, private$genotype, fsep = ",")},
                      #file.path can't handle nulls

                      #getters
                      get_age = function(){private$age},
                      get_mate = function(){private$mate},
                      get_genotype = function(){private$genotype}

                    ), # end public

                    private = list(

                      # fields
                      age = NULL,
                      mate = NULL,
                      genotype = NULL

                    ) # end private

) # end class definition
