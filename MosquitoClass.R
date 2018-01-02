#mosquito class
#Jared Bennett


# This is the class definition for mosquitoes

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
                        private$mate = "NA"
                        private$genotype = genotype

                      }, # end constructor

                      # setters
                      set_age = function(age=NULL){private$age = age},
                      set_mate = function(mate=NULL){private$mate = mate},
                      set_genotype = function(genotype=NULL){private$genotype = genotype},

                      age_one_day = function() {private$age = private$age + 1},
                      print_info = function(){file.path(private$age, private$genotype, private$mate, fsep = "\t")},
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
