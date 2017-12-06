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
                      initialize = function(genotype=NULL, age=NULL, sex=NULL){
                        private$age = age
                        private$sex = sex
                        private$mate = NULL
                        private$genotype = genotype
                        
                      }, # end constructor
                      
                      # setters
                      set_age = function(age=NULL){private$age = age},
                      set_sex = function(sex=NULL){private$sex = sex},
                      set_mate = function(mate=NULL){private$mate = mate},
                      set_genotype = function(genotype=NULL){private$genotype = genotype},
                      
                      hitched = function() Mate(private),
                      grow_up = function() GrowUp(private),
                      
                      #getters
                      get_age = function(){return(private$age)},
                      get_sex = function(){return(private$sex)},
                      get_mate = function(){return(private$mate)},
                      get_genotype = function(){return(private$genotype)}

                    ), # end public
                    
                    private = list(
                      
                      # fields
                      age = NULL,
                      sex = NULL,
                      mate = NULL,
                      genotype = NULL

                    ) # end private
                    
) # end class definition