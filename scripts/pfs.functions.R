# August 2016 - moved all functions to individual .R files in R/


# assuming setwd('C:/Users/Dave/HalfStarted/predictFromSequence/') already done

# http://stackoverflow.com/questions/10291520/reading-all-scripts-and-data-files-from-multiple-folders
file.sources = list.files("R", 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)

sapply(file.sources,source,.GlobalEnv)


# December 2015 
# Decided that base objects are not use-able for whole chromosomes.
# Need to re-write using hits objects

# function 
# description
# Parameters:-
# factorSet
#createBindingFactor <-  function()  {
#  bindingFactor <- list(name=name, type=type, 
#                        profile=profileList, 
 #                       mods=modList)
#  
#  return(bindingFactor)
#}




