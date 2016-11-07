

# General mutating function. 
# Try to alter objects by supplying list of named mutable features and list of ways to mutate them 
#   new values could be vectors, ranges or functions.
# c.f. mutateFactorSet()  , which does insertion/deleation/duplication/rearrangement of sets of factors.
# May need to be replaced by overloaded generic function later


# 2016-11-07 This script was the development bed for the following functions:-
# - recurseListNames()
# - getListNodeByCharacter()
# - assignListNodeByCharacter()
# - mutateBindingFactor()
# - buildMutationSkeleton()


names(XFS)
str(XFS)


# could parse a list of mutable features.  Except only the tips should be mutable?

# Some binding factors carry functions, how to mutate these?  By parameters?  
# parameters need to be within genomic or chromosomal sensible values - or trimmed when applied?
# fairly academic as it is unlikely that any molecular system 'counts'  distances more than several Kb. 
# at that point, interactions are based on larger scale domains and interactions.

# beta is very versatile https://en.wikipedia.org/wiki/Beta_distribution
# but maybe should be more explictit


# need some method to both store a function and mutate the parameters



mut.list <- list()

# but must be able to allocate mutations evenly (or proportionately) to each mutation type
#  therefore, need to be able to count how many valid ways a give binding factor may be mutated.
# One problem is in traversing the tree like structure of a list  to find tips.
#   How many tips in a list.

test.list <- list(a=54, b=c(23,24), c=list(sub.a=34, sub.b=list(subSub.a=10, subSub.b=c(23,19))), d=300)

names(test.list[[3]][[2]])
length(unlist(test.list))

4

mutFuncList <- list(a=function() round(runif(1,min=1, max=100)),
                    c=list(
                        sub.b=list(subSub.b=function() sort(round(runif(2,min=1, max=100))) )
                    ))
                    
length(unlist(mutFuncList))

# unlisted names are not suitable for parsing the list as different configurations could lead to the same unlisted name
# e.g.   a.b$c  and a$b.c both become a.b.c


names(mutFuncList)
str(mutFuncList)


vecList <- recurseListNames(mutFuncList)


recurseListNames(XFS)


#for(i in 

#rapply(mutFuncList, how="replace", f = function(x, stem) paste(stem,names(x))   , stem="STEM")  # Nope
 

# OK, using recurseListNames can get string reference to terminal nodes of an arbitrary list.
# Now need to be able to go in and change the value of a list based on this refernce



getListNodeByCharacter(mutFuncList, vecList[2])



get("mutFuncList[['c']][['sub.b']][['subSub.b']]")  # nope
get("mutFuncList$c$sub.b$subSub.b")                 # nope
eval("mutFuncList[['c']][['sub.b']][['subSub.b']]")  # nope
eval("mutFuncList$c$sub.b$subSub.b") 
get("mutFuncList")["[['c']][['sub.b']][['subSub.b']]"] 
get("mutFuncList['c']['sub.b']['subSub.b']")  # nope    incase [] differed from [[]]

# getListNodeByCharacter() is great, but what I really need to be able to do is assign 
# a new value to the relevant list object !

# The reverse (or "inverse") of a <- get(nam) is assign(nam, a), assigning a to name nam.
# proabbly won't help if can't use get.


# This may be currently impossible in R

# apparently there is already recursive indexing!
# http://stackoverflow.com/questions/1169456/in-r-what-is-the-difference-between-the-and-notations-for-accessing-the
mutFuncList[[c(2,1,1)]]    # Yes!
# can this be assigned?   
mutFuncList[[c(2,1,1)]] <- function(x)  print(x)     # yes!

mutFuncList[[c("c", "sub.b")]]


mutFuncList[[unlist(strsplit(vecList[2], split=","))]]    # boom!


newFuncList <- assignListNodeByCharacter(mutFuncList, function() return(500), charRef=vecList[2])



# OK should now have most parts to build a function to mutate arbitrary lists.



# given a bindingFactor, create a partially matched list with functions to 
# functions that take an argument x should alter the existing value  TODO - not yet implemented
# functions with no argument should generate brand new values from a list or distribution
# new values can themselves be functions...






newList <- mutateBindingFactor(test.list, mutFuncList, verbose=TRUE)
newList <- mutateBindingFactor(test.list, mutFuncList, verbose=TRUE, n.muts=10)
newList
test.list
identical(newList, test.list)



# build mutation model from binding factor

bf.chars <-recurseListNames(bf.spreadRep)
 getListNodeByCharacter(bf.spreadRep, "mods,H3K27me3,offset" )

#mut.bf.spreadRep <- list()
#assignListNodeByCharacter(mut.bf.spreadRep, 200, charRef="mods,H3K27me3,offset" )  # nope




buildMutationSkeleton(bf.chars)   # a full empty (NA list)

buildMutationSkeleton(c("mods,H3K27me3,stateWidth", "mods,H3K27me3,offset" ), values=c(20, 50))

mut.bf.spreadRep <-buildMutationSkeleton(c("mods,H3K27me3,stateWidth", "mods,H3K27me3,offset" ), values=c(20, 50))
mut.bf.spreadRep

# replace a set value with a function
mut.bf.spreadRep <-assignListNodeByCharacter(mut.bf.spreadRep, function() round(runif(1,min=50,max=400)), charRef="mods,H3K27me3,stateWidth")
mut.bf.spreadRep <-assignListNodeByCharacter(mut.bf.spreadRep, function() round(runif(1,min=0,max=1000)), charRef="mods,H3K27me3,offset")
# hist(rnbinom(1000, 10, mu= 200))  an alternative skewed distribution centred around 200


bf.spreadRep.mutated <- mutateBindingFactor(bf.spreadRep, mut.bf.spreadRep, n.muts=5,verbose=TRUE)




