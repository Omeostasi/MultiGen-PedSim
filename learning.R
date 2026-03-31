setClass("Person", 
         slots = c(
           age = "numeric",
         name = "character"))

setGeneric("age", function(x) standardGeneric("age"))

age(Giovanni)

setGeneric("age<-", function(x, value) standardGeneric("age<-"))

setMethod("age", "Person", function(x) x@age)
setMethod("age<-", "Person", function(x,value){
  x@age <- value
  x
})

age(Giovanni) <- 50
age(Giovanni)


# prototype

setClass("Person",
         slots = c(
           age = "numeric",
           name = "character"
         ),
         prototype = list(
           name = NA_character_,
           age = NA_real_
         )
         )

me <- new("Person", name = "Giovanni")
str(Giovanni)

