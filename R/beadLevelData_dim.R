## Returns the dimensions of a beadLeveLData object
## First element is the number of arrays
## Subsequent elements are the number of beads on each array

setMethod("dim", "beadLevelData", function(x) {
    an = sectionNames(x)
    nArrays = length(an)
    nBeads = numBeads(x)
    c("nArrays"=nArrays, "nBeads"=nBeads)
 } )

