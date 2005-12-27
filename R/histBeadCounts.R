"histBeadCounts" <-
function(BLData,array, main=NULL){

#Finds the number of each probe / bead used on the array and plots this.
#In theory this should show a poisson distribution.

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

  probes = sort(unique(BLData$probeID[BLData$probeID[,1] > 0,1]))

counts=vector(length = length(probes))

for(i in 1:length(probes)){

counts[i] = length(BLData$R[BLData$probeID[,array]==probes[i],array])

}


x = seq(1:50)


#noprobes is number of distinct probes found on the array. We discount
#probes with ID 0, those with negative ID and probes with ID 5244 [which
#appears over a 1000 times on each array].


noprobes = length(probes)

#Here lambda is theoretical mean of poisson trials. Approximation to binomial with mean~=50,000 and prob. of success ~= 1/1500

lambda = sum( (BLData$probeID[,array]!=5244) & (BLData$probeID[,array]>0) ) / noprobes



hist(counts[counts>0 & counts<1000], nclass=100, ylim=range(0,120), xlim=range(0,60), xlab="Number of beads", ylab="Frequency", main=main)

#Draw the theoretical Poisson distribution and vertical lines at the expected and observed means.

dens=dpois(0:60, lambda)
lines(0:60, dens*noprobes,col="blue")
abline(v=c(median(counts[counts>0 & counts<1000]), lambda),col=c("black", "blue"))

#par(cex=0.9)
#legend(40,100, c("Observed", "Expected"), c("black", "blue"))


}

