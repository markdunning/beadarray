


"chitest" <-
function(observed, expected){

#Applies standard chi-square formula to a set of expected and observed values

chi = sum((observed-expected)^2/expected)

chi

}




"checkRandomness" <-
function(BLData, array,coords){


#Computes the chi-square statistic to see if the spatial distribution of a set
#of co-ordinates is random. For this we use split the hexagon into 8 sections by 
#firstly splitting into 4 equal quadrants and drawing a circle in the centre. The
#expected number of probes in each section is then proportional to the area of the section
#compared to the area of the whole hexagon.

  #Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

#Check that x,y co-ords are available

  if(is.null(BLData$GrnX)){
    stop("X and Y co-ordinates are not present in this data")
  }


#arrx and ary take all x and y values for beads in the array

arrx <- BLData$GrnX[,array]
arrx<-as.integer(arrx)
arry <- BLData$GrnY[,array]
arry<-as.integer(arry)



xcoords = coords[,1]
ycoords = coords[,2]

#xmax and ymax are maximum x and y co-ordinates of beads in the array.

xmax = max(arrx)-min(arrx)
ymax = max(arry)-min(arry)


xcentre = round(xmax/2)
ycentre = round(ymax/2)



#transx and transy are the x and y displacements from the centre of the hexagon of each point

transx = xcoords - xcentre

transy = ycoords - ycentre


r=matrix(nrow=length(xcoords),ncol=1)
theta=matrix(nrow=length(xcoords),ncol=1)

#Calculate radius of each point in the usual way

r = sqrt(transx^2 + transy^2)

#Angle of each point from the centre
theta = atan(transy/transx)

#quads is going to be the number of probes we find in each of the 8 sections. 

#When r < rmax then point is inside the inner circle
#When r > rmax then point is outside the circle
#If transy > 0 then point is in the upper half of the hexagon
#If transy <0 then point is in the lower half of the hexagon

#If theta > 0 and transy > 0 then point is on right side of the hexagon
#If theta < 0 and transy > 0 then point is on left side of the hexagon


quads = vector(length=8)

quads[1] = length(r[r<xmax/4 & theta >0 & transy > 0])
quads[2] = length(r[r>xmax/4 & theta >0 & transy > 0])
quads[3] = length(r[r < xmax/4 & theta < 0 & transy <0 ])
quads[4] = length(r[r > xmax/4 & theta < 0 & transy <0])
quads[5] = length(r[r< xmax/4 & theta > 0 & transy<0 ])
quads[6] = length(r[r> xmax/4 & theta > 0 & transy<0])
quads[7] = length(r[r<xmax/4 & theta<0 & transy>0])
quads[8] = length(r[r>xmax/4 & theta<0 & transy>0])



#hexarea is area of whole hexagon. This is 3/4 of area of a rectangle with width xmax
#and height ymax

hexarea = xmax * ymax * 0.75

#circlearea is area of circle in middle of hexagon radius xmax/4

circlearea = pi * (xmax/4)^2

#carea is area of each quarter circle
carea = circlearea / 4

#harea is area of all points outside the circle

harea = hexarea / 4 - circlearea / 4


#expexcted_counts is the number of probes we would expect to find in each section if the distribution
#was random (i.e no bias towards a particular section)
#length(coords[,1]) is the number of points we have in total and this is multiplied by the proportional
#area of the section


expected_counts = vector(length =8)

expected_counts[1] = carea/hexarea * length(coords[,1])
expected_counts[2] = harea/hexarea * length(coords[,1])
expected_counts[3] = carea/hexarea * length(coords[,1])
expected_counts[4] = harea/hexarea * length(coords[,1])
expected_counts[5] = carea/hexarea * length(coords[,1])
expected_counts[6] = harea/hexarea * length(coords[,1])
expected_counts[7] = carea/hexarea * length(coords[,1])
expected_counts[8] = harea/hexarea * length(coords[,1])


#We now do a chi-square test

chitest(quads, expected_counts)

}

"findHighestChis" <-
function(BLData, array, limit=14){

  #Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

  probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))
  l = NULL
  for(i in 1:length(probes)){
    c = getProbeCoords(BLData, probes[i], array)   
    r = checkRandomness(BLData, array, c)
    if(r > limit){
      l = c(l, probes[i])
    }
  }
  l
}

"findHighestSD" <-
function(BLData, array, limit=1){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

l = NULL

for(i in 1:length(probes)){

int = getProbeIntensities(BLData, probes[i], array)

if(sd(int, na.rm=TRUE)>limit){

l = c(l, probes[i])

}

}

l

}

"findLowestCounts" <-
function(BLData, array,  limit=24){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
  probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

  l = NULL

  for(i in 1:length(probes)){

    c = getProbeCoords(BLData, probes[i], array)

    if (length(c[,1]) < limit){

      l = c(l, probes[i] )

    }

  }

  l

}

"findMostOutliers" <-
function(BLData, array, limit=5){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

l = NULL

for(i in 1:length(probes)){

o = length(findBeadStatus(BLData, probes[i], array))

if(o > limit){

l = c(l, probes[i])

}

}

l

}

"getDimensions" <-
function(BLData, array){

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
#Check that x,y co-ords are available

  if(is.null(BLData$GrnX)){
    stop("X and Y co-ordinates are not present in this data")
  }
  
arrx <- BLData$GrnX[,array]
arrx<-as.integer(arrx)
arry <- BLData$GrnY[,array]
arry<-as.integer(arry)

#xmax and ymax are maximum co-ordinates

xmax = max(arrx)-min(arrx)
ymax = max(arry)-min(arry)

dims = list(length=2)

dims$xmax = xmax
dims$ymax = ymax

dims
}

"getProbeCoords" <-
function(BLData, probe, array){

#Finds the co-ordinates of all probes with a particular ID

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }
  
#Check that x,y co-ords are available, otherwise stop

  if(is.null(BLData$GrnX)){
    stop("X and Y co-ordinates are not present in this data")
  }
  
arrx <- BLData$GrnX[,array]
arrx<-as.integer(arrx)
arry <- BLData$GrnY[,array]
arry<-as.integer(arry)


#standardise the coordinates to start at 0 


xcoords = (BLData$GrnX[BLData$ProbeID[,array]==probe,array]) - min (arrx) + 1
ycoords = (BLData$GrnY[BLData$ProbeID[,array]==probe,array]) - min (arry) + 1

coords = matrix(nrow=length(xcoords),ncol=2)

coords[,1] = xcoords
coords[,2] = ycoords

coords

}

"histBeadCounts" <-
function(BLData,array, main=NULL){

#Finds the number of each probe / bead used on the array and plots this.
#In theory this should show a poisson distribution.

#Check the object is of class BeadLevelList
  if(class(BLData) != "BeadLevelList"){
    stop("BeadLevelList object required!")
  }

  probes = sort(unique(BLData$ProbeID[BLData$ProbeID[,1] > 0,1]))

counts=vector(length = length(probes))

for(i in 1:length(probes)){

counts[i] = length(BLData$G[BLData$ProbeID[,array]==probes[i],array])

}


x = seq(1:50)


#noprobes is number of distinct probes found on the array. We discount
#probes with ID 0, those with negative ID and 

noprobes = length(probes)

#Here lambda is theoretical mean of poisson trials. Approximation to binomial with mean~=50,000 and prob. of success ~= 1/1500

lambda = sum(BLData$ProbeID[,array]>0)  / noprobes



hist(counts[counts>0 & counts<1000], nclass=100, ylim=range(0,120), xlim=range(0,60), xlab="Number of beads", ylab="Frequency", main=main)

#Draw the theoretical Poisson distribution and vertical lines at the expected and observed means.

dens=dpois(0:60, lambda)
lines(0:60, dens*noprobes,col="blue")
abline(v=c(median(counts[counts>0 & counts<1000]), lambda),col=c("black", "blue"))

#par(cex=0.9)
#legend(40,100, c("Observed", "Expected"), c("black", "blue"))


}

"plotCluster" <-
function(coords,main=NULL){

#Produces a clustering of the distances between a set of given co-ordinates. These
#co-ordinates could represent all probes of the same type on an array, say.

#The 'ward' method is used because it seeks to give us clusters that are
#roughly the same size.

#After clustering we draw the rectangle which gives us 8 distinct clusters. Ideally
#we would like these 8 clusters to represent the 8 sections into which we divide the 
#hexagon. The size of each cluster should give an indication of how well the beads are
#spread.

 
par(cex=0.9)

cl=hclust(dist(coords),method="ward")

par(cex.lab=0.3)

plot(cl, main=main)

r = rect.hclust(cl, k=8)

par(cex=0.6)


}

"plotCoords" <-
function(xmax, ymax, coords, label,...){

plotArray(xmax, ymax,...)

#plotArray draws out the hexagonal shape of the array and divides into 8 sections


if(label){



#On actual arrays, the y co-ordinates are arranged with y=0 at the top of the array, 
#whereas in R we plot with y=0 at the bottom of the plot. Therefore we transform the 
#y co-ordinates before plotting

text(coords[,1], (ymax - coords[,2]),labels = seq(1:length(coords[,1])),...)

}

else{

points(coords[,1], (ymax - coords[,2]),...)

}


}

