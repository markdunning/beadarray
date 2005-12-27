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

  if(is.null(BLData$x)){
    stop("X and Y co-ordinates are not present in this data")
  }


#arrx and ary take all x and y values for beads in the array

arrx <- BLData$x[,array]
arrx<-as.integer(arrx)
arry <- BLData$y[,array]
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

