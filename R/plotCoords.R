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

