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

