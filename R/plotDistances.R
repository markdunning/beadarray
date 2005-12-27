"plotDistances" <-
function(coords){

d = dist(coords)

hist(dist(coords), freq=FALSE,main="")

dens=dnorm(1:1500, mean=650, sd=350)


abline(v=mean(d), col="red")

abline(v=650, col="blue")

lines(1:1500, dens, col="blue")

}

