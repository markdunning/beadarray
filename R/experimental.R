



plotProbe <- function(symbol,rng, tx){
  data(genesymbol)  
  p1 <- autoplot(tx, which=genesymbol[symbol])
  p2 <- rng[which(rng %over% genesymbol[symbol])]
  tracks(p1, autoplot(p2,aes(fill=PROBEQUALITY)))
  
}



