#example usage to calculate DiffScore of array 1 (ref) vs array 2 (cond)
#df = IlluminaDE(BSData$R, QC$AVG.Signal[,11], BSData$BeadStDev, QC$Var[,11], 2,1)


##Function is work in progress!


IlluminaDE = function(E, E.negs, s, s.negs,cond, ref){


  if(length(cond) >1){
  
  I.cond = apply(E[,cond], 1, mean)

  s.cond = apply(E[,cond], 1, sd)
  s2.cond = var(E.negs[ref])
  }

  else{

    I.cond =E[,cond]
    s.cond =s[,cond]
    s2.cond =s.negs[cond]

    
    
  }

  if(length(ref) >1){
    
  I.ref = apply(E[,ref], 1, mean)

  s.ref = apply(E[,ref], 1, sd)
  s2.ref = var(E.negs[ref])
}

  else{

    I.ref=E[,ref]
    s.ref=s[,ref]
    s2.ref=s.negs[ref]
    
  }
  
   lm.cond = rlm(s.cond~I.cond)

  lm.ref = rlm(s.ref~I.ref)

  a.ref = lm.ref$coefficients[1]

  b.ref = lm.ref$coefficients[2]

  a.cond = lm.cond$coefficients[1]

  b.cond = lm.cond$coefficients[2]

  tech.ref = a.ref + (2.5*b.ref*I.ref)

  tech.cond = a.cond + (2.5*b.cond*I.cond)



  denom = sqrt(tech.ref^2 + tech.cond^2 + s2.ref + s2.cond)


  I = abs(I.cond - I.ref)
  
  p = 2*pnorm(I /denom, lower.tail=FALSE)


  mu.cond = mean(I.cond)
  mu.ref = mean(I.ref)

  Diffscore = 10 * sign(I.ref-I.cond)*log10(p)
  
  Diffscore

}





