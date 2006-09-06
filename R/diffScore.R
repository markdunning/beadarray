#example usage to calculate DiffScore of array 1 (ref) vs array 2 (cond)
#df = IlluminaDE(BSData$R, QC$AVG.Signal[,11], BSData$BeadStDev, QC$Var[,11], 2,1)


##Function is work in progress!


##E is the expression matrix from the BeadStudio output file
##E.negs is the Signal columns from the qc file
##s contains the BeadStDev values
##s.negs is the Var columns from the qc value
##cond is the indices of the condition arrays
##ref is the indices of the reference arrays



IlluminaDE = function(E, E.negs, s, s.negs,cond,ref){

  ##Find the columns relating to the condition array
  ##I.cond is the expression values
  ##s.cond is the standard deviation values
  ##sigma_neg.ref is the standard deviations of the negative controls

  
    I.cond =E[,cond]
    s.cond =s[,cond]
    sigma_neg.cond=s.negs[cond]                                                                                                                                                                                                                                                                                                                                                                                                              =s.negs[cond]

    ##Find the columns relating to the reference array
    ##I.cond is the expression values
    ##s.cond is the standard deviation values
    ##sigma_neg.ref is the standard deviations of the negative controls


    
    I.ref=E[,ref]
    s.ref=s[,ref]
    sigma_neg.ref=s.negs[ref]
    

    ##Fit robust linear model
    
   lm.cond = rlm(s.cond~I.cond)

  lm.ref = rlm(s.ref~I.ref)

    ##Get coefficients a and b
    
  a.ref = lm.ref$coefficients[1]

  b.ref = lm.ref$coefficients[2]

  a.cond = lm.cond$coefficients[1]

  b.cond = lm.cond$coefficients[2]

    
    
    
  tech.ref = a.ref + (2.5*b.ref*I.ref)

  tech.cond = a.cond + (2.5*b.cond*I.cond)

  denom = sqrt(tech.ref^2 + tech.cond^2 + sigma_neg.ref^2 + sigma_neg.cond^2)


  I = abs(I.cond - I.ref)
  
  p = 2*pnorm(I /denom, lower.tail=FALSE)


  
  Diffscore = 10 * sign(I.cond-I.ref)*log10(p)
  
  Diffscore

}





