DiffScore = function(BSData, QC, cond, ref){

  Avg_Intensity_Ref = exprs(BSData)[,ref]
  Avg_Intensity_Cond = exprs(BSData)[,cond]

  BeadStDev_Ref = se.exprs(BSData)[,ref]
  BeadStDev_Cond = se.exprs(BSData)[,cond]

  SigmaRef = 2.5*QC$Var[ref,11]
  SigmaCond = 2.5*QC$Var[cond,11]

  lm.ref = rlm(BeadStDev_Ref~Avg_Intensity_Ref)
  lm.cond = rlm(BeadStDev_Cond~Avg_Intensity_Cond)

  Aref = lm.ref$coefficients[1]
  Bref = lm.ref$coefficients[2]

  Acond = lm.cond$coefficients[1]
  Bcond = lm.cond$coefficients[2]
  
  

Denom1 = (Aref +  (2.5*Bref*Avg_Intensity_Ref))^2 + SigmaRef^2

Denom2 = (Acond +  (2.5*Bcond*Avg_Intensity_Cond))^2 + SigmaCond^2

FullDenom = sqrt(Denom1 + Denom2)

Numerator = abs(Avg_Intensity_Cond - Avg_Intensity_Ref)

Diff_Pvalue = 2*pnorm(Numerator/FullDenom, lower.tail=FALSE)

DiffScore = -10*sign(Avg_Intensity_Cond - Avg_Intensity_Ref)*log10(Diff_Pvalue)

DiffScore

}
