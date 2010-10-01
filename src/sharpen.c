#include "beadarray.h"

SEXP illuminaSharpen(SEXP pixelMatrix) {

  SEXP sharpened;
  int imageWidth, imageHeight, i, j, sum;
  int tid, nthreads, chunk = 10000;
  
  #if defined (_OPENMP)
  omp_set_num_threads(2);
  #endif
  
  imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
  imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

  PROTECT(sharpened = allocMatrix(REALSXP, imageHeight, imageWidth));


  #pragma omp parallel shared(nthreads,chunk) private(i, j, sum) 
  {
    #pragma omp for schedule(dynamic,chunk)
    for(i = 0; i < imageHeight; i++) {
        for(j = 0; j < imageWidth; j++) {
            REAL(sharpened)[i + (imageHeight * j)] = INTEGER(pixelMatrix)[i + (imageHeight * j)];
        }
    }
  
    #pragma omp for schedule(dynamic,chunk)
    for(i = 1; i < imageHeight - 1; i++) {
        for(j = 1; j < imageWidth - 1; j++) {     
            sum = INTEGER(pixelMatrix)[i + (imageHeight * (j - 1))] + INTEGER(pixelMatrix)[i + (imageHeight * j) - 1] + INTEGER(pixelMatrix)[i + (imageHeight * (j + 1))] + INTEGER(pixelMatrix)[i + (imageHeight * j) + 1];
      
            REAL(sharpened)[i + (imageHeight * j)] = INTEGER(pixelMatrix)[i + (imageHeight * j)] - (0.5 * (sum - 4 * INTEGER(pixelMatrix)[i + (imageHeight * j)]));
      
        }
    }
  }
   
  UNPROTECT(1);
  return(sharpened);
}
