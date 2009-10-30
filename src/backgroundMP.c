#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <omp.h>

void performCalc(int start, int end, int nbeads, int imageHeight, SEXP pixelMatrix, SEXP coords, double *bg) {
    
    double x, y, newX, newY;
    int i, j, k, count, M[289], tmp;
    
    for(i = start; i < end; i++) {
    
    x = REAL(coords)[i];
    y = REAL(coords)[i + nbeads];

    newX = floor(x);
    newY = floor(y);

    if(newX == x) {
      newX--;
     }
    if(newY == y)
      newY--;
    
    count = 0;

    int start1, start2, end1, end2;
    start1 = (int) newX - 8;
    start2 = (int) newY - 8;
    end1 = (int) newX + 8;
    end2 = (int) newY + 8;

    for(j = start1; j <= end1; j++) {
      tmp = j * imageHeight;
       for(k = start2; k <= end2; k++ ) {
        M[count++] = INTEGER(pixelMatrix)[tmp + k];
      }
    }

    /* use R's internal sorting algorithm */
    R_qsort_int(M, 1, 289);

    bg[i] = (M[0] + M[1] + M[2] + M[3] + M[4]) / 5.0;
    }
}
/*
SEXP backgroundMP(SEXP pixelMatrix, SEXP coords) {

  SEXP background;
  int imageWidth, imageHeight, nbeads, i, j, k;
  double x, y, newX, newY, *bg;
  int tmp; 
  int nthreads, tid;
    
  imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
  imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

  nbeads = INTEGER(getAttrib(coords, R_DimSymbol))[0];
  PROTECT(background = allocVector(REALSXP, nbeads));
  bg = REAL(background);

  omp_set_num_threads(2);
 
  #pragma omp parallel default(shared) private(nthreads, tid)
  {
    Rprintf("Hello World from thread = %d\n", tid);
    if (tid == 0) {
        performCalc(0, nbeads/2, nbeads, imageHeight, pixelMatrix, coords, bg);
    }
    else {
        performCalc(nbeads/2, nbeads, nbeads, imageHeight, pixelMatrix, coords, bg);
    }
  }
  UNPROTECT(1);
  return(background);
} */

void silly(int i) {
 
    int nthreads, tid;
    
//  #pragma omp parallel default(shared) private(nthreads, tid)
//{
    Rprintf("Hello World from thread = %d\n", tid);
//  }
}

SEXP backgroundMP(SEXP pixelMatrix, SEXP coords) {

  SEXP background;
  int imageWidth, imageHeight, nbeads, i, j, k;
  double x, y, newX, newY, *bg;
  int tmp; 
  int nthreads, tid;
    
  imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
  imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

  nbeads = INTEGER(getAttrib(coords, R_DimSymbol))[0];

  omp_set_num_threads(2);
  PROTECT(background = allocVector(REALSXP, nbeads));
  bg = REAL(background);
  for(i = 0; i < nbeads; i++) { bg[i] == 0; }

  //  silly(4);

  UNPROTECT(1);
  return(background);
} 
