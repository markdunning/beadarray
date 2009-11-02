#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <omp.h>

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

void performCalc(int start, int end, int nbeads, int imageHeight, SEXP pixelMatrix, SEXP coords, double *bg, int tid) {
    
    double x, y, newX, newY;
    int i, j, k, count, tmp, M[289];
    
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

SEXP backgroundMP(SEXP pixelMatrix, SEXP coords) {

  SEXP background;
  int imageWidth, imageHeight, nbeads, i,;
  double *bg;
  int start, end; 
  int tid, nthreads, numProcs;
    
  imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
  imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

  nbeads = INTEGER(getAttrib(coords, R_DimSymbol))[0];
  PROTECT(background = allocVector(REALSXP, nbeads));
  bg = REAL(background);

  /* find the number of processors and set the number of threads 
  currently has a maximum value of 4 */
  #if defined (_OPENMP)
  numProcs = min(omp_get_num_procs(), 4);
  omp_set_num_threads(numProcs);
  #endif
  
  /* initialize the background vector with zeros, can probably be removed */
  for(i = 0; i < nbeads; i++) { bg[i] = 0; }
  
  #pragma omp parallel shared(nthreads, nbeads, imageHeight, pixelMatrix, coords) private(tid, start, end)
  {
      
      #if defined (_OPENMP)
      nthreads = omp_get_num_threads();
      tid = omp_get_thread_num();
      #else
      nthreads = 1;
      tid = 0;
      #endif
      
    start = (int) floor(tid * ((double)nbeads / (double)nthreads));
    end = (int) floor((tid + 1) * ((double)nbeads / (double)nthreads));
    Rprintf("%d\t%d\t%d\t%d\t%d\n", tid, nbeads, nthreads, start, end);
    
    performCalc(start, end, nbeads, imageHeight, pixelMatrix, coords, bg, tid);

  }
  UNPROTECT(1);
  return(background);
} 



