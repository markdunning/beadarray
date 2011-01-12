#include "beadarray.h"

/* background */
 
void backgroundCalc(int start, int end, int nbeads, int imageWidth, int imageHeight, SEXP pixelMatrix, SEXP coords, double *bg, int tid, int median) {
    
    double x, y, newX, newY;
    int i, j, k, count, tmp, M[289];
    
    for(i = start; i < end; i++) {
    
        x = REAL(coords)[i];
        y = REAL(coords)[i + nbeads];

        newX = floor(x);
        newY = floor(y);

        if(newX == x) 
            newX--;
        if(newY == y)
            newY--;
        
        count = 0;

        int start1, start2, end1, end2;
        start1 = (int) newX - 8;
        start2 = (int) newY - 8;
        end1 = (int) newX + 8;
        end2 = (int) newY + 8;

        /* catch cases where we're too close to the edge of the image */
        if( (start1 < 0) || (start2 < 0) || (end1 > (imageWidth - 1)) || (y > (imageHeight - 1)) ) {
            //Rprintf("Coordinate out of bounds\n");
            bg[i] = R_NaReal;
        }
        else {
            for(j = start1; j <= end1; j++) {
                tmp = j * imageHeight;
                for(k = start2; k <= end2; k++ ) {
                    M[count++] = INTEGER(pixelMatrix)[tmp + k];
                }
            } 

            /* use R's internal sorting algorithm */
            R_qsort_int(M, 1, 289);
        
            /* should the mean or median of the 5 dimmest be used? */
            if(median) {
                bg[i] = M[2];
            }
            else {
                bg[i] = (M[0] + M[1] + M[2] + M[3] + M[4]) / 5.0;
            }
        }
    }
} 

SEXP illuminaBackground(SEXP pixelMatrix, SEXP coords) {

    SEXP background;
    int imageWidth, imageHeight, nbeads, i;
    double *bg;
    int start, end; 
    int tid, nthreads;
        
    imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
    imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

    nbeads = INTEGER(getAttrib(coords, R_DimSymbol))[0];
    PROTECT(background = allocVector(REALSXP, nbeads));
    bg = REAL(background);

    /* find the number of processors and set the number of threads 
    currently has a maximum value of 6 */
    #if defined (_OPENMP)
    omp_set_num_threads( omp_get_num_procs() );
    #endif
    
    /* initialize the background vector with zeros */
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
        
        backgroundCalc(start, end, nbeads, imageWidth, imageHeight, pixelMatrix, coords, bg, tid, 0);

    }
    UNPROTECT(1);
    return(background);
} 

SEXP medianBackground(SEXP pixelMatrix, SEXP coords) {

    SEXP background;
    int imageWidth, imageHeight, nbeads, i;
    double *bg;
    int start, end; 
    int tid, nthreads;
        
    imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
    imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

    nbeads = INTEGER(getAttrib(coords, R_DimSymbol))[0];
    PROTECT(background = allocVector(REALSXP, nbeads));
    bg = REAL(background);

    /* find the number of processors and set the number of threads  */
    #if defined (_OPENMP)
    omp_set_num_threads( omp_get_num_procs() );
    #endif
    
    /* initialize the background vector with zeros */
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
        
        backgroundCalc(start, end, nbeads, imageWidth, imageHeight, pixelMatrix, coords, bg, tid, 1);

    }
    UNPROTECT(1);
    return(background);
} 


/* foreground */

double matrixMean(SEXP pixelMatrix, int imageHeight, int x, int y) {

  int i, j;
  double result = 0.0;

  for(i = x-1; i <= x+1; i++ ) {
    for(j = y-1; j <= y+1; j++ ) {
      result += REAL(pixelMatrix)[(i * imageHeight) + j];
    }
  }
  return(result / 9.0);
}


SEXP illuminaForeground(SEXP pixelMatrix, SEXP coords) {

    SEXP foreground;
    int imageWidth, imageHeight, nbeads, i;
    double x, y, xc, yc, av[4], w[4], *fg;
    
    /* dimensions of the image */
    imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
    imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];
    /* number of beads we have centres for */
    nbeads = INTEGER(getAttrib(coords, R_DimSymbol))[0];

    PROTECT(foreground = allocVector(REALSXP, nbeads));
    fg = REAL(foreground);
    
    for(i = 0; i < nbeads; i++) {

        x = REAL(coords)[i];
        y = REAL(coords)[i + nbeads];
        
        /* catch cases where coordinates are outside the image dimensions */
        if( (x < 1) || (y < 1) || (x > (imageWidth - 1)) || (y > (imageHeight - 1)) ) {
            //Rprintf("Coordinate out of bounds\n");
            fg[i] = R_NaReal;
        }
        else {
            xc = x - floor(x);
            yc = y - floor(y);

            av[0] = (matrixMean(pixelMatrix, imageHeight, (int) floor(x), (int) floor(y)));
            av[1] = (matrixMean(pixelMatrix, imageHeight, (int) floor(x), (int) floor(y+1)));
            av[2] = (matrixMean(pixelMatrix, imageHeight, (int) floor(x+1), (int) floor(y+1)));
            av[3] = (matrixMean(pixelMatrix, imageHeight, (int) floor(x+1), (int) floor(y)));
            
            w[0] = ((1 - xc) * (1 - yc));
            w[1] = ((1 - xc) * yc);
            w[2] = (xc * yc);
            w[3] = (xc * (1 - yc));

            fg[i] = ((w[0] * av[0]) + (w[1] * av[1]) + (w[2] * av[2]) + (w[3] * av[3]));
        }

    }

    UNPROTECT(1);
    return(foreground);
}

/* sharpening */

SEXP illuminaSharpen(SEXP pixelMatrix) {

  SEXP sharpened;
  int imageWidth, imageHeight, i, j, sum;

  imageHeight = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[0];
  imageWidth = INTEGER(getAttrib(pixelMatrix, R_DimSymbol))[1];

  PROTECT(sharpened = allocMatrix(REALSXP, imageHeight, imageWidth));

    #pragma omp parallel for private(i, j) shared(sharpened) num_threads(2)
    for(i = 0; i < imageHeight; i++) {
        for(j = 0; j < imageWidth; j++) {
            REAL(sharpened)[i + (imageHeight * j)] = INTEGER(pixelMatrix)[i + (imageHeight * j)];
        }
    }
  
    #pragma omp parallel for private(i, j, sum) shared(sharpened) num_threads(2)
    for(i = 1; i < imageHeight - 1; i++) {
        for(j = 1; j < imageWidth - 1; j++) {     
            sum = INTEGER(pixelMatrix)[i + (imageHeight * (j - 1))] + INTEGER(pixelMatrix)[i + (imageHeight * j) - 1] + INTEGER(pixelMatrix)[i + (imageHeight * (j + 1))] + INTEGER(pixelMatrix)[i + (imageHeight * j) + 1];
      
            REAL(sharpened)[i + (imageHeight * j)] = INTEGER(pixelMatrix)[i + (imageHeight * j)] - (0.5 * (sum - 4 * INTEGER(pixelMatrix)[i + (imageHeight * j)]));
      
        }
    }
   
  UNPROTECT(1);
  return(sharpened);
}

