#include "beadarray.h"


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
    int imageWidth, imageHeight, nbeads, i, j,  col, row;
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

    UNPROTECT(1);
    return(foreground);
}
