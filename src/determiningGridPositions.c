#include "beadarray.h"

/* function to round the values from the .locs file into the truncated versions that 
are found in the .txt file 
SEXP roundLocsFileValues(SEXP inputVector) {
    
    int i, vecLength, digits;
    double x;
    
    for(i = 0; i < length(inputVector); i++) {
        x = REAL(inputVector)[i];
        // the precison of the rounding is determined by the integer part of the value 
        if(x >= 10000)
            digits = 2;
        else if (x >= 1000)
            digits = 3;
        else if (x >= 100)
            digits = 4;
        else if (x >= 10)
            digits = 5;
        else
            digits = 6;
        // perform the rounding to the required precision 
        REAL(inputVector)[i] = round(REAL(inputVector)[i]*(pow(10,digits)))/pow(10,digits);
    }
    return(inputVector);
} */

/* takes the list of indicies within the locs file and returns a two column matrix,
with each row the x&y positions within the array grid
params arguments is a vector of length 3 containing the number of rows and columns in a 
segments, and the size of the gap to insert between sgements */
SEXP locsIndicesToGrid(SEXP inputVector, SEXP params) {

    int i, vecLength;
    int nrow, ncol, gap, segIdx;
    int idx, col, row;
    SEXP res;

    vecLength = length(inputVector);
    nrow = INTEGER(params)[0];
    ncol = INTEGER(params)[1];
    gap = INTEGER(params)[2];

    /* allocate the matrix that will be returned */
    PROTECT(res = allocMatrix(INTSXP, vecLength, 2));

    for(i = 0; i < vecLength; i++) {
        idx = INTEGER(inputVector)[i] - 1;
        /* which segment are we in */
        segIdx = idx / (nrow * ncol);

        /* Calculate the column */
        col = abs( ( ( idx % (ncol * nrow) ) / nrow) - ncol );
        
        /* Calculate row */
        row = 2 * (idx % nrow) + 1;
        /* adjust for the grid */
        row = row + ( segIdx * ( (nrow * 2) + gap ) );
        /* if in an odd colum move it up by one */
        if(col % 2 == 0) 
            row = row + 1;

        INTEGER(res)[i] = row;
        INTEGER(res)[i + vecLength] = col;
    }
    UNPROTECT(1);
    return(res);
}


#define TOOSMALL(A) (A <= 0 ? NA_INTEGER : A)
#define TOOLARGE(A, B) (A > B ? NA_INTEGER : A)

SEXP neighboursFromLocs(SEXP params) {

    int i, vecLength;
    int nrow, ncol;
    int idx, col, row;
    SEXP res;
    int tmpInt;


    nrow = INTEGER(params)[0];
    ncol = INTEGER(params)[1];
    vecLength = nrow * ncol;

    /* allocate the matrix that will be returned */
    PROTECT(res = allocMatrix(INTSXP, vecLength, 7));
    int *r = INTEGER(res);

    for(i = 0; i < vecLength; i++) {
        idx = i;

        /* Calculate the column */
        col = abs( ( ( idx % (ncol * nrow) ) / nrow) - ncol );
        
        /* Calculate row */
        row = 2 * (idx % nrow) + 1;

        /* if in an odd colum move it up by one */
        if(col % 2 == 0) 
            row = row + 1;

        tmpInt = 2 * pow(0, col % 2);
        r[i] = i + 1;
        r[i + vecLength] = TOOSMALL(i - nrow + 1);
        r[i + (2 * vecLength)] = TOOSMALL(i - nrow + tmpInt);
        r[i + (3 * vecLength)] = TOOSMALL(i);
        r[i + (4 * vecLength)] = TOOLARGE(i + 2, vecLength);
        r[i + (5 * vecLength)] = TOOLARGE(i + nrow + 1, vecLength);
        r[i + (6 * vecLength)] = TOOLARGE(i + nrow + tmpInt, vecLength);
        
        /* deal with the cases in the top and bottom rows */
        if(row == 1) { // bottom, odd column
            r[i + (2 * vecLength)] = NA_INTEGER;
            r[i + (3 * vecLength)] = NA_INTEGER;
            r[i + (6 * vecLength)] = NA_INTEGER;
        }
        else if(row == 2) { // bottom, even column
            r[i + (3 * vecLength)] = NA_INTEGER;
        }
        else if(row == (2 * nrow) ) { // top, even column
            r[i + (4 * vecLength)] = NA_INTEGER;
        }
        else if(row == (2 * nrow - 1)) { // top, odd column
            r[i + (2 * vecLength)] = NA_INTEGER;
            r[i + (4 * vecLength)] = NA_INTEGER;
            r[i + (6 * vecLength)] = NA_INTEGER;
        }
        
    }
    UNPROTECT(1);
    return(res);
}


SEXP hashLocsToTxtIndices(SEXP neighboursMat, SEXP hash) {
    
    int i, tmp;
    int nrow = INTEGER(getAttrib(neighboursMat, R_DimSymbol))[0];
    int ncol = INTEGER(getAttrib(neighboursMat, R_DimSymbol))[1];
    int len = nrow * ncol;
    
    double *mat = REAL(neighboursMat);
    double *h = REAL(hash);
    
    SEXP res;
 
    /* allocate the matrix that will be returned */
    PROTECT(res = allocMatrix(REALSXP, nrow, ncol));
    double *r = REAL(res);
    
    for(i = 0; i < len; i++) {
        if(!ISNA(mat[i])) {
            tmp = (int) mat[i] - 1;
            r[i] = h[tmp];
        }
        else {
            r[i] = NA_REAL;
        }
    }
    
    UNPROTECT(1);
    return(res);
}