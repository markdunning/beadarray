#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#if defined (_OPENMP)
    #include <omp.h>
#endif

#define min(X, Y) ((X) < (Y) ? (X) : (Y))

#define TRUE 1
#define FALSE 0
