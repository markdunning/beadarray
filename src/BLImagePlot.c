#include <math.h>
#include <R.h>
#include "beadarray.h"

void BLImagePlot(int *nProbes, double *fground, double *yvalues, int *ygrid, double *result, int *nrow){

    int i,j;
	double counter, total;
//	Rprintf("%d\n", (*nrow));
    for(i = 0; i < (*nrow); i++){
	counter = 0;
	total = 0;
//	Rprintf("%d\n", i);
//	Rprintf("Ygrid[i]: %d Ygrid[i+1]: %d\n", ygrid[i], ygrid[i+1]) ;
        for(j = 0; j < (*nProbes); j++){
//            if(yvalues[j] > ygrid[i]){
//			    Rprintf("yvalue: %f, ygrid: %d, j: %d\n", yvalues[j], ygrid[i], j);
//			}
            if((yvalues[j] >= ygrid[i]) & (yvalues[j] < ygrid[i+1])){
//			    Rprintf("here\n");
                total = total + log2(fground[j]);
				counter = counter + 1;
			}
		}
		result[i] = total/counter;
//		counts[i] = counter;
	}
}
