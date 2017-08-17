#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "beadarray.h"

void HULK(double* residuals, int* neighbours, int* nbeads, int* invasions, double *results) {

  int nextslot, firstinv, lastinv;
  int temp, temp2;
  int i, j, k, l;
  int invaded[50000];
  double divide, weight;

  int *done = (int *) R_alloc(*nbeads, sizeof(int));
  memset(done, 0, *nbeads * sizeof(int));

  //for each node...
  for(i = 0; i < *nbeads; i++) {
    nextslot = 1; 
    done[i] = 1;
    invaded[0] = i;
    firstinv = 0;
    lastinv = 0;
    divide = 1.0;
      
    //invasion process for this node
    //----------------
    //invade (invasions) times please
    for(l = 1; l <= *invasions; l++)	{

      weight = 1.0 / (double) (l + 1);
      //weight = 1.0;

      //for each invaded node...
      for(j = firstinv; j <= lastinv; j++) {
	temp = 6*invaded[j]; //first candidate neighbour (cannot pass arrays, hence using a vector form)
	//check its neighbours
	for(k = 0; k <= 5; k++) {
	  //did we invade? if not, let's do it (remembering to avoid 0 slots)
	  temp2 = neighbours[temp]-1; //vector form for neighbour to invade
	  if((temp2 >= 0) && (!done[temp2])){
	    invaded[nextslot] = temp2;
            /* if the residual is > 0 then include it */
	    if(residuals[temp2]) {
	      results[i]  += (residuals[temp2] * weight);
	      divide += weight;
	    }
	    done[temp2] = 1;
	    nextslot++;
	  }
	  temp++;
	}
      }
      firstinv = lastinv + 1;
      lastinv = nextslot - 1;
    }
    
    results[i] /= divide;

    //clear done flags
    for(j = 0; j < nextslot ; j++) {
    	done[invaded[j]] = 0;
    }
  }
  return;
}
