#include <R.h>

void findIndices(int *probeID, int *probeList, int *nrow, int *indices, int *start){

  unsigned int j = 0;
  unsigned int i = (*start)-1;
  unsigned int probe = *probeID;

  while(probeList[i] != probe){
    i++;
  }

  while(probeList[i] == probe){
    indices[j] = i+1;
    j++;
    i++;
  }

  *start = i;

  //  Rprintf("%d\n", i);

  //  for(i=0; i<*nrow; i++){
  //  if(probeList[i] == *probeID){
  //    indices[j] = i+1;
  //  j++;
  //  }
  //}
}
