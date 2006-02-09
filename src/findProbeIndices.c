#include <R.h>

void findIndices(int *probeID, int *probeList, int *nrow, int *indices){
  int i, j;
  j = 0;

  for(i=0; i<*nrow; i++){
    if(probeList[i] == *probeID){
      indices[j] = i+1;
      j++;
    }
  }
}
