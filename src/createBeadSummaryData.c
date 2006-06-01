#include "beadarray.h"
#include <stdlib.h>
#include <math.h>
#include <R.h>

double sd(double *arr, int length, double mean){

    int i;
	double sum, temp;

	sum = 0;

	for(i = 0; i < length; i++){
	    temp = arr[i] - mean;
		sum += pow(temp, 2);
		}

	temp = sum/(length-1);
	return(sqrt(temp));
}

void quicksortDouble(double arr[], int low, int high) {
 int i = low;
 int j = high;
 double y = 0;
 /* compare value */
 double z = arr[(low + high) / 2];

 /* partition */
 do {
  /* find member above ... */
  while(arr[i] < z) i++;

  /* find element below ... */
  while(arr[j] > z) j--;

  if(i <= j) {
   /* swap two elements */
   y = arr[i];
   arr[i] = arr[j]; 
   arr[j] = y;
   i++; 
   j--;
  }
 } while(i <= j);

 /* recurse */
 if(low < j) 
  quicksortDouble(arr, low, j);

 if(i < high) 
  quicksortDouble(arr, i, high); 
}

double median(double *arr, int length){
	int i;
	double j;
	double *tempArr;

	tempArr = (double *)malloc(sizeof(double) * length);

	for(i = 0; i < length; i++){
		  tempArr[i] = arr[i];
	}

	quicksortDouble(tempArr, 0, (length-1));
	i = floor(length/2);

	if(i == (length/2.0)){ //if its of even length take the average of the middle two
		j = (tempArr[i-1]+tempArr[i])/2;
	}
	else { //if odd length
		j = tempArr[i];
	}

	free(tempArr);
	return j;
}

double mad(double *arr, int length){
	int i;
	double temp, med;
	double *arr2;

	med = median(arr, length);

	arr2 = (double *)malloc(sizeof(double) * length);

	for(i = 0;i < length; i++){
		  temp = (arr[i] - med);
		  arr2[i] = fabs(temp);
	}
	temp = 1.4826*(median(arr2, length));
	free(arr2);
	return(temp);
}

int* getProbeIndices(int *probeList, int probeID, int *start, int numBeads){
	 int i;
	 int *indices;

	 indices = (int *)malloc(sizeof(int) * 2);
	 i = *start;

/* Dodgy as hell.  If the probeID isn't present it'll just access memory outside of the probeList 
   array however it still seems to work and returns NaN.  Really should add decent error checking but 
   I'll do it later */
	 while(probeList[i] < probeID){
	 		i++;
		}
	
	 *start = i;

	 while(probeList[i] == probeID){
	 		i++;
		}

		indices[0] = *start;
		indices[1] = i-1;

	return(indices);
}

beadStatusStruct* findBeadStatus(double *intensities, int *probeList, int probeID, int numBeads, int *count, int *start){
	beadStatusStruct *status;
	int *outlierInds, *validInds, *indices;
	double *inten;
	double m, ma;
	int i,j,k, nsize;

	//get the start and end indices for this probeID
	indices = getProbeIndices(probeList, probeID, start, numBeads);

	*start = (indices[1]+1);  //set the next start position
	*count = (indices[1] - indices[0] + 1);  //calculate how many beads there are of this type
	inten = (double *)malloc(sizeof(double) * (*count));

	i = indices[0];
	j = indices[1];
	k = 0;
	
	while(i < (j+1)){
		inten[k] = intensities[i];
		i++;
		k++;
	}

	m = median(inten, (*count));
	ma = mad(inten, (*count));

	status = (beadStatusStruct *)malloc(sizeof(beadStatusStruct));
	validInds = (int *)malloc(sizeof(int));
	outlierInds = (int *)malloc(sizeof(int));

	i = j = k = 0;	

	/* find the indices of valid probes */
	while(k < (*count)){
		if((inten[k] < (m + 3*ma)) && (inten[k] > (m - 3*ma))){
			validInds[i] = (indices[0]+k);
			i++;
			nsize = (sizeof(int) * (i+1));
			validInds = realloc(validInds, nsize);	
		}
		else{
			outlierInds[j] = (indices[0]+k);
			j++;
			nsize = (sizeof(int) * (j+1));
			outlierInds = realloc(outlierInds, nsize);
		}
		k++;
	}

	validInds[i] = -1;
	outlierInds[j] = -1;

	status->validInds = validInds;
	status->outlierInds = outlierInds;

	free(inten);
	free(indices);

	return(status);
}

void findAllOutliers(double *finten, int *binaryStatus, int *probeList, int *probeIDs, int *numProbes, int *numBeads, int *nextStart){

	beadStatusStruct *status; 
	int *valids, *count;
	int k, i, probeID, temp;

	count = (int *)malloc(sizeof(int));
	*count = 0; //initialise counter to record how many beads there are of this type

	for(k = 0; k < (*numProbes); k++){
		  
		  probeID = probeIDs[k];
		  status = findBeadStatus(finten, probeList, probeID, *numBeads, count, nextStart);
		  valids = status->validInds;
//		  outliers = status->outlierInds;

		  i = 0;
		  while(valids[i] != -1){
		  		temp = valids[i];
		  		binaryStatus[temp] = 1;
				i++;
				}
		}
}

void createBeadSummary(double *finten, double *binten, int *probeList, int *probeIDs, int *numProbes, int *numBeads, double *foreground, 
	 						  double *background, double *stdev, int *numValid, int *numOutlier, int *nextStart){
	 
	 beadStatusStruct *status;
	 int *valids, *count, i, k, ind, probeID;
	 double *validInten;
	 double fground, bground;

	 count = (int *)malloc(sizeof(int));
	 *count = 0; //initialise counter to record how many beads there are of this type

	 for(k = 0; k < (*numProbes); k++){

	 probeID = probeIDs[k];

	 status = findBeadStatus(finten, probeList, probeID, *numBeads, count, nextStart);
	 valids = status->validInds;

	 i=0;
	 fground = bground = 0;

	 while(valids[i] != -1){
	 	ind = valids[i];
	 	fground += finten[ind];
		bground += binten[ind];
		i++;
	 }
	 foreground[k] = (fground/i);
	 background[k] = (bground/i);

	 validInten = (double *)malloc(sizeof(double) * i);

	 i = 0;
	 while(valids[i] != -1){
	 	ind = valids[i];
	 	validInten[i] = finten[ind];
		i++;
	}

	stdev[k] = sd(validInten, i, foreground[k]);
	numValid[k] = i;
	numOutlier[k] = ((*count)-i);

	free(validInten);
	free(valids);
	free(status->validInds);
	free(status->outlierInds);
	free(status);
	}
	
	free(count);
}



