#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "beadarray.h"

void kth(double arr[], int low, int high, int k)
{
	/*sorts kth elt into correct position, k counting from 0*/
	int i = low;
	int j = high;
	double y = 0.0;
	/* compare value */
	double z = arr[(low + high) / 2];

	/* partition */
	do
	{
		/* find member above ... */
		while(arr[i] < z) i++;

		/* find element below ... */
		while(arr[j] > z) j--;

		if(i <= j)
		{
			/* swap two elements */
			y = arr[i];
			arr[i] = arr[j]; 
			arr[j] = y;
			i++; 
			j--;
		}
	} while(i <= j);

 if(low >= high){return;}

 //continue in the relevant partition (only)
 if(k >= i)
 {kth(arr, i, high, k);}
 else if(k <= j - low)
 {kth(arr, low, j, k);}

return;
}

double mediansort(double* list, int high)
{
	double median = 0.0;
	if(high % 2 == 0)
	{
		kth(list, 0, high, high/2);
		median = list[high/2];
	}
	else
	{
		kth(list, 0, high, (high-1)/2);
		//want k + 1 also
		double lowest = 0.0;
		for (int i = (high+1)/2; i < high; i++)
		{
			if(list[i] < lowest) {lowest = list[i];}
		}

		median = (list[(high-1)/2] + lowest)/2.0;
	}
	return median;
}

void Neighbours(double* x, double* y, int* no, int* neighbours, double* thresh, double* margin, double* window, int* wcols, int* wrows)
{
	int i,j,k,l;
	int wrow, wcol;
	int best[7], tempint, tempinta, grab;
	int nowindow, nomargin;
	double mwindow = *window + *margin;
	double bestdist[7];
	double curr, temp;
	double left, right, top, bottom;
	double xdiff, ydiff;

        int *margined = (int *) R_alloc(*no, sizeof(int));
	int *windowed = (int *) R_alloc(*no, sizeof(int));

	for(wrow = 1; wrow <= *wrows; wrow++)
	{
	for(wcol = 1; wcol <= *wcols; wcol++)
	{
		//consider a smaller group of nodes
		nowindow = 0;
		nomargin = 0;
		top = *margin + ((double)wrow)*2*(*window);
		bottom = top - 2.0*mwindow;
		right = *margin + ((double)wcol)*2*(*window);
		left = right - 2.0*mwindow;

		//printf("(%lf, %lf, %lf, %lf)", left, right, bottom, top);

		for(i = 0; i < *no; i++)
		{
			if(x[i] > left && x[i] < right && y[i] > bottom && y[i] < top)
			{
				margined[nomargin] = i;
				nomargin++;

				if(x[i] > left + *margin && x[i] < right - *margin && y[i] > bottom + *margin && y[i] < top - *margin)
				{
					windowed[nowindow] = i;
					nowindow++;
				}
			}
		}

		//for each windowed node...
		for(i = 0; i < nowindow; i++)
		{
			//initialise arrays
			for(j = 0; j < 7; j++)
			{
				best[j] = 0;
				bestdist[j] = 99999.9;
			}

			//check prospective neighbours
			for(j = 0; j < nomargin; j++)
			{
				xdiff = x[windowed[i]]-x[margined[j]];
				ydiff = y[windowed[i]]-y[margined[j]];
				curr = xdiff*xdiff + ydiff*ydiff;

				//rearrange to correct slot (to keep ascending order)
				if (curr < bestdist[6])
				{
					best[6] = margined[j] + 1; //translate to R index now
					bestdist[6] = curr;

					k = 6;
					while(k > 0 && bestdist[k-1] > bestdist[k])
					{
						temp = bestdist[k];
						bestdist[k] = bestdist[k-1];
						bestdist[k-1] = temp;

						tempint = best[k];
						best[k] = best[k-1];
						best[k-1] = tempint;
						k--;
					}
				}
			}

			//remove those too far away
			grab = 6;
			for(int j = 6; j > 3; j--)
			{
				if(bestdist[j] > (*thresh)*bestdist[j-1]) {grab = j-1;}
			}

			//write to neighbours
			for(int j = 1; j <= grab; j++)
			{
				tempint = 6*windowed[i] + (j-1);
				neighbours[tempint] = best[j];
			}
		}
	}
	}
	
	//now delete non-mutual links
	for(i = 0; i < *no; i++)
	{
		for(j = 0; j < 6; j++)
		{
			grab = 6*i + j;
			tempinta = neighbours[grab] - 1;
			//check node (tempinta), does it connect back to node i? (remembering to subtract 1 for R index -> C index)

			if(tempinta >= 0)
			{
				l=0;
				tempint = 6*tempinta;

				for(k = 0; k < 6; k++)
				{
					l += (neighbours[tempint] == i + 1);
					tempint++;
				}
				if (l == 0)
				{
					neighbours[grab] = 0;
				}
			}
		}
	}
}

void BGFilter(double* E, double* Etilde, int* neighbours, int* nbeads, int* invasions, int* method)
{
	int i, j, k, l, temp, temp2;
	int nextslot, firstinv, lastinv;
	double median, mad;

	int *invaded = (int *) R_alloc(10 * (*invasions) * (*invasions + 1), sizeof(int) ); //must be big enough to contain IDs of all invaded nodes 
	double *invadedE = (double *) R_alloc(10 * (*invasions) * (*invasions + 1), sizeof(double));
	int *done = (int *) R_alloc(*nbeads, sizeof(int));
	//initialise done matrix to 0
        memset(done, 0, sizeof(int) * *nbeads);

	//for each node...
	for(i = 0; i < *nbeads; i++) // was i <= nbeads - 1
	{
		nextslot = 1;
		done[i] = 1;
		invaded[0] = i;
		invadedE[0] = E[i];
		firstinv = 0;
		lastinv = 0;

		//invasion process for this node
		//----------------
		//invade (invasions) times please
		for(l = 1; l <= *invasions; l++)
		{
			//for each invaded node...
			for(j = firstinv; j <= lastinv; j++)
			{
				temp = 6*invaded[j]; //first candidate neighbour (cannot pass arrays, hence using a vector form)

				//check its neighbours
				for(k = 0; k <=5; k++)
				{
					//did we invade? if not, let's do it (remembering to avoid 0 slots)
                                        if(neighbours[temp] != 0)
                                        {
                                                temp2 = neighbours[temp] - 1;
                                                if(!done[temp2])
                                                {        
							invaded[nextslot] = temp2;
							invadedE[nextslot] = E[temp2];
							done[temp2] = 1;
							nextslot++;
						}
					}
					temp++;
				}
			}
			firstinv = lastinv + 1;
			lastinv = nextslot - 1;
		}

		//choose analysis
		switch(*method)
		{
			case 1:
			//median filter
			Etilde[i] = E[i] - mediansort(invadedE, nextslot - 1);
			break;
			
			case 2:
			//mean filter
			Etilde[i] = 0;
			for(j = 0; j < nextslot; j++)
			{
				Etilde[i] += invadedE[j];
			}
			Etilde[i] = E[i] - Etilde[i]/nextslot;
			break;

			case 3:
			//"MAD filter" (i.e. scale by MAD of nearby beads)
			case 4:
			//"Scale MAD filter" (subtract median and scale automatically)

			//find median first
			median = mediansort(invadedE, nextslot - 1);

			//now find ADs
			for(j = 0; j < nextslot; j++)
			{
				invadedE[j] = fabs(invadedE[j] - median);
			}
			//sort to find MAD
			mad = mediansort(invadedE, nextslot - 1);

			//output (3: mad) or (4: scaled E)
			if(*method == 3) {Etilde[i] = E[i]/mad;}
			else {Etilde[i] = (E[i] - median)/mad;}							
			break;
		}

		//clear done flags
		for(j = 0; j < nextslot ; j++)
		{
			done[invaded[j]] = 0;
		}
	}
	return;
}

void BGFilterWeighted(double* E, double* Etilde, int* neighbours, double* weights, int* nbeads, int* invasions)
{
//this method downweights according to given weights, and also according to distance
	int i, j, k, l, temp, temp2;
	int nextslot, firstinv, lastinv;
	double wtotal = 0;

	int *invaded = (int *) R_alloc(10 * (*invasions) * (*invasions + 1), sizeof(int) ); //must be big enough to contain IDs of all invaded nodes 
	double *invadedE = (double *) R_alloc(10 * (*invasions) * (*invasions + 1),  sizeof(double) );
	double *invadedW = (double *) R_alloc(10 * (*invasions) * (*invasions + 1), sizeof(double) );

	int *done = (int *) R_alloc(*nbeads, sizeof(int));
	//initialise done matrix to 0
        memset(done, 0, *nbeads * sizeof(int));

	//for each node...
	for(i = 0; i < *nbeads; i++) // was i <= nbeads - 1
	{
		nextslot = 1;
		done[i] = 1;
		invaded[0] = i;
		invadedE[0] = E[i];
		invadedW[0] = 1;
		firstinv = 0;
		lastinv = 0;

		//invasion process for this node
		//----------------
		//invade (invasions) times please
		for(l = 1; l <= *invasions; l++)
		{
			//for each invaded node...
			for(j = firstinv; j <= lastinv; j++)
			{
				temp = 6*invaded[j]; //first candidate neighbour (cannot pass arrays, hence using a vector form)

				//check its neighbours
				for(k = 0; k <=5; k++)
				{
					//did we invade? if not, let's do it (remembering to avoid 0 slots)
                                        if(neighbours[temp] != 0)
                                        {
                                                temp2 = neighbours[temp] - 1;
                                                if(!done[temp2])
                                                {        
                                                        invaded[nextslot] = temp2;
                                                        invadedE[nextslot] = E[temp2];
                                                        done[temp2] = 1;
                                                        nextslot++;
                                                }
                                        }
					temp++;
				}
			}
			firstinv = lastinv + 1;
			lastinv = nextslot - 1;
		}

		//weighted mean
		wtotal = 0;
		Etilde[i] = 0;

		for(j = 0; j < nextslot; j++)
		{
			Etilde[i] += invadedE[j]*invadedW[j];
			wtotal += invadedW[j];
		}
		Etilde[i] = Etilde[i]/(wtotal);

		//clear done flags
		for(j = 0; j < nextslot ; j++)
		{
			done[invaded[j]] = 0;
		}
	}
	return;
}

void Flood(int ID, int* neighbours, int CID, int* clusterID, int* clustersize)
{
	int i, temp;

	//convert this node
	clusterID[ID] = CID;
	clustersize[CID] += 1;

	for(i=0; i<=5; i++)
	{
		//index in the neighbours vector
		temp = (6*ID + i);

		//if it's an outlier...
		if(neighbours[temp] > 0)
		{
			//if it hasn't been flooded...
			if(clusterID[neighbours[temp] - 1] == 0)
			{
				//flood!
				Flood(neighbours[temp] - 1, neighbours, CID, clusterID, clustersize);
			}
		}
	}
}

void FloodFill(int* neighbours, int* outlierIDs, int* noutliers, int* clusterID, int* clustersize)
{
	int CID = 1;
	int i, ID;

	//This should return a vector of clusterIDs. (Links to non-outliers, designated as 0 in the mx, are ignored.)
	for(i = 0; i < *noutliers; i++)
	{
		//convert to C indexing
		ID = outlierIDs[i] - 1;

		if(clusterID[ID] != 0)
		{
			//printf("%d:%d ", ID, clusterID[ID]);
		}

		//did we flood this yet? if not, begin flooding from here.
		if(clusterID[ID] == 0)
		{
			//printf("%d ", ID);

			Flood(ID, neighbours, CID, clusterID, clustersize);
			CID++;
		}
	}
}

void DiffuseDefects(int* neighbours, int* IDs, int* nbeads, int* nIDs, int* ncompacts, int* invasions, double* output, double* sig)
{
	//We want to find which beads have a kernel with too many outliers (IDs) in them.
	int i, j, k, l, temp, temp2;
	int nextslot, firstinv, lastinv;
	int invaded[50000];
	int flagged;
	double test = 0.0;
	double p = ((double)*nIDs)/(double)(*nbeads - *ncompacts);

	//initialise done matrix to 0
	int *done = (int *) R_alloc(*nbeads, sizeof(int));
        memset(done, 0, *nbeads * sizeof(int));

	//calculate threshold vals for binomial distribution
	double *choose = (double *) R_alloc(10 * (*invasions) * (*invasions + 1), sizeof(double) );
	int *thresh = (int *) R_alloc(10 * (*invasions) * (*invasions + 1), sizeof(int) );

	choose[0] = 1.0;
	choose[1] = 0.0;
	for(i = 1; i <= 5*(*invasions)*(*invasions + 1); i++)
	{
		//update the choose matrix to be (i)C(j)
		for(j = i-1; j > 0; j--) {choose[j] = choose[j] + choose[j-1];}
		choose[i] = 1.0;

		//under H0, X ~ Bin(i,p)
		//bring x down from the top until P(X > x) > *sig
		test = 0.0;
		j = i;

		while(test < *sig)
		{
			test += choose[j]*pow(p, j)*pow(1.0-p, i-j);
			j--;
		}
		thresh[i] = j;
	}

	//initialise outlier matrix
	int *outlier = (int *) R_alloc(*nbeads, sizeof(int));
	for(i = 0; i < *nbeads; i++)
	{
		outlier[i] = 0;
	}
	for(i = 0; i < *nIDs; i++)
	{
		outlier[IDs[i]-1] = 1;
	}

	//for each node...
	for(i = 0; i <= *nbeads - 1; i++)
	{
		nextslot = 1;
		done[i] = 1;
		invaded[0] = i;
		flagged = outlier[i];
		firstinv = 0;
		lastinv = 0;

		//invasion process for this node
		//----------------
		//invade (invasions) times please
		for(l = 1; l <= *invasions; l++)
		{
			//for each invaded node...
			for(j = firstinv; j <= lastinv; j++)
			{
				temp = 6*invaded[j]; //candidate neighbour (cannot pass arrays, hence using a vector form)
				//check its neighbours
				for(k = 0; k <=5; k++)
				{
					//did we invade? if not, let's do it (remembering to avoid 0 slots)
					if(neighbours[temp] != 0)
					{
						temp2 = neighbours[temp] - 1;
						if(!done[temp2])
						{
							invaded[nextslot] = temp2;
							flagged += outlier[temp2]; //count on if it's an outlier
							done[temp2] = 1;
							nextslot++;
						}
					}
					temp++;
				}
			}
			firstinv = lastinv + 1;
			lastinv = nextslot - 1;
		}

		//check against threshold value
		output[i] = flagged > thresh[nextslot];

		//clear done flags
		for(j = 0; j < nextslot ; j++) {done[invaded[j]] = 0;}
	}
	return;
}

/*void Dilate(int* IDs, int* nIDs, int* neighbours, int* nbeads, int* invasions)*/
/*{*/
/*	int i, j, k, l, temp;*/
/*	int nextslot, firstinv, lastinv;*/
/*	int *done = malloc(*nbeads * sizeof(int));*/

/*	//initialise done matrix*/
/*	for(i = 0; i < *nbeads; i++)*/
/*	{*/
/*		done[i] = 0;*/
/*	}*/
/*	for(i = 0; i < *nIDs; i++)*/
/*	{*/
/*		done[IDs[i]-1] = 1;*/
/*	}*/

/*	nextslot = *nIDs;*/
/*	firstinv = 0;*/
/*	lastinv = *nIDs - 1;*/

/*	//invasion process for this node*/
/*	//----------------*/
/*	//invade (invasions) times please*/
/*	for(l = 1; l <= *invasions; l++)*/
/*	{*/
/*		//for each invaded node...*/
/*		for(j = firstinv; j <= lastinv; j++)*/
/*		{*/
/*			temp = 6*(IDs[j]-1);  //first candidate neighbour (cannot pass arrays, hence using a vector form)*/

/*			//check its neighbours*/
/*			for(k = 0; k <=5; k++)*/
/*			{*/
/*				//did we invade? if not, let's do it (remembering to avoid 0 slots)*/
/*				if(neighbours[temp] != 0 && !done[neighbours[temp]-1])*/
/*				{*/
/*					IDs[nextslot] = neighbours[temp];*/
/*					done[neighbours[temp]-1] = 1;*/
/*					nextslot++;*/
/*				}*/
/*				temp++;*/
/*			}*/
/*		}*/
/*		firstinv = lastinv + 1;*/
/*		lastinv = nextslot - 1;*/
/*	}*/
/*}*/

void Close(int* IDs, int* nIDs, int* neighbours, int* nbeads, int* invasions)
{
	int i, j, k, l, temp;
	int nextslot, firstinv, lastinv;
	int *done = (int *) R_alloc(*nbeads, sizeof(int));
	int *marker = (int *) R_alloc((*invasions + 1), sizeof(int));

	marker[0] = 0;
	//initialise done matrix
        memset(done, 0, *nbeads * sizeof(int));
	for(i = 0; i < *nIDs; i++)
	{
		done[IDs[i]-1] = 1;
	}

	nextslot = *nIDs;
	firstinv = 0;
	lastinv = *nIDs - 1;

	//invasion process for this node
	//----------------
	//invade (invasions) times please
	for(l = 1; l <= *invasions; l++)
	{
		//for each invaded node...
		for(j = firstinv; j <= lastinv; j++)
		{
			temp = 6*(IDs[j]-1);  //first candidate neighbour (cannot pass arrays, hence using a vector form)

			//check its neighbours
			for(k = 0; k <=5; k++)
			{
				//did we invade? if not, let's do it (remembering to avoid 0 slots)
				if(neighbours[temp] != 0 && !done[neighbours[temp]-1])
				{
					IDs[nextslot] = neighbours[temp];
					done[neighbours[temp]-1] = 1;
					nextslot++;
				}
				temp++;
			}
		}
		firstinv = lastinv + 1;
		marker[l] = firstinv; //aids retreat a little bit
		lastinv = nextslot - 1;
	}

	//retreating process for this node
	//----------------
	//retreat (invasions) times please
	for(l = 1; l <= *invasions; l++)
	{
		//for each retreatable node...
		for(j = marker[*invasions - l]; j <= lastinv; j++)
		{
			if(IDs[j] > 0) //only those not yet retreated...
			{
				temp = 6*(IDs[j]-1);  //first candidate neighbour (cannot pass arrays, hence using a vector form)

				//check its neighbours
				for(k = 0; k <=5; k++)
				{
					//should we retreat? (remembering to avoid 0 slots)
					if(neighbours[temp] != 0 && !done[neighbours[temp]-1])
					{
						done[IDs[j] - 1] = 2;
						k = 10;
					}
					temp++;
				}
			}
		}

		//insert code to convert 2s to 0s, and blank dead ID slots
		for(j = marker[*invasions - l]; j <= lastinv; j++)
		{
			if(IDs[j] > 0)
			{
				if(done[IDs[j] - 1] == 2)
				{
					done[IDs[j] - 1] = 0;
					IDs[j] = 0;
				}
			}
		}
	}
}
