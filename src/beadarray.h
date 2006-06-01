//Structure definition
typedef struct
{
 		int *validInds;
		int *outlierInds;
}
beadStatusStruct;

//Function prototypes
void quicksort(int*, int, int); 
int number(int, char*);
void sharpen(int**, int, int);
void asf(int**, int, int);
double matrixMean(int**, int, int);
void calculateBackground(int**, double*, double*, int, int, int, double*, int);
void HIPForeground(int**, double*, double*, int, int, int, double*);
void IlluminaForeground(int**, double*, double*, int, int, int, double*);
