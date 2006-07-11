#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "beadarray.h"

/* A function for opening a tiff file and storing the pixel intensities in 
a matrix in R.
Plenty of code replication between here and the start of readBeadImages.
I'll tidy it up at some point. Now I'm just pleased it seems to work */

SEXP readTIFF(SEXP fileName){

	  SEXP ans;
 	  FILE *fp;
	  char *tif[1];
      int i, j, BeginData, len, ImageWidth, ImageHeight, ImageSize, StripOffset;
      long tag, type, length, offset, records;
      char c2[2], c4[4];    
      int **pixels;

	  PROTECT(fileName = AS_CHARACTER(fileName)); 

	  tif[0] = R_alloc(strlen(CHAR(STRING_ELT(fileName, 0))), sizeof(char)); 
	  strcpy(tif[0], CHAR(STRING_ELT(fileName, 0))); 
	  UNPROTECT(1);

      if ( (fp = fopen(tif[0], "rb")) == NULL)
      {
         Rprintf("Error opening file %s", tif[0]);
         exit(0);
      }

      if (fread(c2, 1, 2, fp) == 2) {
//         Rprintf("Mode: %u %u\n", c2[0], c2[1]);
         }
         
      if (fread(c2, 1, 2, fp) == 2) {
//         Rprintf("Version: %i\n", number(2, c2));
         }
         
      if (fread(c4, 1, 4, fp) == 4) {
         offset = number(4, c4);
//         Rprintf("Offset: %i\n", offset);
         }
         
      if (fseek(fp, offset, SEEK_SET) != 0) {
//         Rprintf("Error in fseek()\n");
         exit(0);
         }
         
   if (fread(c2, 1, 2, fp) == 2) {
      records = number(2, c2);
//      Rprintf("Records: %i\n", records);
   }
   for (i = 1; i <= records; i++) {
      fread(c2, 1, 2, fp);
      tag = number(2, c2);
      fread(c2, 1, 2, fp);
      type = number(2, c2);
      fread(c4, 1, 4, fp);
      length = number(4, c4);
      fread(c4, 1, 4, fp);
      offset = number(4, c4);
  //    Rprintf("%i %i %i %i\n", tag, type, length, offset);
      switch (tag) {
         case (256) :
            ImageWidth = offset;
         case (257) :
            ImageHeight = offset;
         case (273) : 
            StripOffset = offset;
      }
   }
   if (fseek(fp, StripOffset, SEEK_SET) != 0) {
      Rprintf("Error in fseek()\n");
      exit(0);
   }
   if (fread(c4, 1, 4, fp) == 4) {
      BeginData = number(4, c4);
   }
   ImageSize = ImageWidth*ImageHeight*2;
   if (fseek(fp, 0, SEEK_END) != 0) {
      Rprintf("Error in fseek()\n");
      exit(0);
   }
   len = ftell(fp);
   if ((BeginData) + (ImageSize) > len) {
      BeginData = StripOffset; /* it's not a tif */
   }

   rewind(fp);
   
   pixels = malloc(sizeof(int *) * ImageWidth);

      for(i = 0;i < ImageWidth; i++){
           pixels[i] = malloc(sizeof(int) * ImageHeight);
           }
           
//   Rprintf("Reading pixels of %s\n", tif[0]); 
   getPixelIntensities(pixels, fp, ImageWidth, ImageHeight, ImageSize, BeginData);

   PROTECT(ans = allocMatrix(REALSXP, ImageWidth, ImageHeight));
   for(i = 0; i < ImageWidth; i++){
   		 for(j = 0; j < ImageHeight; j++){
		   		 REAL(ans)[i + ImageWidth*j] = pixels[i][j];
		 }
	}
	
	for(i = 0; i < ImageWidth; i++){
         free(pixels[i]);
   	}
   		 free(pixels);
	UNPROTECT(1);
	return(ans);
}
