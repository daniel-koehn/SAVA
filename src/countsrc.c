/*------------------------------------------------------------------------ 
 *   Count number of source positions
 *
 *   D. Koehn
 *   28.10.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void countsrc(){

	/* extern variables */
	extern char	SRC_FILE[STRING_SIZE];
	extern int	NSRC;
	extern FILE	*FP;

	/* local variables */
	FILE	*fpsrc;

	fpsrc = fopen(SRC_FILE,"r");
	
	if (fpsrc==NULL) error(" Source file could not be opened !");
	
	NSRC=0;

	/* read number of source positions */
	fscanf(fpsrc,"%d",&NSRC);
	fprintf(FP," Number of source positions found: %d \n",NSRC);
        fclose(fpsrc);	

}
