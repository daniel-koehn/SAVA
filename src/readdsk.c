/*------------------------------------------------------------------------
 *   Read one single amplitude from file                                   
 *
 *   T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/


float readdsk(FILE *fp_in, int format){

	/* local variables */
	float	amp;


	switch(format){
		case 1 :	/* SU*/ 
			error(" Sorry, SU-format for snapshots not implemented yet. \n");
			break;
		case 2 :	/*ASCII*/
			fscanf(fp_in,"%e\n", &amp); 
			break;
		case 3 :	/* BINARY */
			fread(&amp, sizeof(float), 1, fp_in);
			break;
		default :
			printf(" Don't know the format for the snapshot-data !\n");
			error(" No output was written. ");
	}

	return amp;
}
