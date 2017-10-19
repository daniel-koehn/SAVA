/*------------------------------------------------------------------------
 *   Read one single amplitude from file                                   
 *
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/


void readdsk_array(FILE *fp_in, int na, float *amp, int format){

	int	i;
	float	ampc;

	
	switch(format){
	case 1 :	/* SU*/ 
		error(" Sorry, SU-format for snapshots not implemented yet. \n");
		break;
	case 2 :	/*ASCII*/
		for (i=1;i<=na;i++){
			fscanf(fp_in,"%e\n",&ampc);
			amp[i] = ampc;
		}
		break;
	case 3 :	/* BINARY */
		fread(amp,sizeof(float),na,fp_in);
		break;  
	default :
		printf(" Don't know the format for the snapshot-data !\n");
		error(" No output was written. ");
	}

}
