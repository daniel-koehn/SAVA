/*------------------------------------------------------------------------
 *   Copy vector with amplitudes from one file to another file                                  
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/


void copydsk_array(FILE *fp_out, FILE *fp_in, int na, int format){

	int	i;
	float	amp;
	float	*array=NULL; 
#ifdef _CRAY1
	float	*arrayc=NULL;
	int	type=2, num=na, bitoff=0, stride=1, ierr;
#endif


	switch(format){
	case 1 :	/* SU */
		error(" Sorry, SU-format for snapshots not implemented yet. \n");
		break;
	case 2 :	/* ASCII */
		for (i=1;i<=na;i++)
			fscanf(fp_out,"%e\n",&amp);
			fprintf(fp_in,"%e\n",amp);
		break;
	case 3 :	/* BINARY */
		array = vector(1,na);
	  
		fread(array,sizeof(float),na,fp_out);
#ifdef _CRAY1   
		arrayc = vector(1,na);

		ierr = CRAY2IEG(&type, &num, arrayc, &bitoff, array, &stride);
		fwrite(arrayc,4,na,fp_out);

		free_vector(arrayc,1,na);
#else
		fwrite(array,sizeof(float),na,fp_in);
#endif
		free_vector(array,1,na);
		break;
	default :
		printf(" Don't know the format for the snapshot-data !\n");
		error(" No output was written. ");
	}
}
