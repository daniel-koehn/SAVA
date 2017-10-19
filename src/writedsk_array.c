/*------------------------------------------------------------------------
 *   Write one single amplitude on disk                                   
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


void writedsk_array(FILE *fp_out, int na, float *amp, int format){

	int	i;
#ifdef _CRAY1     
	float	*ampc = NULL;
	int	type=2, num=na, bitoff=0, stride=1, ierr;
#endif


	switch(format){
	case 1 :	/* SU */
		error(" Sorry, SU-format for snapshots not implemented yet. \n");
		break;
	case 2 :	/* ASCII */
		for (i=1;i<=na;i++)
			fprintf(fp_out,"%e\n",amp[i]);
		break;
	case 3 :	/* BINARY */
#ifdef _CRAY1     
		ampc = vector(1,na);
		ierr = CRAY2IEG(&type, &num, ampc, &bitoff, amp, &stride);
		fwrite(ampc,4,na,fp_out);
		free_vector(ampc,1,na);
#else 
		fwrite(amp,sizeof(float),na,fp_out);
#endif
		break;
	default :
		printf(" Don't know the format for the snapshot-data !\n");
		error(" No output was written. ");
	}
}
