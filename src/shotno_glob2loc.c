/*------------------------------------------------------------------------
 *   Estimate local shot number index for given shot
 *
 *   D. Koehn
 *   Kiel, 22.10.2017
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "ISO_struct.h"	/* data structures for isotropic elastic forward modelling */

void shotno_glob2loc(struct acq *acq){

        /* global variables */
	extern int RUN_MULTIPLE_SHOTS, NSRC_LOC, SHOTNO, SHOTNO_LOC;
	extern int FFID;

        /* local variables */
	int l;

	SHOTNO_LOC = 0;
		
        /* find local source position of SHOTNO */
	if(RUN_MULTIPLE_SHOTS && NSRC_LOC){	
	    for(l=1;l<=NSRC_LOC;l++){	    
	        if((*acq).srcpos_loc[l][9]==SHOTNO){
		    SHOTNO_LOC = l;  
		    break;
		}	    
	    }
	}

	FFID = SHOTNO;

}
