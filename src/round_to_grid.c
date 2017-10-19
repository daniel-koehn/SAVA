/*------------------------------------------------------------------------
 *   round sources and receivers to gridpoint        
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

int *irnd_to_grid(float xpos, float *xg, float *xpg, int nx){

	/* local variables */
	float	d, dp, dmin, dpmin;
	int	*ipos=NULL;
	int	i;

	
	ipos = ivector(1,2);
	ipos[1] = 1;
	ipos[2] = 1;
	
	/* find index of closest grid point */
	dmin  = xpos- xg[1];
	if (dmin < 0.0)
		dmin = -dmin;
	dpmin = xpos-xpg[1];
	if (dpmin < 0.0)
		dpmin = -dpmin;

	for (i=2;i<=nx;i++){
		d  = xpos-xg[i];
		if (d < 0.0)
			d = -d;
		dp = xpos-xpg[i];
		if (dp < 0.0)
			dp = -dp;
		
		if (d < dmin){
			dmin = d;
			ipos[1] = i;
		}
		if (dp < dpmin){
			dpmin = dp;
			ipos[2] = i;
		}
		else{
			break;
		}
	}

	return ipos;
}
