/*------------------------------------------------------------------------
 *   round sources and receivers to gridpoint        
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

int irnd_to_grid_max(float xpos, float *xg, int nx1, int nx2, int dx){

	/* local variables */
	int	ipos;
	int	i;


	/* find index of closest grid point */
	ipos = nx2;
	if (xpos < xg[nx2]){
		for (i=nx1;i<=nx2;i+=dx){
			if (xpos < xg[i]){
				ipos = i-1;
				break;
			}
		}
	}


	return ipos;
}
