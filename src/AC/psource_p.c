/*------------------------------------------------------------------------
 *   generate P-wave source or douple couple source at source nodes
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void psource(int nt, float *** p, 
float ** srcpos_loc, float ** signals, int nsrc_loc, 
float * dx, float * dy, float * dz){

	/* extern variables */
	extern float	DT;

	/* local variables */
	int	i, j, k, l, type;


	/* adding source wavelet to stress components at source points */
	for (l=1;l<=nsrc_loc;l++) {
		type = (int)srcpos_loc[l][7];
		if (type==1){
			i   = (int)srcpos_loc[l][1];
			j   = (int)srcpos_loc[l][2];
			k   = (int)srcpos_loc[l][3];

			/* pressure */
			p[i][j][k] += DT*signals[l][nt]*dx[i]*dy[j]*dz[k];
		}
	}
}
