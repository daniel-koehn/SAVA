/*------------------------------------------------------------------------
 *   generate point fource source at source nodes
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void fsource(int nt, struct vector3d ***v, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk,
float ** srcpos_loc, float ** signals, int nsrc, 
float * dx, float * dy, float * dz, float * dxp, float * dyp, float * dzp){

	/* extern variables */
	extern float	DT;

	/* local variables */
	int	i, j, k, l, type;
	float	amp;

	/* adding source wavelet to the particle velocity components 
	   (force source) at the source points */


	for (l=1;l<=nsrc;l++) {
		type = (int)srcpos_loc[l][7];
		if ((type>1) && (type<8)){
			i   = (int)srcpos_loc[l][1];
			j   = (int)srcpos_loc[l][2];
			k   = (int)srcpos_loc[l][3];
			amp = signals[l][nt];
			
			switch (type){
			case 2 :
				v[i][j][k].x += DT*amp*dxp[i]*dy[j]*dz[k]; /* acceleration in x-direction */
				break;
			case 3 :
				v[i][j][k].y += DT*amp*dx[i]*dyp[j]*dz[k]; /* acceleration in y-direction */
				break;
			case 4 :
				v[i][j][k].z += DT*amp*dx[i]*dy[j]*dzp[k]; /* acceleration in z-direction */
				break;
			case 5 :
				v[i][j][k].x += amp*dxp[i]*dy[j]*dz[k]*rhoijpkp[i][j][k]; /* single force in x-direction */
				break;
			case 6 :
				v[i][j][k].y += amp*dx[i]*dyp[j]*dz[k]*rhoipjkp[i][j][k]; /* single force in y-direction */
				break;
			case 7 :
				v[i][j][k].z += amp*dx[i]*dy[j]*dzp[k]*rhoipjpk[i][j][k]; /* single force in z-direction */
				break;
			}
		}
	}
}
