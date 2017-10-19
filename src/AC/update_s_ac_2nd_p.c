/*------------------------------------------------------------------------
 *   updating stress tensor components
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s_ac(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
	struct vector3d ***v, float ***p, float ***div,
	float *** pi, 
	float *** absorb_coeff,
	float * dx, float * dy, float * dz, int infoout){

	/* local variables */
	int	i, j, k;
	float	divv;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern int	AB;
	extern int	OUT_DIV_CURL;
	extern FILE	*FP; 


	/* timing */
	time1 = MPI_Wtime();


	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){

				/* compute derivatives and divergence */
				divv  = (v[i+1][j][k].x-v[i][j][k].x)*dx[i];
				divv += (v[i][j+1][k].y-v[i][j][k].y)*dy[j];
				divv += (v[i][j][k+1].z-v[i][j][k].z)*dz[k];
		
				/* update pressure (p) */
				p[i][j][k] -= pi[i][j][k]*divv;

				/* update div */
				if (OUT_DIV_CURL){
					div[i][j][k] = divv;
				}
			}
		}
	}


	/* absorbing frame */
	if (AB){
		for (i=nx1;i<=nx2;i++){
			for (j=ny1;j<=ny2;j++){
				for (k=nz1;k<=nz2;k++){
					p[i][j][k] *= absorb_coeff[i][j][k];
				}
			}
		}
	}

	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for stress update (elastic): \t %4.2f s.\n",time);

	return time;
}
