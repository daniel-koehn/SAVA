/*------------------------------------------------------------------------
 *   compute FD-coefficients for a staggered grid finite difference scheme of 2nd order accuracy in space
 *   (and second order accuracy in time)
 *   and compute reciprocal value of density
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void fd_coeff(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, 
float * x, float * y, float * z, float * xp, float * yp, float * zp,
float * dx, float *dxp, float * dy, float *dyp, float * dz, float * dzp,
float *** pi, float *** u, float *** uipjk, float *** uijpk, float *** uijkp, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk){

	/* local variables */
	int	i, j, k;

	/* extern variables */
	extern float	DT;
	extern int	MYID;
	extern FILE	*FP;


	/* log information */
	fprintf(FP,"\n **Message from fd_coeff (printed by PE %d):\n",MYID);
	fprintf(FP," Computing FD-coefficients ... \n");

	/* coefficients for spatial derivatives in x-direction */
	for (i=nx1;i<=nx2;i++){
		dx[i]  = 1.0/(x[i+1]-x[i]);
		dxp[i] = 1.0/(xp[i]-xp[i-1]);
	}

	/* coefficients for spatial derivatives in y-direction */
	for (j=ny1;j<=ny2;j++){
		dy[j]  = 1.0/(y[j+1]-y[j]);
		dyp[j] = 1.0/(yp[j]-yp[j-1]);
	}

	/* coefficients for spatial derivatives in z-direction */
	for (k=nz1;k<=nz2;k++){
		dz[k]  = 1.0/(z[k+1]-z[k]);
		dzp[k] = 1.0/(zp[k]-zp[k-1]);
	}

	/* save inverse density values */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				   pi[i][j][k] = DT*   pi[i][j][k];
				    u[i][j][k] = DT*    u[i][j][k];
				uipjk[i][j][k] = DT*uipjk[i][j][k];
				uijpk[i][j][k] = DT*uijpk[i][j][k];
				uijkp[i][j][k] = DT*uijkp[i][j][k];

				if (rhoijpkp[i][j][k] != 0.0)
					rhoijpkp[i][j][k] = DT/rhoijpkp[i][j][k];
				if (rhoipjkp[i][j][k] != 0.0)
					rhoipjkp[i][j][k] = DT/rhoipjkp[i][j][k];
				if (rhoipjpk[i][j][k] != 0.0)
					rhoipjpk[i][j][k] = DT/rhoipjpk[i][j][k];
			}
		}
	}


	/* log information */
	fprintf(FP," Computing FD-coefficients finished \n\n");

}
