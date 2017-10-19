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
float * wxp1, float * wx1, float * wxp2, float * wx2, 
float * wyp1, float * wy1, float * wyp2, float * wy2,
float * wzp1, float * wz1, float * wzp2, float * wz2, 
float *** c1111,  float *** c1122,  float *** c1133,  float *** c1123,  float *** c1113,   float *** c1112, 
	          float *** c2222,  float *** c2233,  float *** c2223,  float *** c1322,   float *** c1222, 
	                            float *** c3333,  float *** c2333,  float *** c1333,   float *** c1233, 
float *** c1123h, float *** c2223h, float *** c2333h, float *** c2323h, float *** c1323ha, float *** c1223ha, 
float *** c1113h, float *** c1322h, float *** c1333h, float *** c1323h, float *** c1313h,  float *** c1213ha, 
float *** c1112h, float *** c1222h, float *** c1233h, float *** c1223h, float *** c1213h,  float *** c1212h, 
float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk){

	/* local variables */
	int	i, j, k;
//	float	dx1, dx2, dx3;
//	float	dy1, dy2, dy3;
//	float	dz1, dz2, dz3;
	float	dx2, dx3;
	float	dy2, dy3;
	float	dz2, dz3;

	/* extern variables */
	extern float	DT;
	extern int	MYID;
	extern FILE	*FP;


	/* log information */
	fprintf(FP,"\n **Message from fd_coeff (printed by PE %d):\n",MYID);
	fprintf(FP," Computing FD-coefficients ... \n");

	/* coefficients for spatial derivatives and wavefield averaging in x-direction */
	for (i=nx1;i<=nx2;i++){
		dx[i]  = 1.0/(x[i+1]-x[i]);
		dxp[i] = 1.0/(xp[i]-xp[i-1]);

		//dx1 = (x[i]-xp[i-1]);
		dx2 = (xp[i]-x[i]);
		dx3 = (x[i+1]-xp[i]);

		wxp1[i] = dx2*dxp[i];
		wx1[i]  = dx3*dx[i];
		wxp2[i] = 1.0 - wxp1[i];
		wx2[i]  = 1.0 - wx1[i];
		//wxp2[i] = dx1*dxp[i];
		//wx2[i]  = dx2*dx[i];
	}

	/* coefficients for spatial derivatives and wavefield averaging in y-direction */
	for (j=ny1;j<=ny2;j++){
		dy[j]  = 1.0/(y[j+1]-y[j]);
		dyp[j] = 1.0/(yp[j]-yp[j-1]);

		//dy1 = (y[j]-yp[j-1]);
		dy2 = (yp[j]-y[j]);
		dy3 = (y[j+1]-yp[j]);

		wyp1[j] = dy2*dyp[j];
		wy1[j]  = dy3*dy[j];
		wyp2[j] = 1.0 - wyp1[j];
		wy2[j]  = 1.0 - wy1[j];
		//wyp2[j] = dy1*dyp[j];
		//wy2[j]  = dy2*dy[j];
	}

	/* coefficients for spatial derivatives and wavefield averaging in z-direction */
	for (k=nz1;k<=nz2;k++){
		dz[k]  = 1.0/(z[k+1]-z[k]);
		dzp[k] = 1.0/(zp[k]-zp[k-1]);
	
		//dz1 = (z[k]-zp[k-1]);
		dz2 = (zp[k]-z[k]);
		dz3 = (z[k+1]-zp[k]);

		wzp1[k] = dz2*dzp[k];
		wz1[k]  = dz3*dz[k];
		wzp2[k] = 1.0 - wzp1[k];
		wz2[k]  = 1.0 - wz1[k];
		//wzp2[k] = dz1*dzp[k];
		//wz2[k]  = dz2*dz[k];
	}

	/* save inverse density values */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				c1111[i][j][k]   = DT*c1111[i][j][k];
				c1122[i][j][k]   = DT*c1122[i][j][k];
				c1133[i][j][k]   = DT*c1133[i][j][k];
				c1123[i][j][k]   = DT*c1123[i][j][k];
				c1113[i][j][k]   = DT*c1113[i][j][k];
				c1112[i][j][k]   = DT*c1112[i][j][k];
				c2222[i][j][k]   = DT*c2222[i][j][k];
				c2233[i][j][k]   = DT*c2233[i][j][k];
				c2223[i][j][k]   = DT*c2223[i][j][k];
				c1322[i][j][k]   = DT*c1322[i][j][k];
				c1222[i][j][k]   = DT*c1222[i][j][k];
				c3333[i][j][k]   = DT*c3333[i][j][k];
				c2333[i][j][k]   = DT*c2333[i][j][k];
				c1333[i][j][k]   = DT*c1333[i][j][k];
				c1233[i][j][k]   = DT*c1233[i][j][k];
				c1123h[i][j][k]  = DT*c1123h[i][j][k];
				c2223h[i][j][k]  = DT*c2223h[i][j][k];
				c2333h[i][j][k]  = DT*c2333h[i][j][k];
				c2323h[i][j][k]  = DT*c2323h[i][j][k];
				c1323ha[i][j][k] = DT*c1323ha[i][j][k];
				c1223ha[i][j][k] = DT*c1223ha[i][j][k];
				c1113h[i][j][k]  = DT*c1113h[i][j][k];
				c1322h[i][j][k]  = DT*c1322h[i][j][k];
				c1333h[i][j][k]  = DT*c1333h[i][j][k];
				c1323h[i][j][k]  = DT*c1323h[i][j][k];
				c1313h[i][j][k]  = DT*c1313h[i][j][k];
				c1213ha[i][j][k] = DT*c1213ha[i][j][k];
				c1112h[i][j][k]  = DT*c1112h[i][j][k];
				c1222h[i][j][k]  = DT*c1222h[i][j][k];
				c1233h[i][j][k]  = DT*c1233h[i][j][k];
				c1223h[i][j][k]  = DT*c1223h[i][j][k];
				c1213h[i][j][k]  = DT*c1213h[i][j][k];
				c1212h[i][j][k]  = DT*c1212h[i][j][k];

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
