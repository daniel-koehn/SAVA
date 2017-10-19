/*------------------------------------------------------------------------
 *   updating particle velocities 
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_v(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
	struct vector3d ***v, float *** p, struct vector3d ***a, 
	float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk,
	float *** absorb_coeff,
	float * dxp, float * dyp, float * dzp, int infoout){

	/* local variables */
	int	i, j, k;
	float	a1, a2, a3;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern float	DT;
	extern int	AB;
	extern int	OUT_ACCEL;
	extern FILE	*FP;


	/* timing */
	time1 = MPI_Wtime();
	

	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				/* compute derivatives */
				a1  = -(p[i][j][k]-p[i-1][j][k])*dxp[i];
				a1 *= rhoijpkp[i][j][k];

				a2  = -(p[i][j][k]-p[i][j-1][k])*dyp[j];
				a2 *= rhoipjkp[i][j][k];

				a3  = -(p[i][j][k]-p[i][j][k-1])*dzp[k];
				a3 *= rhoipjpk[i][j][k];

				/* update particle velocity (vx, vy, vz) */
				v[i][j][k].x += a1;
				v[i][j][k].y += a2;
				v[i][j][k].z += a3;

				/* update components of particle acceleration (ax, ay, az) */
				if (OUT_ACCEL){
					a[i][j][k].x = a1/DT;
					a[i][j][k].y = a2/DT;
					a[i][j][k].z = a3/DT;
				}
			}
		}
	}

	/* modify wavefield at right model boundary to fulfill free surface BC */
	/*if ((POS[1]==NPROCX-1)&&(!(PERIODIC_LR))){
		for (k=nz1;k<=nz2;k++){
			a1  = p[nx2][j][k]*dxp[nx2+1];     dxp[nx2+1] not computed before
			a1 *= 2.0*rho[nx2][j][k];
			v[nx2+1][j][k].x += DT*a1;
		}
	}*/
	/* modify wavefield at front model boundary to fulfill free surface BC */
	/*if ((POS[2]==NPROCY-1)&&(!(PERIODIC_BF))){
		for (k=nz1;k<=nz2;k++){
			a2  = p[i][ny2][k]*dyp[ny2+1];     dyp[ny2+1] not computed before
			a2 *= 2.0*rho[i][ny2][k];
			v[i][ny2+1][k].y += DT*a2;
		}
	}*/
	/* modify wavefield at bottom model boundary to fulfill free surface BC */
	/*if ((POS[3]==NPROCZ-1)&&(!(PERIODIC_TB))){
		for (i=nx1;i<=nx2;i++){
			a3  = p[i][j][nz2]*dzp[nz2+1];     dzp[nz2+1] not computed before
			a3 *= 2.0*rho[i][j][nz2];
			v[i][j][nz2+1].z += DT*a3;
		}
	}*/


	/* absorbing frame */
	if (AB){
		for (i=nx1;i<=nx2;i++){
			for (j=ny1;j<=ny2;j++){
				for (k=nz1;k<=nz2;k++){
					v[i][j][k].x *= absorb_coeff[i][j][k];
					v[i][j][k].y *= absorb_coeff[i][j][k];
					v[i][j][k].z *= absorb_coeff[i][j][k];
				}
			}
		}
	}

	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for particle velocity update: \t %4.2f s.\n",time);

	return time;
}
