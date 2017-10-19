/*------------------------------------------------------------------------
 *   updating stress tensor components
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
	struct vector3d ***v, struct tensor3d ***t, struct divcurl3d ***w,
	float ***pi, float ***u, float ***uipjk, float ***uijpk, float ***uijkp, float ***absorb_coeff,
	float *dx, float *dxp, float *dy, float *dyp, float *dz, float *dzp, int infoout){

	/* local variables */
	int	i, j, k;
	float	divv, a1, a2, a3, a4, a5, a6, a7, a8, a9, d, f;
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
				f = 2.0*u[i][j][k];

				/* compute derivatives and divergence */
				a1 = (v[i+1][j][k].x-v[i][j][k].x)*dx[i];
				a2 = (v[i][j+1][k].y-v[i][j][k].y)*dy[j];
				a3 = (v[i][j][k+1].z-v[i][j][k].z)*dz[k];
				a4 = (v[i][j][k].x-v[i][j-1][k].x)*dyp[j];
				a5 = (v[i][j][k].y-v[i-1][j][k].y)*dxp[i];
				a6 = (v[i][j][k].z-v[i-1][j][k].z)*dxp[i];
				a7 = (v[i][j][k].x-v[i][j][k-1].x)*dzp[k];
				a8 = (v[i][j][k].y-v[i][j][k-1].y)*dzp[k];
				a9 = (v[i][j][k].z-v[i][j-1][k].z)*dyp[j];

				divv = a1+a2+a3;
				d  = (pi[i][j][k]-f)*divv;

				/* update stress tensor (txx, tyy, tzz, tyz, txz, txy) */
				t[i][j][k].xx += (d + f*a1);
				t[i][j][k].yy += (d + f*a2);
				t[i][j][k].zz += (d + f*a3);
				t[i][j][k].yz += uipjk[i][j][k]*(a8 + a9);
				t[i][j][k].xz += uijpk[i][j][k]*(a6 + a7);
				t[i][j][k].xy += uijkp[i][j][k]*(a4 + a5);

				/* update div, curlx, curly and curlz */
				if (OUT_DIV_CURL){
					w[i][j][k].div = divv;
					w[i][j][k].curlx  = a9-a8;
					w[i][j][k].curly  = a7-a6;
					w[i][j][k].curlz  = a5-a4;
				}
			}
		}
	}


	/* absorbing frame */
	if (AB){
		for (i=nx1;i<=nx2;i++){
			for (j=ny1;j<=ny2;j++){
				for (k=nz1;k<=nz2;k++){
					t[i][j][k].xx *= absorb_coeff[i][j][k];
					t[i][j][k].yy *= absorb_coeff[i][j][k];
					t[i][j][k].zz *= absorb_coeff[i][j][k];
					t[i][j][k].yz *= absorb_coeff[i][j][k];
					t[i][j][k].xz *= absorb_coeff[i][j][k];
					t[i][j][k].xy *= absorb_coeff[i][j][k];
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
