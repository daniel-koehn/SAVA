/*------------------------------------------------------------------------
 *   updating strain tensor components, div and curl of particle velocity
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *
 *   D. Köhn, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


double update_e_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
struct vector3d ***v, struct tensor3d ***e, struct divcurl3d ***w,
float * dx, float * dxp, float * dy, float * dyp, float * dz, float * dzp, 
int infoout){

/* local variables */
	int	i, j, k;
	float	a1, a2, a3, a4, a5, a6, a7, a8, a9;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern int	OUT_DIV_CURL;
	extern FILE	*FP;


	/* timing */
	time1 = MPI_Wtime();


	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){

				/* compute derivatives and divergence */
				a1 = (v[i+1][j][k].x-v[i][j][k].x)*dx[i];
				a2 = (v[i][j+1][k].y-v[i][j][k].y)*dy[j];
				a3 = (v[i][j][k+1].z-v[i][j][k].z)*dz[k];
				a4 = (v[i][j][k].x-v[i][j-1][k].x)*dyp[j];
				a7 = (v[i][j][k].x-v[i][j][k-1].x)*dzp[k];
				a5 = (v[i][j][k].y-v[i-1][j][k].y)*dxp[i];
				a8 = (v[i][j][k].y-v[i][j][k-1].y)*dzp[k];
				a6 = (v[i][j][k].z-v[i-1][j][k].z)*dxp[i];
				a9 = (v[i][j][k].z-v[i][j-1][k].z)*dyp[j];

				/* update strain tensor (exx, eyy, ezz, 2eyz, 2exz, 2exy) */
				e[i][j][k].xx = a1;
				e[i][j][k].yy = a2;
				e[i][j][k].zz = a3;
				e[i][j][k].yz = a8+a9;
				e[i][j][k].xz = a6+a7;
				e[i][j][k].xy = a4+a5;

				/* update div, curlx, curly and curlz */
				if (OUT_DIV_CURL){
					w[i][j][k].div   = a1+a2+a3;
					w[i][j][k].curlx = a9-a8;
					w[i][j][k].curly = a7-a6;
					w[i][j][k].curlz = a5-a4;
				}
			}
		}
	}

	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for strain update (elastic): \t %4.2f s.\n",time);

	return time;
}
