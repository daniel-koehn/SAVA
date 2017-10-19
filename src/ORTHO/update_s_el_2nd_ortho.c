/*------------------------------------------------------------------------
 *   updating stress tensor components
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *
 *   D. Köhn, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
 	struct tensor3d ***e, struct tensor3d ***t,
	float *** c1111, float *** c1122, float *** c1133, float *** c2222, float *** c2233, float *** c3333, 
	float *** c2323h, float *** c1313h, float *** c1212h, 
	float *** absorb_coeff,
	int infoout){

	/* local variables */
	int	i, j, k;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern int	AB;
	extern FILE	*FP; 


	/* timing */
	time1 = MPI_Wtime();


	/* update stress tensor (txx, tyy, tzz, tyz, txz, txy) with 
	   strain components that do not require averaging (update for orthotropic medium) */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				t[i][j][k].xx += (e[i][j][k].xx*c1111[i][j][k] + e[i][j][k].yy*c1122[i][j][k] + e[i][j][k].zz*c1133[i][j][k]);
				t[i][j][k].yy += (e[i][j][k].xx*c1122[i][j][k] + e[i][j][k].yy*c2222[i][j][k] + e[i][j][k].zz*c2233[i][j][k]);
				t[i][j][k].zz += (e[i][j][k].xx*c1133[i][j][k] + e[i][j][k].yy*c2233[i][j][k] + e[i][j][k].zz*c3333[i][j][k]);
				t[i][j][k].yz += (e[i][j][k].yz*c2323h[i][j][k]);
				t[i][j][k].xz += (e[i][j][k].xz*c1313h[i][j][k]);
				t[i][j][k].xy += (e[i][j][k].xy*c1212h[i][j][k]);   
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
