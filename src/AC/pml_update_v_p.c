/*------------------------------------------------------------------------
 *   Non-splitting CFS-PML boundary condition at the numerical edge of the grid,
 *   Update of particle velocities
 *
 *   F.H. Drossaert, A. Giannopoulos (2007), A nonsplit complex frequency-shifted
 *   PML based on recursive integration for FDTD modeling of elastic waves. 
 *   Geophysics 72, T9-T17.
 *
 *   F.H. Drossaert, A. Giannopoulos (2007), Complex frequency shifted
 *   convolution PML for FDTD modeling of elastic waves. Wave Motion 44, 593-604.
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double pml_update_v(struct vector3d ***v, float *** p, 
	float *** rhoijpkp, float *** rhoipjkp, float *** rhoipjpk, 
	struct pml *pmlle, struct pml *pmlri, struct pml *pmlba, struct pml *pmlfr, struct pml *pmlto, struct pml *pmlbo, 
	float ***Pp_x_l, float ***Pp_x_r, float ***Pp_y_b, float ***Pp_y_f, float ***Pp_z_t, float ***Pp_z_b, 
	float *dxp, float *dyp, float *dzp, int infoout){

	/* local variables */
	int	i, j, k, l;
	float	p_x, p_y, p_z;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern int	NX[3];
	extern int	PML_LE, PML_RI, PML_BA, PML_FR, PML_TO, PML_BO, FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern FILE	*FP;


	/* timing */
	time1 = MPI_Wtime();


	/* PML update at top boundary */
	if (PML_TO){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=FWTO;k++){
					l = FWTO-k+1;

					/* calculate derivatives */
					p_z = -(p[i][j][k]-p[i][j][k-1])*dzp[k];

					/* calculate convolution operator for PML */
					Pp_z_t[i][j][k] *= pmlto[l].exp;
					Pp_z_t[i][j][k] += p_z;

					/* Apply PML BC */
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlto[l].x1*Pp_z_t[i][j][k] + pmlto[l].x2*p_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pp_z_t[i][j][k] += p_z;
				}
			 }
		}
	}

	/* PML update at bottom boundary */
	if (PML_BO){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
				for (k=NX[2]-FWBO+1;k<=NX[2];k++){
					l = k-NX[2]+FWBO;

					/* calculate derivatives */
					p_z = -(p[i][j][k]-p[i][j][k-1])*dzp[k];

					/* calculate convolution operator for PML */
					Pp_z_b[i][j][k] *= pmlbo[l].exp;
					Pp_z_b[i][j][k] += p_z;

					/* Apply PML BC */
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlbo[l].x1*Pp_z_b[i][j][k] + pmlbo[l].x2*p_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pp_z_b[i][j][k] += p_z;
				}
			}
		}
	}

	/* PML update at back boundary */
	if (PML_BA){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=FWBA;j++){
				l = FWBA-j+1;

				for (k=1;k<=NX[2];k++){

					/* calculate derivatives */
					p_y = -(p[i][j][k]-p[i][j-1][k])*dyp[j];

					/* calculate convolution operator for PML */
					Pp_y_b[i][j][k] *= pmlba[l].exp;
					Pp_y_b[i][j][k] += p_y;

					/* Apply PML BC */
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlba[l].x1*Pp_y_b[i][j][k] + pmlba[l].x2*p_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pp_y_b[i][j][k] += p_y;
				}
			 }
		}
	}

	/* PML update at front boundary */
	if (PML_FR){
		for (i=1;i<=NX[0];i++){
			for (j=NX[1]-FWFR+1;j<=NX[1];j++){
				l = j-NX[1]+FWFR;

				for (k=1;k<=NX[2];k++){

					/* calculate derivatives */
					p_y = -(p[i][j][k]-p[i][j-1][k])*dyp[j];

					/* calculate convolution operator for PML */
					Pp_y_f[i][j][k] *= pmlfr[l].exp;
					Pp_y_f[i][j][k] += p_y;

					/* Apply PML BC */
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlfr[l].x1*Pp_y_f[i][j][k] + pmlfr[l].x2*p_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pp_y_f[i][j][k] += p_y;
				}
			 }
		}
	}

	/* PML update at left boundary */
	if (PML_LE){
		for (i=1;i<=FWLE;i++){
			l = FWLE-i+1;

			for (j=1;j<=NX[1];j++){
				for (k=1;k<=NX[2];k++){

					/* calculate derivatives */
					p_x = -(p[i][j][k]-p[i-1][j][k])*dxp[i];

					/* calculate convolution operator for PML */
					Pp_x_l[i][j][k] *= pmlle[l].exp;
					Pp_x_l[i][j][k] += p_x;
				
					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlle[l].x1*Pp_x_l[i][j][k] + pmlle[l].x2*p_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pp_x_l[i][j][k] += p_x;
				}
			}
		}
	}

	/* PML update at right boundary */
	if (PML_RI){
		for (i=NX[0]-FWRI+1;i<=NX[0];i++){
			l = i-NX[0]+FWRI;

			for (j=1;j<=NX[1];j++){
				for (k=1;k<=NX[2];k++){

					/* calculate derivatives */
					p_x = -(p[i][j][k]-p[i-1][j][k])*dxp[i];

					/* calculate convolution operator for PML */
					Pp_x_r[i][j][k] *= pmlri[l].exp;
					Pp_x_r[i][j][k] += p_x;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlri[l].x1*Pp_x_r[i][j][k] + pmlri[l].x2*p_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pp_x_r[i][j][k] += p_x;
				}
			}
		}
	}


	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for PML particle velocity update: \t %4.2f s.\n",time);

	return time;

}
