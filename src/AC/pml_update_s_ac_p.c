/*------------------------------------------------------------------------
 *   Non-splitting CFS-PML boundary condition at the numerical edge of the grid,
 *   Update of stress tensor components
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

double pml_update_s(struct vector3d ***v, float ***p, float ***pi, 
	struct pml *pmlle, struct pml *pmlri, struct pml *pmlba, struct pml *pmlfr, struct pml *pmlto, struct pml *pmlbo, 
	float ***Pvx_x_l, float ***Pvx_x_r, float ***Pvy_y_b, float ***Pvy_y_f, float ***Pvz_z_t, float ***Pvz_z_b,
	float *dx, float *dy, float *dz, int infoout){

	/* extern variables */
	extern int	NX[3];
	extern int	PML_LE, PML_RI, PML_BA, PML_FR, PML_TO, PML_BO, FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern FILE	*FP;

	/* local variables */
	int	i, j, k, l;
	float	vx_x, vy_y, vz_z;
	double	time=0, time1=0, time2=0;


	/* timing */
	time1 = MPI_Wtime();

	/* PML update at top boundary */
	if (PML_TO){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=FWTO;k++){
					l = FWTO-k+1;

					/* calculate derivatives */
					vz_z = (v[i][j][k+1].z-v[i][j][k].z)*dz[k];

					/* calculate convolution operator for PML */
					Pvz_z_t[i][j][k] *= pmlto[l].expp;
					Pvz_z_t[i][j][k] += vz_z;

					/* Apply PML BC */
					p[i][j][k] -= pi[i][j][k]*(pmlto[l].xp1*Pvz_z_t[i][j][k] + pmlto[l].xp2*vz_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pvz_z_t[i][j][k] += vz_z;
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
					vz_z = (v[i][j][k+1].z-v[i][j][k].z)*dz[k];

					/* calculate convolution operator for PML */
					Pvz_z_b[i][j][k] *= pmlbo[l].expp;
					Pvz_z_b[i][j][k] += vz_z;
					
					/* Apply PML BC */
					p[i][j][k] -= pi[i][j][k]*(pmlbo[l].xp1*Pvz_z_b[i][j][k] + pmlbo[l].xp2*vz_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pvz_z_b[i][j][k] += vz_z;
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
					vy_y = (v[i][j+1][k].y-v[i][j][k].y)*dy[j];

					/* calculate convolution operator for PML */
					Pvy_y_b[i][j][k] *= pmlba[l].expp;
					Pvy_y_b[i][j][k] += vy_y;

					/* Apply PML BC */
					p[i][j][k] -= pi[i][j][k]*(pmlba[l].xp1*Pvy_y_b[i][j][k] + pmlba[l].xp2*vy_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pvy_y_b[i][j][k] += vy_y;
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
					vy_y = (v[i][j+1][k].y-v[i][j][k].y)*dy[j];

					/* calculate convolution operator for PML */
					Pvy_y_f[i][j][k] *= pmlfr[l].expp;
					Pvy_y_f[i][j][k] += vy_y;

					/* Apply PML BC */
					p[i][j][k] -= pi[i][j][k]*(pmlfr[l].xp1*Pvy_y_f[i][j][k] + pmlfr[l].xp2*vy_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pvy_y_f[i][j][k] += vy_y;
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
					vx_x = (v[i+1][j][k].x-v[i][j][k].x)*dx[i];

					/* calculate convolution operator for PML */
					Pvx_x_l[i][j][k] *= pmlle[l].expp;
					Pvx_x_l[i][j][k] += vx_x;

					/* Apply PML BC */
					p[i][j][k] -= pi[i][j][k]*(pmlle[l].xp1*Pvx_x_l[i][j][k] + pmlle[l].xp2*vx_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pvx_x_l[i][j][k] += vx_x;
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
					vx_x = (v[i+1][j][k].x-v[i][j][k].x)*dx[i];

					/* calculate convolution operator for PML */
					Pvx_x_r[i][j][k] *= pmlri[l].expp;
					Pvx_x_r[i][j][k] += vx_x;

					/* Apply PML BC */
					p[i][j][k] -= pi[i][j][k]*(pmlri[l].xp1*Pvx_x_r[i][j][k] + pmlri[l].xp2*vx_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pvx_x_r[i][j][k] += vx_x;
				}
			}
		}
	}


	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for PML stress tensor update: \t %4.2f s.\n",time);

	return time;
}
