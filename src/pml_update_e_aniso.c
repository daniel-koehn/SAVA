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

double pml_update_e(struct vector3d ***v, struct tensor3d ***e, 
	struct pml *pmlle, struct pml *pmlri, struct pml *pmlba, struct pml *pmlfr, struct pml *pmlto, struct pml *pmlbo, 
	struct vector3d ***Pv_x_l, struct vector3d ***Pv_x_r, struct vector3d ***Pv_y_b, struct vector3d ***Pv_y_f, struct vector3d ***Pv_z_t, struct vector3d ***Pv_z_b,
	float *dx, float *dxp, float *dy, float *dyp, float *dz, float *dzp, int infoout){

	/* extern variables */
	extern int	NX[3];
	extern int	PML_LE, PML_RI, PML_BA, PML_FR, PML_TO, PML_BO, FWLE, FWRI, FWBA, FWFR, FWTO, FWBO;
	extern FILE	*FP;

	/* local variables */
	int	i, j, k, l;
	float	vx_x, vy_x, vz_x, vx_y, vy_y, vz_y, vx_z, vy_z, vz_z;
	double	time=0, time1=0, time2=0;


	/* timing */
	time1 = MPI_Wtime();

	/* PML update at left boundary */
	if (PML_LE){
		for (i=1;i<=FWLE;i++){
			l = FWLE-i+1;

			for (j=1;j<=NX[1];j++){
				for (k=1;k<=NX[2];k++){

					/* calculate derivatives */
					vx_x = (v[i+1][j][k].x-v[i][j][k].x)*dx[i];
					vy_x = (v[i][j][k].y-v[i-1][j][k].y)*dxp[i];
					vz_x = (v[i][j][k].z-v[i-1][j][k].z)*dxp[i];

					/* calculate convolution operator for PML */
					Pv_x_l[i][j][k].x *= pmlle[l].expp;
					Pv_x_l[i][j][k].x += vx_x;
					Pv_x_l[i][j][k].y *= pmlle[l].exp;
					Pv_x_l[i][j][k].y += vy_x;
					Pv_x_l[i][j][k].z *= pmlle[l].exp;
					Pv_x_l[i][j][k].z += vz_x;

					/* Apply PML BC */
					e[i][j][k].xx += (pmlle[l].xp1*Pv_x_l[i][j][k].x + pmlle[l].xp2*vx_x);
					e[i][j][k].xz += (pmlle[l].x1 *Pv_x_l[i][j][k].z + pmlle[l].x2 *vz_x);
					e[i][j][k].xy += (pmlle[l].x1 *Pv_x_l[i][j][k].y + pmlle[l].x2 *vy_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pv_x_l[i][j][k].x += vx_x;
					//Pv_x_l[i][j][k].y += vy_x;
					//Pv_x_l[i][j][k].z += vz_x;
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
					vy_x = (v[i][j][k].y-v[i-1][j][k].y)*dxp[i];
					vz_x = (v[i][j][k].z-v[i-1][j][k].z)*dxp[i];

					/* calculate convolution operator for PML */
					Pv_x_r[i][j][k].x *= pmlri[l].expp;
					Pv_x_r[i][j][k].x += vx_x;
					Pv_x_r[i][j][k].y *= pmlri[l].exp;
					Pv_x_r[i][j][k].y += vy_x;
					Pv_x_r[i][j][k].z *= pmlri[l].exp;
					Pv_x_r[i][j][k].z += vz_x;

					/* Apply PML BC */
					e[i][j][k].xx += (pmlri[l].xp1*Pv_x_r[i][j][k].x + pmlri[l].xp2*vx_x);
					e[i][j][k].xz += (pmlri[l].x1 *Pv_x_r[i][j][k].z + pmlri[l].x2 *vz_x);
					e[i][j][k].xy += (pmlri[l].x1 *Pv_x_r[i][j][k].y + pmlri[l].x2 *vy_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pv_x_r[i][j][k].x += vx_x;
					//Pv_x_r[i][j][k].y += vy_x;
					//Pv_x_r[i][j][k].z += vz_x;
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
					vx_y = (v[i][j][k].x-v[i][j-1][k].x)*dyp[j];
					vy_y = (v[i][j+1][k].y-v[i][j][k].y)*dy[j];
					vz_y = (v[i][j][k].z-v[i][j-1][k].z)*dyp[j];

					/* calculate convolution operator for PML */
					Pv_y_b[i][j][k].x *= pmlba[l].exp;
					Pv_y_b[i][j][k].x += vx_y;
					Pv_y_b[i][j][k].y *= pmlba[l].expp;
					Pv_y_b[i][j][k].y += vy_y;
					Pv_y_b[i][j][k].z *= pmlba[l].exp;
					Pv_y_b[i][j][k].z += vz_y;

					/* Apply PML BC */
					e[i][j][k].yy += (pmlba[l].xp1*Pv_y_b[i][j][k].y + pmlba[l].xp2*vy_y);
					e[i][j][k].yz += (pmlba[l].x1 *Pv_y_b[i][j][k].z + pmlba[l].x2 *vz_y);
					e[i][j][k].xy += (pmlba[l].x1 *Pv_y_b[i][j][k].x + pmlba[l].x2 *vx_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pv_y_b[i][j][k].x += vx_y;
					//Pv_y_b[i][j][k].y += vy_y;
					//Pv_y_b[i][j][k].z += vz_y;
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
					vx_y = (v[i][j][k].x-v[i][j-1][k].x)*dyp[j];
					vy_y = (v[i][j+1][k].y-v[i][j][k].y)*dy[j];
					vz_y = (v[i][j][k].z-v[i][j-1][k].z)*dyp[j];

					/* calculate convolution operator for PML */
					Pv_y_f[i][j][k].x *= pmlfr[l].exp;
					Pv_y_f[i][j][k].x += vx_y;
					Pv_y_f[i][j][k].y *= pmlfr[l].expp;
					Pv_y_f[i][j][k].y += vy_y;
					Pv_y_f[i][j][k].z *= pmlfr[l].exp;
					Pv_y_f[i][j][k].z += vz_y;

					/* Apply PML BC */
					e[i][j][k].yy += (pmlfr[l].xp1*Pv_y_f[i][j][k].y + pmlfr[l].xp2*vy_y);
					e[i][j][k].yz += (pmlfr[l].x1 *Pv_y_f[i][j][k].z + pmlfr[l].x2 *vz_y);
					e[i][j][k].xy += (pmlfr[l].x1 *Pv_y_f[i][j][k].x + pmlfr[l].x2 *vx_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pv_y_f[i][j][k].x += vx_y;
					//Pv_y_f[i][j][k].y += vy_y;
					//Pv_y_f[i][j][k].z += vz_y;
				}
			}
		}
	}

	/* PML update at top boundary */
	if (PML_TO){
		for (i=1;i<=NX[0];i++){
			for (j=1;j<=NX[1];j++){
				for (k=1;k<=FWTO;k++){
					l = FWTO-k+1;

					/* calculate derivatives */
					vx_z = (v[i][j][k].x-v[i][j][k-1].x)*dzp[k];
					vy_z = (v[i][j][k].y-v[i][j][k-1].y)*dzp[k];
					vz_z = (v[i][j][k+1].z-v[i][j][k].z)*dz[k];

					/* calculate convolution operator for PML */
					Pv_z_t[i][j][k].x *= pmlto[l].exp;
					Pv_z_t[i][j][k].x += vx_z;
					Pv_z_t[i][j][k].y *= pmlto[l].exp;
					Pv_z_t[i][j][k].y += vy_z;
					Pv_z_t[i][j][k].z *= pmlto[l].expp;
					Pv_z_t[i][j][k].z += vz_z;

					/* Apply PML BC */
					e[i][j][k].zz += (pmlto[l].xp1*Pv_z_t[i][j][k].z + pmlto[l].xp2*vz_z);
					e[i][j][k].yz += (pmlto[l].x1 *Pv_z_t[i][j][k].y + pmlto[l].x2 *vy_z);
					e[i][j][k].xz += (pmlto[l].x1 *Pv_z_t[i][j][k].x + pmlto[l].x2 *vx_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pv_z_t[i][j][k].x += vx_z;
					//Pv_z_t[i][j][k].y += vy_z;
					//Pv_z_t[i][j][k].z += vz_z;
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
					vx_z = (v[i][j][k].x-v[i][j][k-1].x)*dzp[k];
					vy_z = (v[i][j][k].y-v[i][j][k-1].y)*dzp[k];
					vz_z = (v[i][j][k+1].z-v[i][j][k].z)*dz[k];

					/* calculate convolution operator for PML */
					Pv_z_b[i][j][k].x *= pmlbo[l].exp;
					Pv_z_b[i][j][k].x += vx_z;
					Pv_z_b[i][j][k].y *= pmlbo[l].exp;
					Pv_z_b[i][j][k].y += vy_z;
					Pv_z_b[i][j][k].z *= pmlbo[l].expp;
					Pv_z_b[i][j][k].z += vz_z;

					/* Apply PML BC */
					e[i][j][k].zz += (pmlbo[l].xp1*Pv_z_b[i][j][k].z + pmlbo[l].xp2*vz_z);
					e[i][j][k].yz += (pmlbo[l].x1 *Pv_z_b[i][j][k].y + pmlbo[l].x2 *vy_z);
					e[i][j][k].xz += (pmlbo[l].x1 *Pv_z_b[i][j][k].x + pmlbo[l].x2 *vx_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pv_z_b[i][j][k].x += vx_z;
					//Pv_z_b[i][j][k].y += vy_z;
					//Pv_z_b[i][j][k].z += vz_z;
				}
			}
		}
	}


	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for PML strain tensor update: \t %4.2f s.\n",time);

	return time;
}
