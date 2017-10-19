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

double pml_update_v(struct vector3d ***v, struct tensor3d ***t,
	float ***rhoijpkp, float ***rhoipjkp, float ***rhoipjpk,
	struct pml *pmlle, struct pml *pmlri, struct pml *pmlba, struct pml *pmlfr, struct pml *pmlto, struct pml *pmlbo, 
	struct vector3d ***Pt_x_l, struct vector3d ***Pt_x_r, struct vector3d ***Pt_y_b, struct vector3d ***Pt_y_f, struct vector3d ***Pt_z_t, struct vector3d ***Pt_z_b,
	float *dx, float *dxp, float *dy, float *dyp, float *dz, float *dzp, int infoout){

	/* local variables */
	int	i, j, k, l;
	float	txx_x, txy_x, txz_x, txy_y, tyy_y, tyz_y, txz_z, tyz_z, tzz_z;
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
					tzz_z = (t[i][j][k].zz-t[i][j][k-1].zz)*dzp[k];
					tyz_z = (t[i][j][k+1].yz-t[i][j][k].yz)*dz[k];
					txz_z = (t[i][j][k+1].xz-t[i][j][k].xz)*dz[k];

					/* calculate convolution operator for PML */
					Pt_z_t[i][j][k].x *= pmlto[l].expp;
					Pt_z_t[i][j][k].x += txz_z;
					Pt_z_t[i][j][k].y *= pmlto[l].expp;
					Pt_z_t[i][j][k].y += tyz_z;
					Pt_z_t[i][j][k].z *= pmlto[l].exp;
					Pt_z_t[i][j][k].z += tzz_z;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlto[l].xp1*Pt_z_t[i][j][k].x + pmlto[l].xp2*txz_z);
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlto[l].xp1*Pt_z_t[i][j][k].y + pmlto[l].xp2*tyz_z);
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlto[l].x1 *Pt_z_t[i][j][k].z + pmlto[l].x2 *tzz_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pt_z_t[i][j][k].x += txz_z;
					//Pt_z_t[i][j][k].y += tyz_z;
					//Pt_z_t[i][j][k].z += tzz_z;
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
					tzz_z = (t[i][j][k].zz-t[i][j][k-1].zz)*dzp[k];
					tyz_z = (t[i][j][k+1].yz-t[i][j][k].yz)*dz[k];
					txz_z = (t[i][j][k+1].xz-t[i][j][k].xz)*dz[k];
					
					/* calculate convolution operator for PML */
					Pt_z_b[i][j][k].x *= pmlbo[l].expp;
					Pt_z_b[i][j][k].x += txz_z;
					Pt_z_b[i][j][k].y *= pmlbo[l].expp;
					Pt_z_b[i][j][k].y += tyz_z;
					Pt_z_b[i][j][k].z *= pmlbo[l].exp;
					Pt_z_b[i][j][k].z += tzz_z;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlbo[l].xp1*Pt_z_b[i][j][k].x + pmlbo[l].xp2*txz_z);
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlbo[l].xp1*Pt_z_b[i][j][k].y + pmlbo[l].xp2*tyz_z);
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlbo[l].x1 *Pt_z_b[i][j][k].z + pmlbo[l].x2 *tzz_z);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pt_z_b[i][j][k].x += txz_z;
					//Pt_z_b[i][j][k].y += tyz_z;
					//Pt_z_b[i][j][k].z += tzz_z;
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
					tyy_y = (t[i][j][k].yy-t[i][j-1][k].yy)*dyp[j];
					tyz_y = (t[i][j+1][k].yz-t[i][j][k].yz)*dy[j];
					txy_y = (t[i][j+1][k].xy-t[i][j][k].xy)*dy[j];

					/* calculate convolution operator for PML */
					Pt_y_b[i][j][k].x *= pmlba[l].expp;
					Pt_y_b[i][j][k].x += txy_y;
					Pt_y_b[i][j][k].y *= pmlba[l].exp;
					Pt_y_b[i][j][k].y += tyy_y;
					Pt_y_b[i][j][k].z *= pmlba[l].expp;
					Pt_y_b[i][j][k].z += tyz_y;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlba[l].xp1*Pt_y_b[i][j][k].x + pmlba[l].xp2*txy_y);
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlba[l].x1 *Pt_y_b[i][j][k].y + pmlba[l].x2 *tyy_y);
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlba[l].xp1*Pt_y_b[i][j][k].z + pmlba[l].xp2*tyz_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pt_y_b[i][j][k].x += txy_y;
					//Pt_y_b[i][j][k].y += tyy_y;
					//Pt_y_b[i][j][k].z += tyz_y;
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
					tyy_y = (t[i][j][k].yy-t[i][j-1][k].yy)*dyp[j];
					tyz_y = (t[i][j+1][k].yz-t[i][j][k].yz)*dy[j];
					txy_y = (t[i][j+1][k].xy-t[i][j][k].xy)*dy[j];

					/* calculate convolution operator for PML */
					Pt_y_f[i][j][k].x *= pmlfr[l].expp;
					Pt_y_f[i][j][k].x += txy_y;
					Pt_y_f[i][j][k].y *= pmlfr[l].exp;
					Pt_y_f[i][j][k].y += tyy_y;
					Pt_y_f[i][j][k].z *= pmlfr[l].expp;
					Pt_y_f[i][j][k].z += tyz_y;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlfr[l].xp1*Pt_y_f[i][j][k].x + pmlfr[l].xp2*txy_y);
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlfr[l].x1 *Pt_y_f[i][j][k].y + pmlfr[l].x2 *tyy_y);
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlfr[l].xp1*Pt_y_f[i][j][k].z + pmlfr[l].xp2*tyz_y);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pt_y_f[i][j][k].x += txy_y;
					//Pt_y_f[i][j][k].y += tyy_y;
					//Pt_y_f[i][j][k].z += tyz_y;
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
					txx_x = (t[i][j][k].xx-t[i-1][j][k].xx)*dxp[i];
					txz_x = (t[i+1][j][k].xz-t[i][j][k].xz)*dx[i];
					txy_x = (t[i+1][j][k].xy-t[i][j][k].xy)*dx[i];

					/* calculate convolution operator for PML */
					Pt_x_l[i][j][k].x *= pmlle[l].exp;
					Pt_x_l[i][j][k].x += txx_x;
					Pt_x_l[i][j][k].y *= pmlle[l].expp;
					Pt_x_l[i][j][k].y += txy_x;
					Pt_x_l[i][j][k].z *= pmlle[l].expp;
					Pt_x_l[i][j][k].z += txz_x;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlle[l].x1 *Pt_x_l[i][j][k].x + pmlle[l].x2 *txx_x);
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlle[l].xp1*Pt_x_l[i][j][k].y + pmlle[l].xp2*txy_x);
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlle[l].xp1*Pt_x_l[i][j][k].z + pmlle[l].xp2*txz_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pt_x_l[i][j][k].x += txx_x;
					//Pt_x_l[i][j][k].y += txy_x;
					//Pt_x_l[i][j][k].z += txz_x;
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
					txx_x = (t[i][j][k].xx-t[i-1][j][k].xx)*dxp[i];
					txz_x = (t[i+1][j][k].xz-t[i][j][k].xz)*dx[i];
					txy_x = (t[i+1][j][k].xy-t[i][j][k].xy)*dx[i];

					/* calculate convolution operator for PML */
					Pt_x_r[i][j][k].x *= pmlri[l].exp;
					Pt_x_r[i][j][k].x += txx_x;
					Pt_x_r[i][j][k].y *= pmlri[l].expp;
					Pt_x_r[i][j][k].y += txy_x;
					Pt_x_r[i][j][k].z *= pmlri[l].expp;
					Pt_x_r[i][j][k].z += txz_x;

					/* Apply PML BC */
					v[i][j][k].x += rhoijpkp[i][j][k]*(pmlri[l].x1 *Pt_x_r[i][j][k].x + pmlri[l].x2 *txx_x);
					v[i][j][k].y += rhoipjkp[i][j][k]*(pmlri[l].xp1*Pt_x_r[i][j][k].y + pmlri[l].xp2*txy_x);
					v[i][j][k].z += rhoipjpk[i][j][k]*(pmlri[l].xp1*Pt_x_r[i][j][k].z + pmlri[l].xp2*txz_x);

					/* calculate convolution operator for PML (only PML version 2) */
					//Pt_x_r[i][j][k].x += txx_x;
					//Pt_x_r[i][j][k].y += txy_x;
					//Pt_x_r[i][j][k].z += txz_x;
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
