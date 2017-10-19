/*------------------------------------------------------------------------
 *   updating stress tensor components
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *
 *   D. Köhn, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "fd3dtricl.h"

double update_s_el(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
	struct tensor3d ***e, struct tensor3d ***t,
	float *** c1111,  float *** c1122,  float *** c1133,  float *** c1123,  float *** c1113,   float *** c1112, 
			  float *** c2222,  float *** c2233,  float *** c2223,  float *** c1322,   float *** c1222, 
					    float *** c3333,  float *** c2333,  float *** c1333,   float *** c1233, 
	float *** c1123h, float *** c2223h, float *** c2333h, float *** c2323h, float *** c1323ha, float *** c1223ha, 
	float *** c1113h, float *** c1322h, float *** c1333h, float *** c1323h, float *** c1313h,  float *** c1213ha, 
	float *** c1112h, float *** c1222h, float *** c1233h, float *** c1223h, float *** c1213h,  float *** c1212h, 
	float *** absorb_coeff,
	float * wxp1, float * wx1, float * wxp2, float * wx2, 
	float * wyp1, float * wy1, float * wyp2, float * wy2, 
	float * wzp1, float * wz1, float * wzp2, float * wz2,
	int infoout){

	/* local variables */
	int	i, j, k;
	float	wa, wb, wc, wd, ea;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern int	AB;
	extern FILE	*FP; 


	/* timing */
	time1 = MPI_Wtime();


	/* update stress tensor (txx, tyy, tzz, tyz, txz, txy) with 
	   strain components that do not require averaging (update for orthorhombic medium) */
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
	
	/* averaging of strain components (exx, eyy, ezz, 2eyz, 2exz, 2exy) and 
	   update stress tensor (txx, tyy, tzz, tyz, txz, txy) with strain components 
	   that require averaging (additional update for triclinic medium) */
	//av_wavefield(e.yz, ea, 15, wxp1, wxp2, wy1, wy2, wz1, wz2); /* left, front */
	for (i=nx1;i<=nx2;i++){
		wa = wxp1[i];
		wb = wxp2[i];
		for (j=ny1;j<=ny2;j++){
			wc  = wy2[j]*wb;
			wd  = wy2[j]*wa;
			wa *= wy1[j];
			wb *= wy1[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wa*e[i][j-1][k].yz;
				ea += wb*e[i+1][j-1][k].yz;
				ea += wc*e[i+1][j][k].yz;
				ea += wd*e[i][j][k].yz;

				t[i][j][k].xz += ea*c1323h[i][j][k];
			}
		}
	}
	//av_wavefield(e.yz, ea, 14, wxp1, wxp2, wy1, wy2, wz1, wz2); /* left, down */
	for (i=nx1;i<=nx2;i++){
		wa = wxp1[i];
		wb = wxp2[i];
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				ea  = wz1[k]*wa*e[i][j][k-1].yz;
				ea += wz1[k]*wb*e[i+1][j][k-1].yz;
				ea += wz2[k]*wb*e[i+1][j][k].yz;
				ea += wz2[k]*wa*e[i][j][k].yz;

				t[i][j][k].xy += ea*c1223h[i][j][k];
			}
		}
	}
	
	//av_wavefield(e.xz, ea, 18, wx1, wx2, wyp1, wyp2, wz1, wz2); /* right, back */
	for (i=nx1;i<=nx2;i++){
		wa = wx1[i];
		wb = wx2[i];
		for (j=ny1;j<=ny2;j++){
			wc  = wy2[j]*wb;
			wd  = wy2[j]*wa;
			wa *= wy1[j];
			wb *= wy1[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wa*e[i-1][j][k].xz;
				ea += wb*e[i][j][k].xz;
				ea += wc*e[i][j+1][k].xz;
				ea += wd*e[i-1][j+1][k].xz;

				t[i][j][k].yz += ea*c1323ha[i][j][k]; /*???*/
			}
		}
	}
	//av_wavefield(e.xz, ea, 13, wx1, wx2, wyp1, wyp2, wz1, wz2); /* back, down */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			wa = wyp1[j];
			wb = wyp2[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wz1[k]*wa*e[i][j][k-1].xz;
				ea += wz1[k]*wb*e[i][j+1][k-1].xz;
				ea += wz2[k]*wb*e[i][j+1][k].xz;
				ea += wz2[k]*wa*e[i][j][k].xz;

				t[i][j][k].xy += ea*c1213h[i][j][k];
			}
		}
	}

	//av_wavefield(e.xy, ea, 17, wx1, wx2, wy1, wy2, wzp1, wzp2); /* right, up */
	for (i=nx1;i<=nx2;i++){
		wa = wx1[i];
		wb = wx2[i];
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				ea  = wzp1[k]*wa*e[i-1][j][k].xy;
				ea += wzp1[k]*wb*e[i][j][k].xy;
				ea += wzp2[k]*wb*e[i][j][k+1].xy;
				ea += wzp2[k]*wa*e[i-1][j][k+1].xy;

				t[i][j][k].yz += ea*c1223ha[i][j][k]; /*???*/
			}
		}
	}
	//av_wavefield(e.xy, ea, 16, wx1, wx2, wy1, wy2, wzp1, wzp2); /* front, up */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			wa = wy1[j];
			wb = wy2[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wzp1[k]*wa*e[i][j-1][k].xy;
				ea += wzp1[k]*wb*e[i][j][k].xy;
				ea += wzp2[k]*wb*e[i][j][k+1].xy;
				ea += wzp2[k]*wa*e[i][j-1][k+1].xy;

				t[i][j][k].xz += ea*c1213ha[i][j][k]; /*???*/
			}
		}
	}

	//av_wavefield(e.xx, ea, 7, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* back, up */
	//av_wavefield(e.yy, ea, 7, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* back, up */
	//av_wavefield(e.zz, ea, 7, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* back, up */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			wa = wyp1[j];
			wb = wyp2[j];
			for (k=nz1;k<=nz2;k++){
				wc  = wzp2[k]*wb;
				wd  = wzp2[k]*wa;
				wa *= wzp1[k];
				wb *= wzp1[k];

				ea  = wa*e[i][j][k].xx;
				ea += wb*e[i][j+1][k].xx;
				ea += wc*e[i][j+1][k+1].xx;
				ea += wd*e[i][j][k+1].xx;

				t[i][j][k].yz += ea*c1123h[i][j][k];

				ea  = wa*e[i][j][k].yy;
				ea += wb*e[i][j+1][k].yy;
				ea += wc*e[i][j+1][k+1].yy;
				ea += wd*e[i][j][k+1].yy;

				t[i][j][k].yz += ea*c2223h[i][j][k];

				ea  = wa*e[i][j][k].zz;
				ea += wb*e[i][j+1][k].zz;
				ea += wc*e[i][j+1][k+1].zz;
				ea += wd*e[i][j][k+1].zz;

				t[i][j][k].yz += ea*c2333h[i][j][k];
			}
		}
	}
	//av_wavefield(e.xx, ea, 8, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* left, up */
	//av_wavefield(e.yy, ea, 8, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* left, up */
	//av_wavefield(e.zz, ea, 8, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* left, up */
	for (i=nx1;i<=nx2;i++){
		wa = wxp1[i];
		wb = wxp2[i];
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				wc  = wzp2[k]*wb;
				wd  = wzp2[k]*wa;
				wa *= wzp1[k];
				wb *= wzp1[k];
			  
				ea  = wa*e[i][j][k].xx;
				ea += wb*e[i+1][j][k].xx;
				ea += wc*e[i+1][j][k+1].xx;
				ea += wd*e[i][j][k+1].xx;

				t[i][j][k].xz += ea*c1113h[i][j][k];

				ea  = wa*e[i][j][k].yy;
				ea += wb*e[i+1][j][k].yy;
				ea += wc*e[i+1][j][k+1].yy;
				ea += wd*e[i][j][k+1].yy;

				t[i][j][k].xz += ea*c1322h[i][j][k];

				ea  = wa*e[i][j][k].zz;
				ea += wb*e[i+1][j][k].zz;
				ea += wc*e[i+1][j][k+1].zz;
				ea += wd*e[i][j][k+1].zz;

				t[i][j][k].xz += ea*c1333h[i][j][k];
			}
		}
	}
	//av_wavefield(e.xx, ea, 9, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* left, back */
	//av_wavefield(e.yy, ea, 9, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* left, back */
	//av_wavefield(e.zz, ea, 9, wxp1, wxp2, wyp1, wyp2, wzp1, wzp2); /* left, back */
	for (i=nx1;i<=nx2;i++){
		wa = wxp1[i];
		wb = wxp2[i];
		for (j=ny1;j<=ny2;j++){
			wc  = wyp2[j]*wb;
			wd  = wyp2[j]*wa;
			wa *= wyp1[j];
			wb *= wyp1[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wa*e[i][j][k].xx;
				ea += wb*e[i+1][j][k].xx;
				ea += wc*e[i+1][j+1][k].xx;
				ea += wd*e[i][j+1][k].xx;

				t[i][j][k].xy += ea*c1112h[i][j][k];

				ea  = wa*e[i][j][k].yy;
				ea += wb*e[i+1][j][k].yy;
				ea += wc*e[i+1][j+1][k].yy;
				ea += wd*e[i][j+1][k].yy;

				t[i][j][k].xy += ea*c1222h[i][j][k];

				ea  = wa*e[i][j][k].zz;
				ea += wb*e[i+1][j][k].zz;
				ea += wc*e[i+1][j+1][k].zz;
				ea += wd*e[i][j+1][k].zz;

				t[i][j][k].xy += ea*c1233h[i][j][k];
			}
		}
	}

	//av_wavefield(e.yz, ea, 10, wx1, wx2, wy1, wy2, wz1, wz2); /* front, down */
	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			wa = wy1[j];
			wb = wy2[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wz1[k]*wa*e[i][j-1][k-1].yz;
				ea += wz1[k]*wb*e[i][j][k-1].yz;
				ea += wz2[k]*wb*e[i][j][k].yz;
				ea += wz2[k]*wa*e[i][j-1][k].yz;

				t[i][j][k].xx += ea*c1123[i][j][k];
				t[i][j][k].yy += ea*c2223[i][j][k];
				t[i][j][k].zz += ea*c2333[i][j][k];
			}
		}
	}

	//av_wavefield(e.xz, ea, 11, wx1, wx2, wy1, wy2, wz1, wz2); /* right, down */
	for (i=nx1;i<=nx2;i++){
		wa = wx1[i];
		wb = wx2[i];
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				ea  = wz1[k]*wa*e[i-1][j][k-1].xz;
				ea += wz1[k]*wb*e[i][j][k-1].xz;
				ea += wz2[k]*wb*e[i][j][k].xz;
				ea += wz2[k]*wa*e[i-1][j][k].xz;

				t[i][j][k].xx += ea*c1113[i][j][k];
				t[i][j][k].yy += ea*c1322[i][j][k];
				t[i][j][k].zz += ea*c1333[i][j][k];
			}
		}
	}

	//av_wavefield(e.xy, ea, 12, wx1, wx2, wy1, wy2, wz1, wz2); /* right, front */
	for (i=nx1;i<=nx2;i++){
		wa = wx1[i];
		wb = wx2[i];
		for (j=ny1;j<=ny2;j++){
			wc  = wy2[j]*wb;
			wd  = wy2[j]*wa;
			wa *= wy1[j];
			wb *= wy1[j];
			for (k=nz1;k<=nz2;k++){
				ea  = wa*e[i-1][j-1][k].xy;
				ea += wb*e[i][j-1][k].xy;
				ea += wc*e[i][j][k].xy;
				ea += wd*e[i-1][j][k].xy;

				t[i][j][k].xx += ea*c1112[i][j][k];
				t[i][j][k].yy += ea*c1222[i][j][k];
				t[i][j][k].zz += ea*c1233[i][j][k];
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
					t[i][j][k].yz *= absorb_coeff[i][j][k];
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
