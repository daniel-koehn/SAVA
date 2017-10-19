/*------------------------------------------------------------------------
 *   updating stress tensor components
 *   by a staggered grid finite difference scheme of 2nd order accuracy in space
 *   and second order accuracy in time
 *   (viscoelastic version using one relaxation mechanism) 
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

double update_s_ve(int nx1, int nx2, int ny1, int ny2, int nz1, int nz2,
	struct vector3d ***v, struct tensor3d ***t, struct tensor3d ***r, struct divcurl3d ***w,
	float ***pi, float ***u, float ***uipjk, float ***uijpk, float ***uijkp,
	float ***taup, float ***taus, float ***tausipjk, float ***tausijpk, float ***tausijkp, float * eta, 
	float ***absorb_coeff, 
	float *dx, float *dxp, float *dy, float *dyp, float *dz, float *dzp, int infoout){

	/* local variables */
	int	i, j, k;
	float	divv, a1, a2, a3, a4, a5, a6, a7, a8, a9, b, c, d, e, f, g;
	double	time=0, time1=0, time2=0;

	/* extern variables */
	extern int	AB;
	extern int	OUT_DIV_CURL;
	extern FILE	*FP;


	/* timing */
	time1 = MPI_Wtime();


	c = (2.0-eta[1])/(2.0+eta[1]);
	g = c-1.0;

	for (i=nx1;i<=nx2;i++){
		for (j=ny1;j<=ny2;j++){
			for (k=nz1;k<=nz2;k++){
				e = 2.0*u[i][j][k]*taus[i][j][k];
				f = 2.0*u[i][j][k]*(1.0+taus[i][j][k]);

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
				d = ((pi[i][j][k]*(1.0+taup[i][j][k]))-f)*divv;

				/* update stress tensor (txx, txy, txz, tyy, tyz, tzz) */
				t[i][j][k].xx += (d + f*a1 + 0.5*r[i][j][k].xx);
				t[i][j][k].yy += (d + f*a2 + 0.5*r[i][j][k].yy);
				t[i][j][k].zz += (d + f*a3 + 0.5*r[i][j][k].zz);
				t[i][j][k].yz += (uipjk[i][j][k]*(1.0 + tausipjk[i][j][k])*(a8 + a9) + 0.5*r[i][j][k].yz);
				t[i][j][k].xz += (uijpk[i][j][k]*(1.0 + tausijpk[i][j][k])*(a6 + a7) + 0.5*r[i][j][k].xz);
				t[i][j][k].xy += (uijkp[i][j][k]*(1.0 + tausijkp[i][j][k])*(a4 + a5) + 0.5*r[i][j][k].xy);
				
				

				/* update div, curlx, curly and curlz */
				if (OUT_DIV_CURL){
					w[i][j][k].div   = divv;
					w[i][j][k].curlx = a9-a8;
					w[i][j][k].curly = a7-a6;
					w[i][j][k].curlz = a5-a4;
				}

				/* update memory variables */
				b = (pi[i][j][k]*taup[i][j][k] - e)*divv;
				r[i][j][k].xx *= c;
				r[i][j][k].xx += g*(b + e*a1);
				r[i][j][k].yy *= c;
				r[i][j][k].yy += g*(b + e*a2);
				r[i][j][k].zz *= c;
				r[i][j][k].zz += g*(b + e*a3);
				r[i][j][k].yz *= c;
				r[i][j][k].yz += g*(uipjk[i][j][k]*tausipjk[i][j][k]*(a8 + a9));
				r[i][j][k].xz *= c;
				r[i][j][k].xz += g*(uijpk[i][j][k]*tausijpk[i][j][k]*(a6 + a7));
				r[i][j][k].xy *= c;
				r[i][j][k].xy += g*(uijkp[i][j][k]*tausijkp[i][j][k]*(a4 + a5));

				/* add updated memory variables to stress tensor */
				t[i][j][k].xx += 0.5*r[i][j][k].xx;
				t[i][j][k].yy += 0.5*r[i][j][k].yy;
				t[i][j][k].zz += 0.5*r[i][j][k].zz;
				t[i][j][k].yz += 0.5*r[i][j][k].yz;
				t[i][j][k].xz += 0.5*r[i][j][k].xz;
				t[i][j][k].xy += 0.5*r[i][j][k].xy;
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
		fprintf(FP," Real time for stress update (viscoelastic): \t %4.2f s.\n",time);

	return time;
}

