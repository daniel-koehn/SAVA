/*------------------------------------------------------------------------
 *   store amplitudes (particle velocities and/or pressure and/or curl and div) 
 *    at receiver positions into arrays
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo(int lsamp, int ntr_loc, int **recpos_loc, 
float **sectionvx, float **sectionvy, float **sectionvz, float **sectionax, float **sectionay, float ** sectionaz, 
float **sectiondiv, float **sectioncurlx, float **sectioncurly, float **sectioncurlz, float **sectionp,
float **sectiontxx, float **sectiontxy, float ** sectiontxz, float ** sectiontyy, float ** sectiontyz, float ** sectiontzz, 
struct vector3d ***v, struct tensor3d ***t, struct vector3d ***a, struct divcurl3d ***w){

	/* extern variables */
	extern int	NDT, SEISMO;

	/* local variables */
	int	itr, ins, nxrec, nyrec, nzrec;


	ins = lsamp/NDT;
	for (itr=1;itr<=ntr_loc;itr++){
		nxrec = recpos_loc[itr][1];
		nyrec = recpos_loc[itr][2];
		nzrec = recpos_loc[itr][3];

		if (SEISMO & 16){
			/* stress components */
			sectiontxx[itr][ins] = t[nxrec][nyrec][nzrec].xx;	/* (xp,yp,zp) */
			sectiontyy[itr][ins] = t[nxrec][nyrec][nzrec].yy;	/* (xp,yp,zp) */
			sectiontzz[itr][ins] = t[nxrec][nyrec][nzrec].zz;	/* (xp,yp,zp) */
			sectiontyz[itr][ins] = t[nxrec][nyrec][nzrec].yz;	/* (xp,y,z) */
			sectiontxz[itr][ins] = t[nxrec][nyrec][nzrec].xz;	/* (x,yp,z) */
			sectiontxy[itr][ins] = t[nxrec][nyrec][nzrec].xy;	/* (x,y,zp) */
		}
		if (SEISMO & 8){
			/* particle acceleration */
			sectionax[itr][ins] = a[nxrec][nyrec][nzrec].x;	/* (x,yp,zp) */
			sectionay[itr][ins] = a[nxrec][nyrec][nzrec].y;	/* (xp,y,zp) */
			sectionaz[itr][ins] = a[nxrec][nyrec][nzrec].z;	/* (xp,yp,z) */
		}
		if (SEISMO & 4){
			/* div and curl */
			sectiondiv[itr][ins]   = w[nxrec][nyrec][nzrec].div;	/* (xp,yp,zp) */
			sectioncurlx[itr][ins] = w[nxrec][nyrec][nzrec].curlx;	/* (xp,y,z) */
			sectioncurly[itr][ins] = w[nxrec][nyrec][nzrec].curly;	/* (x,yp,z) */
			sectioncurlz[itr][ins] = w[nxrec][nyrec][nzrec].curlz;	/* (x,y,zp) */
		}
		if (SEISMO & 2){
			/* pressure */
			sectionp[itr][ins] = (-t[nxrec][nyrec][nzrec].xx-t[nxrec][nyrec][nzrec].yy-t[nxrec][nyrec][nzrec].zz)/3.0;	/* (xp,yp,zp) */
		}
		if (SEISMO & 1){
			/* particle velocity */     
			sectionvx[itr][ins] = v[nxrec][nyrec][nzrec].x;	/* (x,yp,zp) */
			sectionvy[itr][ins] = v[nxrec][nyrec][nzrec].y;	/* (xp,y,zp) */
			sectionvz[itr][ins] = v[nxrec][nyrec][nzrec].z;	/* (xp,yp,z) */
		}
	}
}
