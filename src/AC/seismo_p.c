/*------------------------------------------------------------------------
 *   store amplitudes (particle velocities and/or pressure and/or div) 
 *    at receiver positions into arrays
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void seismo(int lsamp, int **recpos_loc, 
float **sectionvx, float **sectionvy, float **sectionvz, float **sectionp, 
float **sectionax, float **sectionay, float ** sectionaz, float **sectiondiv, 
struct vector3d ***v, float ***p, struct vector3d ***a, float ***diverg){

	/* extern variables */
	extern int	NDT, SEISMO, NTR_LOC;

	/* local variables */
	int	itr, ins, nxrec, nyrec, nzrec;


	ins = lsamp/NDT;
	for (itr=1;itr<=NTR_LOC;itr++){
		nxrec = recpos_loc[itr][1];
		nyrec = recpos_loc[itr][2];
		nzrec = recpos_loc[itr][3];

		if (SEISMO & 8){
			/* particle acceleration */
			sectionax[itr][ins] = a[nxrec][nyrec][nzrec].x;	/* (x,yp,zp) */
			sectionay[itr][ins] = a[nxrec][nyrec][nzrec].y;	/* (xp,y,zp) */
			sectionaz[itr][ins] = a[nxrec][nyrec][nzrec].z;	/* (xp,yp,z) */
		}
		if (SEISMO & 4){
			/* div */
			sectiondiv[itr][ins] = diverg[nxrec][nyrec][nzrec];	/* (xp,yp,zp) */
		}
		if (SEISMO & 2){
			/* pressure */
			sectionp[itr][ins] = p[nxrec][nyrec][nzrec];	/* (xp,yp,zp) */
		}
		if (SEISMO & 1){
			/* particle velocity */     
			sectionvx[itr][ins] = v[nxrec][nyrec][nzrec].x;	/* (x,yp,zp) */
			sectionvy[itr][ins] = v[nxrec][nyrec][nzrec].y;	/* (xp,y,zp) */
			sectionvz[itr][ins] = v[nxrec][nyrec][nzrec].z;	/* (xp,yp,z) */
		}
	}
}
