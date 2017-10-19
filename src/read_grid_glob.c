/*------------------------------------------------------------------------
 *   compute receiver positions or read them from external file        
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void read_grid_glob(FILE *fp, char *dx_file, int dim, float *xpg){

	/* extern variables */
	extern int	NXG[3];
	extern float	REFMOD[3];
	extern char	AXIS[3];

	/* local variables */
	FILE		*fpx;
	int		n, nd, nsum, i, j;
	int		c=0;
	const int	nxg1 = NXG[dim]+1;
	float		*xg=NULL;
	float		sum;
	float		d;


	/* allocate memory for grid point positions */
	xg   = vector(0,nxg1);

	fprintf(fp,"\n Reading grid-spacing in %c-direction from file: \n\t%s\n",AXIS[dim],dx_file);
	fpx = fopen(dx_file,"r");
	if (fpx==NULL)
		error(" File with grid-spacing could not be opened!");

	/* find number of entries in file with grid-spacing */
	n = 0;
	while ((c=fgetc(fpx)) != EOF) 
		if (c=='\n') 
			n++;
	rewind(fpx);

	/* compute vertical vectors x and xp */
	nsum  = 0;
	sum   = REFMOD[dim];
	xg[1] = sum;
	for (i=1;i<=n;i++){
		fscanf(fpx,"%d %f\n",&nd,&d);
		for (j=nsum+1;j<=min(nsum+nd,NXG[dim]);j++){
			xpg[j]  = sum + (j-nsum-0.5)*d;
			xg[j+1] = sum + (j-nsum)*d;
		}
		nsum += nd;
		sum  += nd*d;
	}
	fclose(fpx);

	for (j=nsum+1;j<=NXG[dim];j++){
		xpg[j]  = sum + (j-nsum-0.5)*d;
		xg[j+1] = sum + (j-nsum)*d;
	}
	xg[0]     =  REFMOD[dim] +  xg[NXG[dim]] -  xg[nxg1];
	xpg[0]    =  REFMOD[dim] + xpg[NXG[dim]] -  xg[nxg1];
	xpg[nxg1] = -REFMOD[dim] +  xg[nxg1]     + xpg[1];

	/* deallocate memory for grid point positions */
	free_vector(xg,0,nxg1);

}
