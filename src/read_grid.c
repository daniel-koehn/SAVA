/*------------------------------------------------------------------------
 *   compute receiver positions or read them from external file        
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void read_grid(FILE *fp, char *dx_file, int dim, float *x, float *xp, float *xg, float *xpg){

	/* extern variables */
	extern int	NX[3], NXG[3];
	extern int	POS[3], MYID;
	extern float	REFMOD[3], DX, DY, DZ;
	extern char	AXIS[3];

	/* local variables */
	FILE		*fpx;
	int		n, nd, nsum, i, j;
	const int	nxg1 = NXG[dim]+1, nxpos = POS[dim]*NX[dim];
	int		c=0;
	float		sum;
	float		d;


	fprintf(fp,"\n **Message from read_grid (printed by PE %d):\n",MYID);
	if (!(MYID)){
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
		xg[0]     =  REFMOD[dim] +  xg[NXG[dim]]- xg[nxg1];
		xpg[0]    =  REFMOD[dim] + xpg[NXG[dim]]- xg[nxg1];
		xpg[nxg1] = -REFMOD[dim] +  xg[nxg1]    +xpg[1];
		
		fprintf(fp,"\n Grid-spacing in %c-direction read successfully\n",AXIS[dim]);
		fprintf(fp,"\n Grid point positions (index / %c / %cp): \n",(char)(AXIS[dim]),AXIS[dim]);
		for (i=1;i<=NXG[dim];i++){
			fprintf(fp,"%d   %f   %f\n",i,xg[i],xpg[i]);
		}
		fprintf(fp,"\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast( &xg[0],NXG[dim]+2,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&xpg[0],NXG[dim]+2,MPI_FLOAT,0,MPI_COMM_WORLD);

	/* compute local vectors for x- and z-grid point positions */
	for (i=0;i<=NX[dim]+1;i++){
		x[i]  =  xg[nxpos + i];
		xp[i] = xpg[nxpos + i];
	}
	fprintf(fp,"\n Computation of local grid positions finished.\n\n");

}
