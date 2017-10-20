/*------------------------------------------------------------------------
 *   compute receiver positions or read them from external file        
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

int **receiver(FILE *fp, float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg){

	/* extern variables */
	extern int	NXG[3];
	extern char	REC_FILE[STRING_SIZE];
	extern float	REFREC[3], XREC1[3], XREC2[3];
	extern int	READREC, DXREC[3];
	extern int	MYID, NTR;

	/* local variables */
	float	xrec, yrec, zrec;
	int	**recpos=NULL;
	int	itr, i, j, k;
	int	c=0;
	int	*iposx1=NULL, *iposx2=NULL, *iposy1=NULL, *iposy2=NULL, *iposz1=NULL, *iposz2=NULL;
	FILE	*fprec;


	if (!(MYID)){
		fprintf(fp,"\n **Message from function receiver (printed by PE %d):\n",MYID);
		if (READREC){     /* read receiver positions from file */
			fprintf(fp,"\n Reading receiver positions from file: \n\t%s\n",REC_FILE);
			fprec = fopen(REC_FILE,"r");
			if (fprec==NULL)
				error(" Receiver file could not be opened !");

			NTR = 0;
			while ((c=fgetc(fprec)) != EOF)
				if (c=='\n') 
					++(NTR);
			rewind(fprec);
     
			recpos = imatrix(1,NTR,1,7);
			for (itr=1;itr<=NTR;itr++){
				fscanf(fprec,"%f%f%f\n",&xrec, &yrec, &zrec);

				xrec += REFREC[0];
				yrec += REFREC[1];
				zrec += REFREC[2];
				iposx1 = irnd_to_grid(xrec,xg,xpg,NXG[0]);
				iposy1 = irnd_to_grid(yrec,yg,ypg,NXG[1]);
				iposz1 = irnd_to_grid(zrec,zg,zpg,NXG[2]);

				recpos[itr][1] = iposx1[2];
				recpos[itr][2] = iposy1[2];
				recpos[itr][3] = iposz1[2];
				recpos[itr][4] = iposx1[2];
				recpos[itr][5] = iposy1[2];
				recpos[itr][6] = iposz1[2];
				recpos[itr][7] = itr;

				free_ivector(iposx1,1,2);
				free_ivector(iposy1,1,2);
				free_ivector(iposz1,1,2);
			}

			fclose(fprec);
		}
		else{   /* receiver array with parameters from input file */
			xrec  = XREC1[0] + REFREC[0];
			yrec  = XREC1[1] + REFREC[1];
			zrec  = XREC1[2] + REFREC[2];
			iposx1 = irnd_to_grid(xrec,xg,xpg,NXG[0]);
			iposy1 = irnd_to_grid(yrec,yg,ypg,NXG[1]);
			iposz1 = irnd_to_grid(zrec,zg,zpg,NXG[2]);

			xrec  = XREC2[0] + REFREC[0];
			yrec  = XREC2[1] + REFREC[1];
			zrec  = XREC2[2] + REFREC[2];
			iposx2 = irnd_to_grid(xrec,xg,xpg,NXG[0]);
			iposy2 = irnd_to_grid(yrec,yg,ypg,NXG[1]);
			iposz2 = irnd_to_grid(zrec,zg,zpg,NXG[2]);

			/* number of receiver positions*/
			NTR = 0;
			for (i=iposx1[2];i<=iposx2[2];i+=DXREC[0]){
				for (j=iposy1[2];j<=iposy2[2];j+=DXREC[1]){
					for (k=iposz1[2];k<=iposz2[2];k+=DXREC[2]){
						++(NTR);
					}
				}
			}

			recpos = imatrix(1,NTR,1,7);
			itr = 0;
			for (i=iposx1[2];i<=iposx2[2];i+=DXREC[0]){
				for (j=iposy1[2];j<=iposy2[2];j+=DXREC[1]){
					for (k=iposz1[2];k<=iposz2[2];k+=DXREC[2]){
						itr++;
						recpos[itr][1] = i;
						recpos[itr][2] = j;
						recpos[itr][3] = k;
						recpos[itr][4] = i;
						recpos[itr][5] = j;
						recpos[itr][6] = k;
						recpos[itr][7] = itr;
					}
				}
			}

			free_ivector(iposx1,1,2);
			free_ivector(iposy1,1,2);
			free_ivector(iposz1,1,2);
			free_ivector(iposx2,1,2);
			free_ivector(iposy2,1,2);
			free_ivector(iposz2,1,2);
		}
	}

	/* send global receiver positions to all CPUs*/
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(&NTR,1,MPI_INT,0,MPI_COMM_WORLD);
	if (MYID)
		recpos = imatrix(1,NTR,1,7);
	MPI_Bcast(&recpos[1][1],NTR*7,MPI_INT,0,MPI_COMM_WORLD);

	if (!(MYID)){
		fprintf(fp," Number of receiver positions found: %i\n",NTR);
		fprintf(fp," Receiver positions (in gridpoints) in the global model-system:\n");
		fprintf(fp," x	\ty	\tz	\txp	\typ	\tzp \n");
		for (itr=1;itr<=NTR;itr++)
			fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\n",xg[recpos[itr][1]],yg[recpos[itr][2]],zg[recpos[itr][3]],xpg[recpos[itr][4]],ypg[recpos[itr][5]],zpg[recpos[itr][6]]);
		fprintf(fp,"\n\n");
	}

	return recpos;
}
