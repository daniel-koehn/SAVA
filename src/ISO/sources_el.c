/*------------------------------------------------------------------------ 
 *   Reading (distributed) source positions, timeshift, centre frequency 
 *   and amplitude from SRC_FILE.
 *
 *   O. Hellwig, S. Wenk 
 *  ----------------------------------------------------------------------*/

#include "fd.h"

float **sources(float *xg, float *yg, float *zg, float *xpg, float *ypg, float *zpg){

	/* extern variables */
	extern int	NXG[3];
	extern float	TS;
	extern float	REFSRC[3];
	extern char	SRC_FILE[STRING_SIZE];
	extern int	MYID, NSRC;
	extern FILE	*FP;

	/* local variables */
	int	l;
	int	*iposx, *iposy, *iposz;
	float	**srcpos=NULL;
	float	xsrc, ysrc, zsrc, tshift, fc, amp, type, nsamp;
	FILE	*fpsrc;


	if (!(MYID)){	/* read source positions from file */
		fprintf(FP,"\n **Message from function sources (printed by PE %d):\n",MYID);
		fprintf(FP," Reading source positions, time-shift, centre frequency \n");
		fprintf(FP," and amplitude from file: %s\n",SRC_FILE);
		fpsrc = fopen(SRC_FILE,"r");
	
		if (fpsrc==NULL) 
			error(" Source file could not be opened !");
		NSRC=0;

		/* read number of source positions */
		fscanf(fpsrc,"%d",&NSRC);
		fprintf(FP," Number of source positions found: %d \n",NSRC);

		iposx  = ivector(1,2);
		iposy  = ivector(1,2);
		iposz  = ivector(1,2);
		srcpos = matrix(1,NSRC,1,8);
		for (l=1;l<=NSRC;l++){
			fscanf(fpsrc,"%f%f%f%f%f%f%f%f",&xsrc, &ysrc, &zsrc, &tshift, &fc, &amp, &type, &nsamp);

			xsrc += REFSRC[0];
			ysrc += REFSRC[1];
			zsrc += REFSRC[2];
			iposx = irnd_to_grid(xsrc,xg,xpg,NXG[0]);
			iposy = irnd_to_grid(ysrc,yg,ypg,NXG[1]);
			iposz = irnd_to_grid(zsrc,zg,zpg,NXG[2]);
			switch ((int)type){
				case 1: /* explosion source (xp,yp,zp) */
					srcpos[l][1] = (float)iposx[2];
					srcpos[l][2] = (float)iposy[2];
					srcpos[l][3] = (float)iposz[2];
					break;
				case 2: /* acceleration in x-direction (x,yp,zp) */
					srcpos[l][1] = (float)iposx[1];
					srcpos[l][2] = (float)iposy[2];
					srcpos[l][3] = (float)iposz[2];
					break;
				case 3: /* acceleration in y-direction (xp,y,zp) */
					srcpos[l][1] = (float)iposx[2];
					srcpos[l][2] = (float)iposy[1];
					srcpos[l][3] = (float)iposz[2];
					break;
				case 4: /* acceleration in z-direction (xp,yp,z) */
					srcpos[l][1] = (float)iposx[2];
					srcpos[l][2] = (float)iposy[2];
					srcpos[l][3] = (float)iposz[1];
					break;
				case 5: /* force in x-direction (x,yp,zp) */
					srcpos[l][1] = (float)iposx[1];
					srcpos[l][2] = (float)iposy[2];
					srcpos[l][3] = (float)iposz[2];
					break;
				case 6: /* force in y-direction (xp,y,zp) */
					srcpos[l][1] = (float)iposx[2];
					srcpos[l][2] = (float)iposy[1];
					srcpos[l][3] = (float)iposz[2];
					break;
				case 7: /* force in z-direction (xp,yp,z) */
					srcpos[l][1] = (float)iposx[2];
					srcpos[l][2] = (float)iposy[2];
					srcpos[l][3] = (float)iposz[1];
					break;
				case 8:	/* double couple in x-y-plane (x,y,zp) */
					srcpos[l][1] = (float)iposx[1];
					srcpos[l][2] = (float)iposy[1];
					srcpos[l][3] = (float)iposz[2];
					break;
				case 9:	/* double couple in y-z-plane (xp,y,z) */
					srcpos[l][1] = (float)iposx[2];
					srcpos[l][2] = (float)iposy[1];
					srcpos[l][3] = (float)iposz[1];
					break;
				case 10:/* double couple in x-z-plane (x,yp,z) */
					srcpos[l][1] = (float)iposx[1];
					srcpos[l][2] = (float)iposy[2];
					srcpos[l][3] = (float)iposz[1];
					break;
				default: 
					error("Wrong source type specified in source file!");
			}
			srcpos[l][4] = tshift; /* time shift */
			srcpos[l][5] = fc;     /* source center frequency */
			srcpos[l][6] = amp;    /* amplitude */
			srcpos[l][7] = type;   /* source type */
			srcpos[l][8] = nsamp;  /* number of samples in source signal file */
		}
		free_ivector(iposx,1,2);
		free_ivector(iposy,1,2);
		free_ivector(iposz,1,2);

		fclose(fpsrc);

		/* Compute maximum frequency */
		for (l=1;l<=NSRC;l++)
			if (srcpos[l][5]>fc)
				fc = srcpos[l][5];
		fprintf(FP," Maximum frequency defined in %s: %6.2e Hz\n",SRC_FILE,fc);
		TS = 1.0/fc;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&NSRC,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&TS,1,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (MYID)
		srcpos = matrix(1,NSRC,1,8);
	MPI_Bcast(&srcpos[1][1],NSRC*8,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (!(MYID)){
		fprintf(FP," x \t\t y \t\t z \t\t tshift \t\t fc \t\t amp \t\t type \t\t nsamp\n");
		for (l=1;l<=NSRC;l++){
			switch ((int)type){
				case 1: /* explosion source (xp,yp,zp) */
					xsrc = xpg[(int)srcpos[l][1]];
					ysrc = ypg[(int)srcpos[l][2]];
					zsrc = zpg[(int)srcpos[l][3]];
					break;
				case 2: /* acceleration in x-direction (x,yp,zp) */
					xsrc =  xg[(int)srcpos[l][1]];
					ysrc = ypg[(int)srcpos[l][2]];
					zsrc = zpg[(int)srcpos[l][3]];
					break;
				case 3: /* acceleration in y-direction (xp,y,zp) */
					xsrc = xpg[(int)srcpos[l][1]];
					ysrc =  yg[(int)srcpos[l][2]];
					zsrc = zpg[(int)srcpos[l][3]];
					break;
				case 4: /* acceleration in z-direction (xp,yp,z) */
					xsrc = xpg[(int)srcpos[l][1]];
					ysrc = ypg[(int)srcpos[l][2]];
					zsrc =  zg[(int)srcpos[l][3]];
					break;
				case 5: /* force in x-direction (x,yp,zp) */
					xsrc =  xg[(int)srcpos[l][1]];
					ysrc = ypg[(int)srcpos[l][2]];
					zsrc = zpg[(int)srcpos[l][3]];
					break;
				case 6: /* force in y-direction (xp,y,zp) */
					xsrc = xpg[(int)srcpos[l][1]];
					ysrc =  yg[(int)srcpos[l][2]];
					zsrc = zpg[(int)srcpos[l][3]];
					break;
				case 7: /* force in z-direction (xp,yp,z) */
					xsrc = xpg[(int)srcpos[l][1]];
					ysrc = ypg[(int)srcpos[l][2]];
					zsrc =  zg[(int)srcpos[l][3]];
					break;
				case 8: /* double couple in x-y-plane (x,y,zp) */
					xsrc =  xg[(int)srcpos[l][1]];
					ysrc =  yg[(int)srcpos[l][2]];
					zsrc = zpg[(int)srcpos[l][3]];
					break;
				case 9: /* double couple in y-z-plane (xp,y,z) */
					xsrc = xpg[(int)srcpos[l][1]];
					ysrc =  yg[(int)srcpos[l][2]];
					zsrc =  zg[(int)srcpos[l][3]];
					break;
				case 10: /* double couple in x-z-plane (x,yp,z) */
					xsrc =  xg[(int)srcpos[l][1]];
					ysrc = ypg[(int)srcpos[l][2]];
					zsrc =  zg[(int)srcpos[l][3]];
					break;
			}
			fprintf(FP," %6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%6.2e\t%1.2e\t%6.2e\n",
				xsrc,ysrc,zsrc,srcpos[l][4],srcpos[l][5],srcpos[l][6],srcpos[l][7],srcpos[l][8]);
		}
		fprintf(FP,"\n\n");
	}

	return srcpos;
}
