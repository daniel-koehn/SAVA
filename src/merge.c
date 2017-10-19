/*------------------------------------------------------------------------
 *   merge snapshots files written by the different processes to 
 *   a single file
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void merge(int nsnap, int type, float *x, float *y, float *z){

	/* extern variables */
	extern char	SNAP_FILE[STRING_SIZE];
	extern int	SNAP_FORMAT, NPROCX[3];
	extern int	NX[3], NXG[3], IDX[3];
	extern int	FFID;
	extern FILE	*FP;
	extern float	XSNAPMIN, XSNAPMAX, YSNAPMIN, YSNAPMAX, ZSNAPMIN, ZSNAPMAX;

	/* local variables */
	char	file[STRING_SIZE], mfile[STRING_SIZE], outfile[STRING_SIZE], ext[12];
	FILE	*fp[NPROCX[0]][NPROCX[1]][NPROCX[2]], *fpout;
	int	i, j, ip, jp, kp, n;
	int	ix, iy, iz;
	int	na;
	int	ixmin, ixmax, iymin, iymax, izmin, izmax;
	int	ixsnapmin, ixsnapmax, iysnapmin, iysnapmax, izsnapmin, izsnapmax;
//	float	a, amax;


	switch(type){
	case 1: 
		fprintf(FP," xx-component of stress (txx)");
		sprintf(ext,"_txx");
		break;
	case 2: 
		fprintf(FP," xy-component of stress (txy)");
		sprintf(ext,"_txy");
		break;
	case 3: 
		fprintf(FP," xz-component of stress (txz)");
		sprintf(ext,"_txz");
		break;
	case 4: 
		fprintf(FP," yy-component of stress (tyy)");
		sprintf(ext,"_tyy");
		break;
	case 5: 
		fprintf(FP," yz-component of stress (tyz)");
		sprintf(ext,"_tyz");
		break;
	case 6: 
		fprintf(FP," zz-component of stress (tzz)");
		sprintf(ext,"_tzz");
		break;
	case 7: 
		fprintf(FP," x-component of particle acceleration (ax)");
		sprintf(ext,"_ax");
		break;
	case 8: 
		fprintf(FP," y-component of particle acceleration (ay)");
		sprintf(ext,"_ay");
		break;
	case 9: 
		fprintf(FP," z-component of particle acceleration (az)");
		sprintf(ext,"_az");
		break;
	case 10: 
		fprintf(FP," x-component of S-wavefield (curlx)");
		sprintf(ext,"_curlx");
		break;
	case 11: 
		fprintf(FP," y-component of S-wavefield (curly)");
		sprintf(ext,"_curly");
		break;
	case 12: 
		fprintf(FP," z-component of S-wavefield (curlz)");
		sprintf(ext,"_curlz");
		break; 
	case 13: 
		fprintf(FP," P-wavefield (div)");
		sprintf(ext,"_div");
		break;
	case 14: 
		fprintf(FP," pressure (p)");
		sprintf(ext,"_p");
		break;
	case 15: 
		fprintf(FP," x-component of particle velocity (vx)");
		sprintf(ext,"_vx");
		break;
	case 16: 
		fprintf(FP," y-component of particle velocity (vy)");
		sprintf(ext,"_vy");
		break;
	case 17: 
		fprintf(FP," z-component of particle velocity (vz)");
		sprintf(ext,"_vz");
		break;     
	default: 
		error(" merge: cannot find snapfiles! ");
		break;
	}

	switch(SNAP_FORMAT){
	case 1: 
		strcat(ext,".su");
		break;
	case 2: 
		strcat(ext,".asc");
		break;
	case 3: 
		strcat(ext,".bin");
		break;
	}

	sprintf(mfile,"%s%.4i%s",SNAP_FILE,FFID,ext);
	fprintf(FP," (files: %s.??? ).\n",mfile);

	sprintf(outfile,"%s%.4i%s",SNAP_FILE,FFID,ext);
	fprintf(FP,"\n writing merged snapshot file to  %s \n",outfile);
	fpout=fopen(outfile,"w");

	fprintf(FP," Opening snapshot files: %s.??? ",mfile);

	/* loop over processes */
	ixmin = irnd_to_grid_min(XSNAPMIN,x,1,NXG[0],IDX[0]);
	ixmax = irnd_to_grid_max(XSNAPMAX,x,1,NXG[0],IDX[0]);
	iymin = irnd_to_grid_min(YSNAPMIN,y,1,NXG[1],IDX[1]);
	iymax = irnd_to_grid_max(YSNAPMAX,y,1,NXG[1],IDX[1]);
	izmin = irnd_to_grid_min(ZSNAPMIN,z,1,NXG[2],IDX[2]);
	izmax = irnd_to_grid_max(ZSNAPMAX,z,1,NXG[2],IDX[2]);

	/* open snapshot files */
	for (ip=0;ip<NPROCX[0]; ip++){
		ix        = ip*NX[0];
		ixsnapmin = max(ixmin,ix+1);
		ixsnapmax = min(ixmax,ix+NX[0]);

		if (ixsnapmin <= ixsnapmax){
			for (jp=0;jp<NPROCX[1]; jp++){
				iy        = jp*NX[1];
				iysnapmin = max(iymin,iy+1);
				iysnapmax = min(iymax,iy+NX[1]);

				if (iysnapmin <= iysnapmax){
					for (kp=0;kp<NPROCX[2]; kp++){
						iz        = kp*NX[2];
						izsnapmin = max(izmin,iz+1);
						izsnapmax = min(izmax,iz+NX[2]);

						if (izsnapmin <= izsnapmax){
							sprintf(file,"%s.%i.%i.%i",mfile,ip,jp,kp);
							fp[ip][jp][kp]=fopen(file,"r");
							if (fp[ip][jp][kp]==NULL) 
								error("merge: can't read snapfile !");
						}
					}
				}
			}
		}
	}

	fprintf(FP," ... finished. \n");
	fprintf(FP," Copying...");

	/* loop over time steps */
//	amax  = 0.0;
	for (n=0;n<nsnap; n++){
		/* loop over processes */
		for (ip=0;ip<NPROCX[0]; ip++){
			ix        = ip*NX[0];
			ixsnapmin = max(ixmin,ix+1);
			ixsnapmax = min(ixmax,ix+NX[0]);
			/* loop over grid points */
			for (i=ixsnapmin;i<=ixsnapmax;i+=IDX[0]){

				for (jp=0;jp<NPROCX[1]; jp++){
					iy        = jp*NX[1];
					iysnapmin = max(iymin,iy+1);
					iysnapmax = min(iymax,iy+NX[1]);
					/* loop over grid points */
					for (j=iysnapmin;j<=iysnapmax;j+=IDX[1]){

						/* loop over processes */
						for (kp=0;kp<NPROCX[2]; kp++){
							iz        = kp*NX[2];
							izsnapmin = max(izmin,iz+1);
							izsnapmax = min(izmax,iz+NX[2]);

							if (izsnapmax>=izsnapmin){
								/* number of samples to write */
								na = (int)((izsnapmax-izsnapmin)/IDX[2])+1;

								copydsk_array(fp[ip][jp][kp],fpout,na,SNAP_FORMAT);
							}
						}
					}
				}
			}
		}
	}

	fprintf(FP," ... finished. \n");

	/* close snapshot files */
	for (ip=0;ip<NPROCX[0]; ip++){
		ix        = ip*NX[0];
		ixsnapmin = max(ixmin,ix+1);
		ixsnapmax = min(ixmax,ix+NX[0]);

		if (ixsnapmin <= ixsnapmax){
			for (jp=0;jp<NPROCX[1]; jp++){
				iy        = jp*NX[1];
				iysnapmin = max(iymin,iy+1);
				iysnapmax = min(iymax,iy+NX[1]);

				if (iysnapmin <= iysnapmax){
					for (kp=0;kp<NPROCX[2]; kp++){
						iz        = kp*NX[2];
						izsnapmin = max(izmin,iz+1);
						izsnapmax = min(izmax,iz+NX[2]);

						if (izsnapmin <= izsnapmax){
							fclose(fp[ip][jp][kp]);
						}
					}
				}
			}
		}
	}

	/*if (SNAP_FORMAT==3){
		fprintf(FP," Use \n");
		fprintf(FP," xmovie n1=%d n2=%d  < %s loop=1 label1=Y label2=X title=%%g d1=%f d2=%f f1=%f f2=%f clip=%e \n",
			(int)((izmax-izmin)/IDX[2])+1,(int)((ixmax-ixmin)/IDX[0])+1,outfile,1.0,1.0,1.0,1.0,amax/100.0);
		fprintf(FP," to play the movie. \n");
	
	}*/

	/*if (SNAP_FORMAT==3){
		fprintf(FP," Use \n");
		fprintf(FP," xmovie n1=%d n2=%d  < %s loop=1 label1=Y label2=X title=%%g d1=%f d2=%f f1=%f f2=%f clip=%e \n",
			(int)((izmax-izmin+1)/IDX[2])+1,(int)((ixmax-ixmin+1)/IDX[0])+1,outfile,1.0,1.0,1.0,1.0,amax/100.0);
		fprintf(FP," to play the movie. \n");
	
	}*/


}


