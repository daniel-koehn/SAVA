/*------------------------------------------------------------------------
 *   Write 2D snapshot for current timestep  to file                                   
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void snap(FILE *fp, int nt, int nsnap, struct vector3d ***v, float ***p, struct vector3d ***a, float ***diverg, float *x, float *y, float *z){

	/* 	different data formats of output:
		SNAP_FORMAT=1  :  SU (native byte order)
		SNAP_FORMAT=2  :  ASCII
		SNAP_FORMAT=3  :  BINARY (native byte order)
		
		output is:
		SNAP=0:       no snapshots
		SNAP=SNAP+1:  add particle-velocities
		SNAP=SNAP+2:  add pressure (hydrophones)
		SNAP=SNAP+4:  add curl and div
		SNAP=SNAP+8:  add particle-accelerations
		SNAP=SNAP+16: add stress components*/

	/* local variables */
	int	i, j, k, l;
	int	na;
	float	*array1=NULL, *array2=NULL, *array3=NULL;
	char	snapfile1[STRING_SIZE], ext[8], wm[2];
	FILE	*fp1;

	/* extern variables */
	extern float	DT;
	extern char	SNAP_FILE[STRING_SIZE];
	extern int	SNAP_FORMAT, SNAP;
	extern int	IXSNAPMIN[3], IXSNAPMAX[3];
	extern int	MYID, POS[3], IDX[3];
	extern int	FFID;


	switch(SNAP_FORMAT){
	case 1:
		sprintf(ext,".su");
		break;
	case 2:
		sprintf(ext,".asc");
		break;
	case 3:
		sprintf(ext,".bin");
		break;
	}

	if (nsnap==1) 
		sprintf(wm,"w");
	else 
		sprintf(wm,"a");

	fprintf(fp,"\n\n PE %d is writing snapshot-data at T=%fs to \n",MYID,nt*DT);

	/* number of samples to write */
	na = ((int)((IXSNAPMAX[0]-IXSNAPMIN[0])/IDX[0])+1)*((int)((IXSNAPMAX[1]-IXSNAPMIN[1])/IDX[1])+1)*((int)((IXSNAPMAX[2]-IXSNAPMIN[2])/IDX[2])+1);

	if ((IXSNAPMAX[0]>=IXSNAPMIN[0]) && (IXSNAPMAX[1]>=IXSNAPMIN[1]) && (IXSNAPMAX[2]>=IXSNAPMIN[2])){

		if (SNAP & 8){
			/* particle acceleration */
			array1 = vector(1,na);
			array2 = vector(1,na);
			array3 = vector(1,na);

			l=0;
			for (i=IXSNAPMIN[0];i<=IXSNAPMAX[0];i+=IDX[0]){
				for (j=IXSNAPMIN[1];j<=IXSNAPMAX[1];j+=IDX[1]){
					for (k=IXSNAPMIN[2];k<=IXSNAPMAX[2];k+=IDX[2]){
						l++;
						array1[l] = a[i][j][k].x;
						array2[l] = a[i][j][k].y;
						array3[l] = a[i][j][k].z;
					}
				}
			}

			sprintf(snapfile1,"%s%.4i_ax%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array1[1],SNAP_FORMAT);
			fclose(fp1);


			sprintf(snapfile1,"%s%.4i_ay%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array2[1],SNAP_FORMAT);
			fclose(fp1);


			sprintf(snapfile1,"%s%.4i_az%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array3[1],SNAP_FORMAT);
			fclose(fp1);


			free_vector(array1,1,na);
			free_vector(array2,1,na);
			free_vector(array3,1,na);
		}
		if (SNAP & 4){
			/* div */
			array1 = vector(1,na);

			l=0;
			for (i=IXSNAPMIN[0];i<=IXSNAPMAX[0];i+=IDX[0]){
				for (j=IXSNAPMIN[1];j<=IXSNAPMAX[1];j+=IDX[1]){
					for (k=IXSNAPMIN[2];k<=IXSNAPMAX[2];k+=IDX[2]){
						l++;
						array1[l] = diverg[i][j][k];
					}
				}
			}

			sprintf(snapfile1,"%s%.4i_div%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array1[1],SNAP_FORMAT);
			fclose(fp1);

			free_vector(array1,1,na);
		}
		if (SNAP & 2){
			/* pressure */
			array1 = vector(1,na);

			l=0;
			for (i=IXSNAPMIN[0];i<=IXSNAPMAX[0];i+=IDX[0]){
				for (j=IXSNAPMIN[1];j<=IXSNAPMAX[1];j+=IDX[1]){
					for (k=IXSNAPMIN[2];k<=IXSNAPMAX[2];k+=IDX[2]){
						l++;
						array1[l] = p[i][j][k];
					}
				}
			}

			sprintf(snapfile1,"%s%.4i_p%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array1[1],SNAP_FORMAT);
			fclose(fp1);

			free_vector(array1,1,na);
		}
		if (SNAP & 1){
			/* particle velocity */
			array1 = vector(1,na);
			array2 = vector(1,na);
			array3 = vector(1,na);

			l=0;
			for (i=IXSNAPMIN[0];i<=IXSNAPMAX[0];i+=IDX[0]){
				for (j=IXSNAPMIN[1];j<=IXSNAPMAX[1];j+=IDX[1]){
					for (k=IXSNAPMIN[2];k<=IXSNAPMAX[2];k+=IDX[2]){
						l++;
						array1[l] = v[i][j][k].x;
						array2[l] = v[i][j][k].y;
						array3[l] = v[i][j][k].z;
					}
				}
			}

			sprintf(snapfile1,"%s%.4i_vx%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array1[1],SNAP_FORMAT);
			fclose(fp1);


			sprintf(snapfile1,"%s%.4i_vy%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array2[1],SNAP_FORMAT);
			fclose(fp1);


			sprintf(snapfile1,"%s%.4i_vz%s.%i.%i.%i",SNAP_FILE,FFID,ext,POS[0],POS[1],POS[2]);
			fprintf(fp,"%s\n\n",snapfile1);

			fp1 = fopen(snapfile1,wm);
			writedsk_array(fp1,na,&array3[1],SNAP_FORMAT);
			fclose(fp1);


			free_vector(array1,1,na);
			free_vector(array2,1,na);
			free_vector(array3,1,na);
		}
	}
}


