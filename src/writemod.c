/*------------------------------------------------------------------------
 *   write local model to file              
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


void writemod(char modfile[STRING_SIZE], float *** array, int format){

	/* different data formats of output:
	   SNAP_FORMAT=1  :  SU (native byte order)
	   SNAP_FORMAT=2  :  ASCII
	   SNAP_FORMAT=3  :  BINARY (native byte order) */


	/* extern variables */
	extern int	MYID, NX[3], POS[3], MODEL_IDX[3];
	extern FILE	*FP;

	/* local variables */
	int	i, j, k, l;
	int	na;
	float	*array1=NULL;
	char	file[STRING_SIZE];
	FILE	*fpmod;



	fprintf(FP,"\n **Message from function writemod (printed by PE %d):\n",MYID);

	/* number of samples to write */
	na = ((int)((NX[0]-1)/MODEL_IDX[0])+1)*((int)((NX[1]-1)/MODEL_IDX[1])+1)*((int)((NX[2]-1)/MODEL_IDX[2])+1);

	array1 = vector(1,na);

	l = 0;
	for (i=1;i<=NX[0];i+=MODEL_IDX[0]){
		for (j=1;j<=NX[1];j+=MODEL_IDX[1]){
			for (k=1;k<=NX[2];k+=MODEL_IDX[2]){
				l++;
				array1[k] = array[i][j][k];
			}
		}
	}
	
	fprintf(FP," Writing model to \n");
	sprintf(file,"%s.%i.%i.%i",modfile,POS[0],POS[1],POS[2]);
	fprintf(FP,"\t%s\n", file);
	fprintf(FP," Note: the spatial sampling is the same as for the snapshot data!\n");
	fprintf(FP," x-direction (IDX): every %d grid point\n", MODEL_IDX[0]);
	fprintf(FP," y-direction (IDY): every %d grid point\n", MODEL_IDX[1]);
	fprintf(FP," z-direction (IDZ): every %d grid point\n", MODEL_IDX[2]);

	fpmod = fopen(file,"w");
	writedsk_array(fpmod,na,&array1[1],format);
	fclose(fpmod);

}


