/*------------------------------------------------------------------------
 *   Computation of local receiver coordinates
 *   (within each subgrid)
 *
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/
 
 #include "fd.h"
int **splitrec(int *ntr_loc, int **recpos, int ntr, int *recswitch){

	/* extern variables */
	extern int	NX[3], POS[3], MYID;
	extern FILE	*FP;

	/* local variables */
	int	a, b, c;
	int	i, j, k;
	int	*recpos_dummy, **recpos_local=NULL;


	fprintf(FP,"\n **Message from splitrec (printed by PE %d):\n",MYID);
	recpos_dummy = ivector(1,ntr);

	i = 0;
	for (j=1;j<=ntr;j++){
		recswitch[j] = 0;
		a = (recpos[j][1]-1)/NX[0];
		b = (recpos[j][2]-1)/NX[1];
		c = (recpos[j][3]-1)/NX[2];

		if ((POS[0]==a)&&(POS[1]==b)&&(POS[2]==c)){
			recswitch[j] = 1;
			i++;
			recpos_dummy[i] = j;
		}
	}
   
	if (i)
		recpos_local = imatrix(1,i,1,7);
	for (k=1;k<=i;k++){
		j = recpos_dummy[k];

		recpos_local[k][1] = ((recpos[j][1]-1)%NX[0])+1;
		recpos_local[k][2] = ((recpos[j][2]-1)%NX[1])+1;
		recpos_local[k][3] = ((recpos[j][3]-1)%NX[2])+1;
		recpos_local[k][4] = ((recpos[j][4]-1)%NX[0])+1;
		recpos_local[k][5] = ((recpos[j][5]-1)%NX[1])+1;
		recpos_local[k][6] = ((recpos[j][6]-1)%NX[2])+1;
		recpos_local[k][7] =   recpos[j][7];
	}

	free_ivector(recpos_dummy,1,ntr);

	fprintf(FP," Splitting of receivers from global to local grids finished.\n");
	fprintf(FP," Number of local receiver positions: %d \n",i);
	if (i){
		fprintf(FP," Table of local receiver positions (in gridpoints):\n");
		fprintf(FP," MYID \t x \t y \t z  \t xp \t yp \t zp\n");
	}
	for (j=1;j<=i;j++)
		fprintf(FP," %d \t %d \t %d \t %d \t %d \t %d \t %d\n",MYID,recpos_local[j][1],recpos_local[j][2],recpos_local[j][3],recpos_local[j][4],recpos_local[j][5],recpos_local[j][6]);
	fprintf(FP,"\n");

	*ntr_loc = i;

	return recpos_local;

}
