/*------------------------------------------------------------------------
 *   Computation of local source coordinates
 * 
 *   T. Bohlen, modified by O. Hellwig
 *  ----------------------------------------------------------------------*/
 
 #include "fd.h"
 
float **splitsrc(float **srcpos){

	/* extern variables */
	extern int	NX[3], POS[3], MYID, NSRC, NSRC_LOC;
	extern FILE	*FP;

	/* local variables */
	int	a, b, c;
	int	i, j, k;
	float	sum;
	float	**srcpos_dummy, **srcpos_local=NULL;


	fprintf(FP,"\n **Message from splitsrc (printed by PE %d):\n",MYID);
	srcpos_dummy = matrix(1,NSRC,1,2);

	i   = 0;
	sum = 0.0;
	for (j=1;j<=NSRC;j++) {
		a = (srcpos[j][1]-1)/NX[0];
		b = (srcpos[j][2]-1)/NX[1];
		c = (srcpos[j][3]-1)/NX[2];

		if ((POS[0]==a)&&(POS[1]==b)&&(POS[2]==c)){
			i++;
			srcpos_dummy[i][1] = (float)j;
			srcpos_dummy[i][2] = sum;
		}
		sum += srcpos[j][8];
	}
   
	if (i)
		srcpos_local = matrix(1,i,1,9);
	for (k=1;k<=i;k++){
		j   = (int)srcpos_dummy[k][1];
		sum =      srcpos_dummy[k][2];

		srcpos_local[k][1] = (float)(((int)(srcpos[j][1]-1)%NX[0])+1);
		srcpos_local[k][2] = (float)(((int)(srcpos[j][2]-1)%NX[1])+1);
		srcpos_local[k][3] = (float)(((int)(srcpos[j][3]-1)%NX[2])+1);
		srcpos_local[k][4] = srcpos[j][4];
		srcpos_local[k][5] = srcpos[j][5];
		srcpos_local[k][6] = srcpos[j][6];
		srcpos_local[k][7] = srcpos[j][7];
		srcpos_local[k][8] = srcpos[j][8];
		srcpos_local[k][9] = j;
	}

	free_matrix(srcpos_dummy,1,NSRC,1,2);

	fprintf(FP," Splitting of source positions from global to local grids finished.\n");
	fprintf(FP," Number of local source positions: %d \n",i);
	if (i){
		fprintf(FP," Table of local source positions (in gridpoints), time-shift, centre frequency and amplitude:\n");
		fprintf(FP," MYID\t  x\t  y\t  z\t  tshift\t  fc\t  amp\t  type\t  nts\n");
	}
	for (j=1;j<=i;j++)
		fprintf(FP," %3d\t%4.0f\t%4.0f\t%4.0f\t%4.0f\t%4.0f\t%4.0f\t%6.2f\t%6.2f\n",MYID,
			srcpos_local[j][1],srcpos_local[j][2],srcpos_local[j][3],srcpos_local[j][4],
			srcpos_local[j][5],srcpos_local[j][6],srcpos_local[j][7],srcpos_local[j][8]);
	fprintf(FP,"\n");

	NSRC_LOC = i;

	return srcpos_local;
}
