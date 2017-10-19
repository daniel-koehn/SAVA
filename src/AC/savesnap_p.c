/*------------------------------------------------------------------------
 *   merging snapshot files
 *
 *   O. Hellwig
 *  ---------------------------------------------------------------------- */

#include "fd.h"
#include "fd3d.h"

void savesnap(FILE *fp, float * xpg, float * ypg, float * zpg){

	/* extern variables */
	extern int	MYID;
	extern int	SNAP;
	extern float	TSNAP1, TSNAP2, TSNAPINC;
	extern char	SNAP_FILE[STRING_SIZE];
	

	/* local variables */
	int	nsnap;
	char	outfiles[STRING_SIZE];

	/* log information */
	fprintf(fp,"\n **Message from function savesnap (printed by PE %d):\n",MYID);

	/* number of snapshots */
	nsnap = iround((TSNAP2-TSNAP1)/TSNAPINC)+1;


	if (SNAP & 8){
		/* particle acceleration */
		merge(nsnap,7,xpg,ypg,zpg);
		merge(nsnap,8,xpg,ypg,zpg);
		merge(nsnap,9,xpg,ypg,zpg);	
	}
	if (SNAP & 4){
		/* div */
		merge(nsnap,13,xpg,ypg,zpg);
	}
	if (SNAP & 2){
		/* pressure */
		merge(nsnap,14,xpg,ypg,zpg);
	}
	if (SNAP & 1){
		/* particle velocity */
		merge(nsnap,15,xpg,ypg,zpg);
		merge(nsnap,16,xpg,ypg,zpg);
		merge(nsnap,17,xpg,ypg,zpg);
	}
	
	/* cleanup temporary files */
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(outfiles,"rm %s*.bin.*",SNAP_FILE);
	printf("%s\n",outfiles);
	system(outfiles);

}
