/*------------------------------------------------------------------------
 *   write values of stress tensor field at the edges of the
 *   local grid into buffer arrays and  exchange between
 *   processes.
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


double exchange_s(float ***p, 
float ** sbuf_rig_to_lef, float ** sbuf_fro_to_bac, float ** sbuf_bot_to_top, 
float ** rbuf_rig_to_lef, float ** rbuf_fro_to_bac, float ** rbuf_bot_to_top, 
MPI_Request * request, 
int infoout){

	/* extern variables */
	extern int	NX[3];
	extern FILE	*FP;

	/* MPI variables */
	MPI_Status	status[6];

	/* local variables */
	int		i, j, k;
	double		time=0, time1=0, time2=0;


	/* timing */
	time1 = MPI_Wtime();


	/* exchange bottom-top */
	for (j=1;j<=NX[1];j++){
		for (i=1;i<=NX[0];i++){
			/* storage of bottom of local volume into buffer */
			sbuf_bot_to_top[i][j] = p[i][j][NX[2]];
		}
	}
	/* exchange front-back */
	for (k=1;k<=NX[2];k++){
		for (i=1;i<=NX[0];i++){
			/* storage of front of local volume into buffer */
			sbuf_fro_to_bac[i][k] = p[i][NX[1]][k];
		}
	}
	/* exchange right-left */
	for (k=1;k<=NX[2];k++){
		for (j=1;j<=NX[1];j++){
			/* storage of right edge of local volume into buffer */
			sbuf_rig_to_lef[j][k] = p[NX[0]][j][k];
		}
	}

	MPI_Startall(6,&request[0]);
	MPI_Waitall(6, &request[0],&status[0]);

	/* exchange bottom-top */
	for (j=1;j<=NX[1];j++){
		for (i=1;i<=NX[0];i++){
			p[i][j][0] = rbuf_bot_to_top[i][j];
		}
	}
	/* exchange front-back */
	for (k=1;k<=NX[2];k++){
		for (i=1;i<=NX[0];i++){
			p[i][0][k] = rbuf_fro_to_bac[i][k];
		}
	}
	/* exchange right-left */
	for (k=1;k<=NX[2];k++){
		for (j=1;j<=NX[1];j++){
			p[0][j][k] = rbuf_rig_to_lef[j][k];
		}
	}


	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for stress exchange: \t\t %4.2f s.\n",time);

	return time;
}
