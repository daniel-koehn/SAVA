/*------------------------------------------------------------------------
 *   write values of particle velocities at the edges of the
 *   local grid into buffer arrays and  exchange between
 *   processes.
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/

#include "fd.h"


double exchange_v(struct vector3d ***v,
float ** sbuf_lef_to_rig, float ** sbuf_bac_to_fro, float ** sbuf_top_to_bot, 
float ** rbuf_lef_to_rig, float ** rbuf_bac_to_fro, float ** rbuf_top_to_bot, 
MPI_Request * request, 
int infoout){

	/* extern variables */
	extern int	NX[3];
	extern FILE	*FP;

	/* MPI variables */
	MPI_Status	status[6];

	/* local variables */
	int		i, j, k;
	const int	nx1 = NX[0]+1, ny1 = NX[1]+1, nz1 = NX[2]+1;
	double		time=0, time1=0, time2=0;


	/* timing */
	time1 = MPI_Wtime();


	/* exchange top-bottom */
	for (i=1;i<=NX[0];i++){
		for (j=1;j<=NX[1];j++){
			/* storage of top of local volume into buffer */
			sbuf_top_to_bot[i][j] = v[i][j][1].z;
		}
	}
	/* exchange back-front */
	for (i=1;i<=NX[0];i++){
		for (k=1;k<=NX[2];k++){
			/* storage of back of local volume into buffer */
			sbuf_bac_to_fro[i][k] = v[i][1][k].y;
		}
	}
	/* exchange left-right */
	for (j=1;j<=NX[1];j++){
		for (k=1;k<=NX[2];k++){
			/* storage of left edge of local volume into buffer */
			sbuf_lef_to_rig[j][k] = v[1][j][k].x;
		}
	}

	MPI_Startall(6,&request[0]);
	MPI_Waitall(6, &request[0],&status[0]);

	/* exchange top-bottom */
	for (i=1;i<=NX[0];i++){
		for (j=1;j<=NX[1];j++){
			v[i][j][nz1].z = rbuf_top_to_bot[i][j];
		}
	}
	/* exchange back-front */
	for (i=1;i<=NX[0];i++){
		for (k=1;k<=NX[2];k++){
			v[i][ny1][k].y = rbuf_bac_to_fro[i][k];
		}
	}
	/* exchange left-right */
	for (j=1;j<=NX[1];j++){
		for (k=1;k<=NX[2];k++){
			v[nx1][j][k].x = rbuf_lef_to_rig[j][k];
		}
	}


	/* timing */
	time2 = MPI_Wtime();
	time  = time2-time1;
	if (infoout)
		fprintf(FP," Real time for particle velocity exchange: \t %4.2f s.\n",time);

	return time;
}
